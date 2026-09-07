# 6 · Linearization — operator-level chain rule

> **For:** retrieval / inversion developers; anyone who needs `∂R/∂x` for a parameter `x`. The runnable workflow is at [User Guide → Compute Jacobians](../jacobians.md); this page explains the *why*.
>
> **Prev:** [5 · Surfaces](05_surfaces.md) · **Next:** [7 · Architecture-Agnostic Code](07_architecture.md)

The MOM solver in [Concepts/04](04_mom_solver.md) computes
``\mathbf{R}, \mathbf{T}, \mathbf{J}`` for a given atmospheric state.
*Retrievals* need ``\partial \mathbf{R}/\partial \mathbf{x}`` for every
parameter ``\mathbf{x}`` they're trying to estimate — aerosol optical
depths, refractive indices, size distributions, gas VMRs, surface BRDF
parameters. This page explains how vSmartMOM computes those derivatives:
operator-level analytic chain rule on the adding-doubling formulas, with
ForwardDiff used *only* upstream at the optical-property boundary.

## The three-tier Jacobian

```
   State parameters x:                ForwardDiff
   τ_ref, n_r, n_i, μ_logr, σ_logr,    ──────(upstream)──────►
   VMR, BRDF, ...

           CoreScatteringOpticalPropertiesLin per layer
                  (τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺)              ── AD boundary ──

                  ┌────────────────────────────────────────────┐
                  │ Hand-coded chain rule (analytic RT kernel) │
                  │ Sanghavi 2014 App. C, Eqs C.8-C.21         │
                  └────────────────────────────────────────────┘
                                    │
                                    ▼
                              ∂R/∂x ,  ∂T/∂x
```

The Jacobian is split into two zones with a clean boundary between them.

**Upstream — the AD zone** (Mie cross-sections, gas absorption, atmospheric
profile, surface BRDF parameters). Here ForwardDiff `Dual` numbers carry
parameter derivatives forward, *or* analytic derivatives are computed by
hand for the cases where Mie series are AD-hostile (refractive index, size
distribution). The output of this zone is one
`CoreScatteringOpticalPropertiesLin` per atmospheric layer:

```julia
# src/CoreRT/types_lin.jl:119–149
struct CoreScatteringOpticalPropertiesLin{T1,T2,T3} <: AbstractOpticalPropertiesLin
    τ̇::T1        # ∂τ/∂x
    ϖ̇::T2        # ∂ϖ_0/∂x
    Ż⁺⁺::T3      # ∂Z⁺⁺/∂x
    Ż⁻⁺::T3      # ∂Z⁻⁺/∂x
end
```

The aliases `OpticalPropertyJacobian = CoreScatteringOpticalPropertiesLin`
make this the explicit AD-boundary handoff struct.

Retrieval-specific Jacobian plans may project this handoff onto a compact
active parameter basis before it enters the RT kernels. A plan stores named
physical parameter keys, one local-to-global map per spectral band, and the
native optical-property columns needed by that band. The downstream analytic
formulas are unchanged; only the trailing tangent dimension becomes smaller.
See [Compute Jacobians → Retrieval-Selected Jacobians](../jacobians.md#retrieval-selected-jacobians).

Selection is performed during the moment-invariant aerosol/gas cache build,
not after a full mixed tangent has been formed. Each native column is first
assigned to pressure, a component-local aerosol block, or the gas block.
Only those compact blocks enter optical assembly. Elemental/doubling
propagation may use a still smaller local optical basis; adding and output
accumulation use the band-local retrieval count. Fixed forward physics is retained even when its
tangent is omitted.

**Downstream — the pure-`FT` zone** (elemental, doubling, interaction).
Here every kernel is plain `Float32` or `Float64` and every derivative is
computed by a hand-coded tangent-linear partner of the forward kernel.
``\mathbf{R}, \mathbf{T}, \mathbf{J}`` and their derivatives
``\dot{\mathbf{R}}, \dot{\mathbf{T}}, \dot{\mathbf{J}}`` propagate together
through the same adding-doubling sequence.

For an external solar direction, the tangent-linear solver also carries
rectangular solar-column derivatives
``\dot R_0^{-+},\dot R_0^{+-},\dot T_0^{++},\dot T_0^{--}`` with layout
`(NquadN,nStokes,nSpec,nParams)`. They contain the local ``\dot\tau``,
``\dot\varpi`` and ``\dot Z_0`` contributions. The derivative of direct-beam
extinction above the layer is added only when the columns are contracted:

```math
\dot J_0 = \dot O_0 F_0e^{-\tau_a/\mu_0}
          + O_0F_0e^{-\tau_a/\mu_0}
            \left(-\frac{\dot\tau_a}{\mu_0}\right),
```

where ``O_0`` denotes the applicable ``R_0`` or ``T_0`` column operator.
Keeping these terms separate prevents the solar direction from re-entering
the diffuse angular operator during linearization.

The split has three benefits:

1. **The hot loop stays pure-`FT`.** No `Dual{T,V,N}` arithmetic in the
   inner kernels — those would multiply work by `1+N`.
2. **Analytic derivatives are numerically stable** through batched matrix
   inversion. ForwardDiff through `batch_inv!` works (and is supported on
   GPU), but the analytic chain rule on `(E − R·R)⁻¹` is closed-form
   ``\partial A^{-1} = -A^{-1}\,\partial A\,A^{-1}`` and avoids accumulating
   AD round-off.
3. **The chain rule on adding-doubling is closed-form.** Sanghavi 2014
   App. C derives compact expressions for the tangent-linear adding/doubling
   updates. They're roughly the same shape as the forward updates — same
   matrix structure, computed against the same operands.

## Why this is fast: the matrix inversion is reused

The forward solver computes the geometric-series inverse
``\mathbf{G}=(\mathbf{E}-\mathbf{R}\mathbf{R})^{-1}`` at each doubling
step. Its tangent follows directly from the inverse rule:

```math
\dot{\mathbf{G}} = \mathbf{G}\,\dot{(\mathbf{R}\mathbf{R})}\,\mathbf{G}.
```

The positive sign comes from differentiating ``\mathbf{E}-\mathbf{R}\mathbf{R}``.
The same forward inverse is reused for every supplied tangent direction.
Doubling uses local optical directions when they are fewer than retrieval
columns. Default matrix adding still propagates matrix tangents for each
retrieval column; the opt-in source path below propagates vectors. The total cost depends on the number of
requested columns, operator size, spectral batch, precision and backend;
there is no universal ratio to the forward-only runtime.

Keep the **core-optics solve** separate from the **physical-parameter solve**
when evaluating the development target of combined forward + Jacobians below
2× forward. The former starts with supplied ``(τ,ϖ,Z^{++},Z^{-+})`` and tangent
directions, excluding Mie, optical mixing and upstream derivative construction.
It still includes every elemental, doubling and adding step. A phase matrix is
not a scalar parameter: specify the phase perturbation basis and layer count.
Three stored elemental partial arrays do not represent the dense phase-matrix
Jacobian after multiple scattering. The dedicated core-properties benchmark in
`docs/dev_notes/jacobian_batched/` records its derivative basis explicitly.

The default `jacobian_basis=:auto` chooses a local optical basis for a
single-band solve whenever it reduces the atmospheric tangent dimension.
Use `:physical` to retain direct physical-column propagation for comparison,
or `:local` to force the factored path. Multi-band concatenated solves retain
the physical path; `:local` rejects them explicitly.

For solar illumination over Lambertian surfaces, the opt-in
`jacobian_adding=:source` path converts local doubled matrix tangents into
source perturbations at the full-column incident fields. It contracts those
vectors to retrieval columns and propagates them through fixed forward
operators. This removes retrieval-sized matrix tangents from atmospheric
adding and surface closure; source vectors and outputs still grow with the
number of requested columns. Use it with `jacobian_basis=:local`. It supports
endpoint observers only and retains external solar's TOA-only contract.
Prescribed or retrievable `SurfaceSIF` can accompany the solar beam. SIF
amplitude and slope derivatives enter as boundary source vectors; they add no
matrix directions to atmospheric doubling. The incident fields include emitted
light, so atmospheric and surface-albedo Jacobians retain its multiple scattering.
Only reflected sunlight receives the direct-beam attenuation derivative at the
surface; emitted SIF is transported from that boundary through the atmosphere.
The default `jacobian_adding=:matrix` keeps the established general path.
The new path caches local layer tangents for a backward illumination pass,
so memory still grows with layer count and spectral batch size.


Sanghavi et al. (2014), (C.25)–(C.26), factor microphysical derivatives through
``(τ_i,ω_i,Z_i,f_i)``; their truncation symbol ``β_i`` is the code's ``f^t``,
not a Greek expansion coefficient. Both the scalar truncation chain and the
normalized truncated-Greek tangents include ``df^t`` upstream.

For the current fixed Rayleigh phase, write scattering fractions
``α_i=S_i/S``, where ``S_i=τ_iϖ_i`` and ``S`` includes Rayleigh scattering.
Differentiating the mixture definitions (C.22)–(C.24) gives

```math
Z=Z_r+\sum_i α_i(Z_i-Z_r),\qquad
\dot Z=\sum_i[\dot α_i(Z_i-Z_r)+α_i\dot Z_i].
```

Thus at most `2 + 5Naer` directions suffice: optical depth, albedo, and for each
aerosol its phase difference from Rayleigh plus four truncated Mie phase
tangents. A retrieval plan removes unselected Mie directions before workspace
allocation, leaving `2 + Naer + N_selected_microphysics`. Three aerosols with
fixed microphysics need only five local directions. Their forward phases and
mixture-weight derivatives are retained. All layers share these phase
directions at each Fourier order.
Small scalar coefficient arrays, prepared once before the Fourier loop,
map the local directions to retrieval columns. Gas and profile parameters
introduce no additional phase directions under this fixed phase model.

The elemental and doubling kernels propagate these **complete matrix
directions**. After doubling, a scalar-weighted contraction produces the
physical operator tangents for adding. This is valid by linearity of the
tangent map at fixed forward state; it does not apply an elementwise phase
partial after multiple scattering. Above-layer solar attenuation is held
fixed during local propagation and appended once as
``-J\dot τ_{above}/μ_0`` after contraction.

An adjoint is therefore an option, not a prerequisite for sub-2× performance.
The scalar [Sanghavi et al. (2013)](https://doi.org/10.1016/j.jqsrt.2012.10.021)
§4.6 measured 1.29–1.98× for 14–16 parameters with tangent linearization;
the vector [Sanghavi et al. (2014)](https://doi.org/10.1016/j.jqsrt.2013.09.004)
§5.3 measured about 5× for ten parameters. These are implementation-specific
measurements, not a state-count-independent complexity bound.

The batched propagation path evaluates product-rule terms across wavelength
and parameter together. It fuses ``\dot{A}B+A\dot{B}`` on supported small GPU
operators and reuses scratch arrays through doubling and general layer
interaction. Embedded-solar Lambertian source tangents also batch all
wavelengths and parameters into matrix products. Larger GPU operators retain the reference BLAS path until their
performance is validated. CPU propagation uses in-place BLAS products.
Selecting only the retrieval's requested Jacobian columns avoids propagating
unneeded derivatives through either path.

| Approach | Work required |
|---|---|
| Analytic operator derivatives | Reused forward inverses plus products for the active parameter columns |
| ForwardDiff through the RT kernel | Values and derivative slabs flow through kernel arithmetic; cost depends on chunk size and backend |
| Forward finite differences | One baseline solve plus one perturbed solve per parameter |

Benchmark the actual retrieval workload, including model construction when it
changes between iterations. Matrix-inversion reuse is an algorithmic property;
end-to-end speedup is a measured result. Refining the upstream AD boundary can
improve construction cost independently of the RT propagation kernels.

## The chain rule on adding-doubling

The full derivation is Sanghavi 2014 App. C. The structure (skipping arithmetic):

| Eq. | What it says | Source file |
|---|---|---|
| (C.5)–(C.7) | Differentiation rules for matrix products and inverses | (foundation; used everywhere below) |
| (C.8)–(C.10) | Infinitesimal elemental/thermal derivatives; current finite-δ solar kernels differentiate SF2023-II (10)–(11) using the same calculus | `elemental_fused_lin.jl` |
| (C.11)–(C.16) | Doubling/adding derivatives — same shape as the forward Eqs (23)–(28), tangent-linear | `doubling_lin.jl`, `interaction_lin.jl` |
| (C.17)–(C.20) | D-matrix symmetry on derivatives — halves the linearized doubling cost | `doubling_lin.jl` |
| (C.21) | Final assembled derivative form — written directly by `get_elem_rt_fused!` / `get_elem_rt_SFI_fused!` into the `ap_*` arrays during the elemental step | `elemental_fused_lin.jl` |
| (C.22)–(C.24) | Layer-averaged ``\bar{\tau}``, ``\bar{\varpi}_0``, ``\bar{\mathbf{Z}}`` definitions | `compEffectiveLayerProperties.jl` (forward); `compEffectiveLayerProperties_lin.jl` (lin) |
| (C.25)–(C.26) | Chain rule from the elemental SS variables to the microphysical parameters ``(n_r, n_i, μ_{logr}, σ_{logr})`` | `elemental_fused_lin.jl` + `Scattering/types_lin.jl` |
| (C.27)–(C.31) | δ derivatives (elemental thickness from `N_doubl`) | `compEffectiveLayerProperties_lin.jl` |
| (C.32)–(C.39) | ``\bar{\varpi}_0`` and ``\bar{\mathbf{Z}}`` derivatives (post-truncation) | same |
| (C.40) | ``\dot{\mathbf{Z}}_m`` from generalized spherical harmonics | `compute_Z_matrices_lin.jl` |
| (C.41)–(C.42) | ``\dot{\beta}^*`` and ``\dot{\mathbf{B}}_l^*`` for the truncated case | `delta_m_truncation_lin.jl` |

The `_lin.jl` files in `src/CoreRT/CoreKernel/` are tangent-linear partners
of `elemental.jl`, `doubling.jl`, `interaction.jl`. Each forward kernel has
a `_lin` companion that computes the derivative quantities alongside the
forward ones in a single sweep.

## ParameterLayout — naming the Jacobian columns

The Jacobian matrix has rows indexed by `(VZA, Stokes, λ)` (the radiance
output) and columns indexed by *retrieval parameters*. The column ordering
is centralized in [`src/CoreRT/parameter_layout.jl:1–67`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/parameter_layout.jl#L1-L67):

```julia
struct ParameterLayout
    n_atmosphere::Int      # = 1   (p_surf)
    aerosol_params::Int     # = 7   (τ_ref, n_r, n_i, μ_logr, σ_logr, profile location/width)
    n_aerosols::Int
    n_gases::Int
    n_surface::Int
    n_sif::Int              # = 2 for [SIF755, slope], otherwise 0
    n_canopy::Int
end

n_total(layout)             # total number of Jacobian columns
aerosol_range(layout, iaer) # column indices for aerosol iaer (length 7)
gas_range(layout)           # column indices for gas VMRs
gas_profile_range(layout, igas, Nz) # all layers for one gas
gas_layer_index(layout, igas, iz, Nz) # one gas/layer column
surface_range(layout)       # column indices for surface params
sif_range(layout)           # [SIF755, slope], referenced to 755 nm
```

Always use these accessors instead of hand-writing arithmetic like
`1 + 7*NAer + NGasSpecies*Nz + NSurf`. Column 1 is `p_surf`; the seven aerosol
parameters per mode follow:

| # | Parameter | Meaning |
|---|---|---|
| 1 | `τ_ref` | aerosol optical depth at reference wavelength |
| 2 | `n_r`   | real refractive index |
| 3 | `n_i`   | imaginary refractive index |
| 4 | `μ_logr` | log of median radius (`LogNormal.μ`) |
| 5 | `σ_logr` | log of geometric standard deviation (`LogNormal.σ`) |
| 6 | `p₀` / `z₀` | profile location for pressure-form `Normal` / altitude-form `LogNormal` |
| 7 | `σ_p` / `σ₀` | corresponding profile width |

A run with two aerosol modes, four gases, `Nz` layers, and one surface
parameter has `1 + 2·7 + 4·Nz + 1` Jacobian columns. Gas columns use
TOA-to-BOA layer order within each species.

## What goes through ForwardDiff vs analytic

Strategy by parameter type (as currently implemented and recommended for
new parameters):

| Parameter | Recommended path | Reason |
|---|---|---|
| Lambertian albedo | analytic | trivial: ``\partial r/\partial \rho = 1/\pi`` |
| RPV / Ross-Li / Cox-Munk BRDF | ForwardDiff | low-dimensional, simple surface code |
| `τ_ref` (aerosol OD) | analytic | trivial: ``\partial \tau/\partial \tau_\mathrm{ref} = \tau/\tau_\mathrm{ref}`` |
| profile location/width (`p₀, σ_p` or `z₀, σ₀`) | analytic | already in `atmo_prof_lin.jl` |
| `n_r, n_i` (refractive index) | analytic Mie | Mie series is AD-hostile (recurrences) |
| `μ_logr, σ_logr` (size distribution) | analytic Mie | same |
| Gas VMR scaling | analytic | ``\partial \tau_\mathrm{abs}/\partial \mathrm{VMR} = \sigma`` |
| Surface pressure | ForwardDiff | affects many code paths (Rayleigh, profile, absorption) |
| Temperature profile | ForwardDiff (future) | affects absorption cross-sections nonlinearly |

For parameters in the "ForwardDiff" rows, the *upstream* code carries
`Dual{T,V,N}` numbers; the resulting `CoreScatteringOpticalPropertiesLin`
arrays are real-valued (the Dual partials having been extracted at the
boundary) and the analytic RT chain rule takes over. This is the
**hybrid AD** pattern from [Concepts/07](07_architecture.md), and it's why
ForwardDiff Duals need to flow through `batched_mul` and `batch_inv!` on
GPU — both for the upstream Mie/profile/surface code (when the user picks
ForwardDiff for those) and for forward-mode AD of higher-level functions
involving `rt_run`.

## Code anchors

| Concept | Source |
|---|---|
| Linearized RT entry | [`src/CoreRT/rt_run_lin.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/rt_run_lin.jl) |
| Linearized model construction | [`src/CoreRT/tools/lin_model_from_parameters.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/tools/lin_model_from_parameters.jl) |
| Elemental derivatives | [`src/CoreRT/CoreKernel/elemental_lin.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/elemental_lin.jl) |
| Doubling derivatives | [`src/CoreRT/CoreKernel/doubling_lin.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/doubling_lin.jl) |
| Interaction derivatives | [`src/CoreRT/CoreKernel/interaction_lin.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/interaction_lin.jl) |
| Elemental tangent chain rule | `src/CoreRT/CoreKernel/elemental_fused_lin.jl` |
| Local-basis contraction after doubling | `src/CoreRT/CoreKernel/local_jacobian.jl` |
| Local optical coefficients and shared phase basis | `src/CoreRT/LayerOpticalProperties/local_jacobian_cache.jl`, `local_jacobian.jl` |
| Three-core-variable lin type | [`src/CoreRT/types_lin.jl:119–149`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/types_lin.jl#L119-L149) |
| Optical-property Jacobian boundary | `src/CoreRT/types.jl::CoreScatteringOpticalPropertiesLin` |
| ParameterLayout | [`src/CoreRT/parameter_layout.jl:1–67`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/parameter_layout.jl#L1-L67) |
| Mie linearization (analytic) | [`src/Scattering/types_lin.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/types_lin.jl) |
| Linearized atmospheric profile | [`src/CoreRT/tools/atmo_prof_lin.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/tools/atmo_prof_lin.jl) |
| Linearized δ-M truncation | [`src/CoreRT/LayerOpticalProperties/delta_m_truncation_lin.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/LayerOpticalProperties/delta_m_truncation_lin.jl) |
| Linearized Cox-Munk | [`src/CoreRT/Surfaces/coxmunk_surface_lin.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/Surfaces/coxmunk_surface_lin.jl) |

See also [User Guide → Compute Jacobians](../jacobians.md) for the runnable
workflow, [Tutorial_Jacobians](../tutorials/Tutorial_Jacobians.md) for a
hands-on example with finite-difference cross-checks, and
[Tutorial_HybridAD](../tutorials/Tutorial_HybridAD.md) for ForwardDiff
across the analytic RT kernel.

## References

- **Sanghavi et al. (2014)**, JQSRT **133**:412–433, [doi:10.1016/j.jqsrt.2013.09.004](https://doi.org/10.1016/j.jqsrt.2013.09.004), App. C. **Primary linearization reference.**
- **Sanghavi et al. (2013)**, *Linearization of a scalar matrix operator method radiative transfer model*, JQSRT **116**:1–16, [doi:10.1016/j.jqsrt.2012.10.021](https://doi.org/10.1016/j.jqsrt.2012.10.021). (Scalar predecessor; same structure, simpler derivation.)
- Hasekamp & Landgraf (2005), *Linearization of vector radiative transfer with respect to aerosol properties and its use in satellite remote sensing*, JGR **110**:D04203. (Independent linearized vector RT for comparison.)
- Crib sheet: `docs/dev_notes/theory_references.md` §H.
