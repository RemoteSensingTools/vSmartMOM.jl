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
   τ_ref, n_r, n_i, r_m, σ_g,    ──────(upstream)──────►
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

**Downstream — the pure-`FT` zone** (elemental, doubling, interaction).
Here every kernel is plain `Float32` or `Float64` and every derivative is
computed by a hand-coded tangent-linear partner of the forward kernel.
``\mathbf{R}, \mathbf{T}, \mathbf{J}`` and their derivatives
``\dot{\mathbf{R}}, \dot{\mathbf{T}}, \dot{\mathbf{J}}`` propagate together
through the same adding-doubling sequence.

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

> **A combined forward + linearized run costs *less than 2×* a forward-only
> run.** Not `1 + N_params×forward`. Not even `2× forward`. **Less than
> 2×.** This is the production-impact property of operator-level analytic
> linearization in MOM, and it's the headline reason vSmartMOM is usable
> inside Levenberg-Marquardt loops with hundreds of state vector elements.

The dominant cost in the forward solver is the batched matrix inversion at
each doubling step (the geometric-series factor
``\mathbf{G} = (\mathbf{E} - \mathbf{R}\mathbf{R})^{-1}``) and at each
inhomogeneous-adding step (the two `T*_inv` factors in
[Concepts/04 § Adding](04_mom_solver.md#adding--interaction)).
``N_\mathrm{quad}^3`` per spectral point per layer per Fourier moment,
LU-decomp + back-substitution, dispatched to CUBLAS / KA-LU.

The linearized partner needs the derivative of that inverse:

```math
\dot{\mathbf{G}} = -\mathbf{G}\,\dot{(\mathbf{R}\mathbf{R})}\,\mathbf{G}.
```

It **reuses the already-computed** ``\mathbf{G}``. No second LU, no second
back-substitution. The marginal cost of the linearized step is two extra
batched *matmuls* (cheap, well-parallelized, full GPU bandwidth) per core
parameter ``c \in \{\tau, \varpi, \mathbf{Z}\}``, not another inversion.

Compare:

| Approach | Cost ratio (lin+fwd : fwd) |
|---|---|
| **vSmartMOM analytic linearization (this codebase)** | **< 2×** |
| ForwardDiff `Dual{T,V,N}` through the kernels | `(1 + N) ×` (every operation, including the inversion, picks up `N` partial slabs) |
| Finite differences | `(1 + N_state) ×` (one full forward per parameter) |

For a typical OCO-2 retrieval with ``N_\mathrm{state} \approx 30``, this is
the difference between a Jacobian costing ~2 forward runs and one costing
~30. Inside a Levenberg-Marquardt loop that hits the forward model ~10
times to converge, the savings are dispositive — 20 forward equivalents vs
300+.

That property holds at the *kernel* level (one doubling iteration, one
adding step) and aggregates through the entire RT solve. It is preserved
on GPU because both ``\mathbf{G}`` and the matmuls dispatch through the
same `batched_mul` / `batch_inv!` interface — see [Concepts/07](07_architecture.md).

::: info Status — refined AD→analytic boundary in progress
The current production codebase runs the production-fast linearized path
described above; a **cleaner AD-upstream → hand-coded-analytic-downstream
boundary** is in active development. The goal is to make the seam between
the two zones (the `CoreScatteringOpticalPropertiesLin` handoff struct,
δ-M chain-rule expansion, microphysical parameter forwarding) more
idiomatic — fewer manual chain-rule expansions in upstream code,
cleaner extension points for adding new state-vector parameters. None of
this changes the *production performance* claim above (the matrix-inversion
reuse is the kernel-level guarantee and stays put); it cleans up the
ergonomics of the AD boundary. See `docs/dev_notes/` for the work-in-progress
plans.
:::

## The chain rule on adding-doubling

The full derivation is Sanghavi 2014 App. C. The structure (skipping arithmetic):

| Eq. | What it says | Source file |
|---|---|---|
| (C.5)–(C.7) | Differentiation rules for matrix products and inverses | (foundation; used everywhere below) |
| (C.8)–(C.10) | Elemental derivatives ``\dot{\mathbf{T}}_\delta``, ``\dot{\mathbf{R}}_\delta``, ``\dot{\mathbf{J}}_\delta`` w.r.t. the three core layer variables ``(\tau, \varpi_0, \mathbf{Z})`` | `elemental_lin.jl` |
| (C.11)–(C.16) | Doubling/adding derivatives — same shape as the forward Eqs (23)–(28), tangent-linear | `doubling_lin.jl`, `interaction_lin.jl` |
| (C.17)–(C.20) | D-matrix symmetry on derivatives — halves the linearized doubling cost | `doubling_lin.jl` |
| (C.21) | Final assembled derivative form — what `lin_added_layer_all_params!` produces | `lin_added_layer_all_params.jl` |
| (C.22)–(C.24) | Layer-averaged ``\bar{\tau}``, ``\bar{\varpi}_0``, ``\bar{\mathbf{Z}}`` definitions | `compEffectiveLayerProperties.jl` (forward); `compEffectiveLayerProperties_lin.jl` (lin) |
| (C.25)–(C.26) | Chain rule from the elemental SS variables to the microphysical parameters ``(n_r, n_i, r_m, \sigma)`` | `lin_added_layer_all_params.jl` + `Scattering/types_lin.jl` |
| (C.27)–(C.31) | δ derivatives (elemental thickness from `N_doubl`) | `compEffectiveLayerProperties_lin.jl` |
| (C.32)–(C.39) | ``\bar{\varpi}_0`` and ``\bar{\mathbf{Z}}`` derivatives (post-truncation) | same |
| (C.40) | ``\dot{\mathbf{Z}}_m`` from generalized spherical harmonics | `compute_Z_matrices.jl` (linearized variant) |
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
    aerosol_params::Int     # = 7   (τ_ref, n_r, n_i, r_m, σ_g, p₀, σ_p)
    n_aerosols::Int
    n_gases::Int
    n_surface::Int
end

n_total(layout)             # total number of Jacobian columns
aerosol_range(layout, iaer) # column indices for aerosol iaer (length 7)
gas_range(layout)           # column indices for gas VMRs
surface_range(layout)       # column indices for surface params
```

Always use these accessors instead of hand-writing arithmetic like
`7*NAer + NGas + 1`. The seven aerosol parameters per mode are:

| # | Parameter | Meaning |
|---|---|---|
| 1 | `τ_ref` | aerosol optical depth at reference wavelength |
| 2 | `n_r`   | real refractive index |
| 3 | `n_i`   | imaginary refractive index |
| 4 | `r_m`   | median radius of size distribution |
| 5 | `σ_g`   | geometric standard deviation of size distribution |
| 6 | `p₀`    | central pressure of vertical Gaussian |
| 7 | `σ_p`   | width of vertical Gaussian |

A run with two aerosol modes, four gases, and one surface parameter has a
Jacobian with `2·7 + 4 + 1 = 19` columns, and `aerosol_range(layout, 2)`
returns `8:14`.

## What goes through ForwardDiff vs analytic

Strategy by parameter type (as currently implemented and recommended for
new parameters):

| Parameter | Recommended path | Reason |
|---|---|---|
| Lambertian albedo | analytic | trivial: ``\partial r/\partial \rho = 1/\pi`` |
| RPV / Ross-Li / Cox-Munk BRDF | ForwardDiff | low-dimensional, simple surface code |
| `τ_ref` (aerosol OD) | analytic | trivial: ``\partial \tau/\partial \tau_\mathrm{ref} = \tau/\tau_\mathrm{ref}`` |
| `p₀, σ_p` (aerosol height) | analytic | already in `atmo_prof_lin.jl` |
| `n_r, n_i` (refractive index) | analytic Mie | Mie series is AD-hostile (recurrences) |
| `r_m, σ_g` (size distribution) | analytic Mie | same |
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
| Chain-rule expansion to all parameters | [`src/CoreRT/CoreKernel/lin_added_layer_all_params.jl:1–100`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/lin_added_layer_all_params.jl#L1-L100) |
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
