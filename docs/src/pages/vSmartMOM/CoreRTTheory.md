# Core RT Theory

This page maps the matrix-operator radiative-transfer method to the code paths
that implement it in vSmartMOM. It is meant as a reader's map: enough theory to
identify the equations, and enough source context to find the solver core.

## Solver Spine

For each spectral band, Fourier moment, and atmospheric layer, the forward
elastic path follows one sequence:

```julia
elemental!(...)     # build a thin added layer: r, t, j
doubling!(...)      # double that layer to the requested optical thickness
interaction!(...)   # add it to the accumulated composite atmosphere
```

This sequence is executed in `src/CoreRT/CoreKernel/rt_kernel.jl`. The top-level
loop over bands, layers, and Fourier moments is in `src/CoreRT/rt_run.jl`.

The code uses lowercase fields for a newly added homogeneous layer and uppercase
fields for the accumulated composite layer:

| Concept | Added layer fields | Composite layer fields |
|---|---|---|
| Reflection | `r⁻⁺`, `r⁺⁻` | `R⁻⁺`, `R⁺⁻` |
| Transmission | `t⁺⁺`, `t⁻⁻` | `T⁺⁺`, `T⁻⁻` |
| Source vector | `j₀⁺`, `j₀⁻` | `J₀⁺`, `J₀⁻` |

## Vector RTE

The method starts from the plane-parallel vector radiative-transfer equation
for a Stokes vector ``\mathbf{L}`` with a direct solar beam:

```math
\mu\frac{d\mathbf{L}}{d\tau}
= -\mathbf{L}
+ (1-\varpi_0)\mathbf{B}(T)
+ \frac{\varpi_0}{4\pi}\mathbf{Z}\mathbf{S}_0 e^{-\tau/\mu_0}
+ \frac{\varpi_0}{4\pi}
\int_0^{2\pi}\int_{-1}^{1}
\mathbf{Z}(\mu,\phi;\mu',\phi')\mathbf{L}(\mu',\phi')
\,d\mu'\,d\phi' .
```

The vector formulation and linearization follow Sanghavi et al. (2014). The
scalar predecessor is Sanghavi et al. (2013).

### Polarization State

vSmartMOM carries Stokes subsets through fixed polarization types: `Stokes_I`,
`Stokes_IQ`, `Stokes_IQU`, and `Stokes_IQUV`. The phase-matrix moments,
quadrature stream blocks, surface matrices, and kernel work arrays are sized
from that type, so scalar and vector runs share the same operator structure
without runtime branching over unused Stokes components.

## Discretization

After Fourier decomposition in azimuth and quadrature discretization in polar
angle, each Fourier moment becomes a finite matrix problem for upward and
downward streams:

```math
\frac{d}{d\tau}
\begin{bmatrix}
\mathbf{l}_m^+ \\
\mathbf{l}_m^-
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{A}^{++}_m & \mathbf{A}^{+-}_m \\
\mathbf{A}^{-+}_m & \mathbf{A}^{--}_m
\end{bmatrix}
\begin{bmatrix}
\mathbf{l}_m^+ \\
\mathbf{l}_m^-
\end{bmatrix}.
```

The quadrature state is `QuadPoints`, built by `rt_set_streams` in
`src/CoreRT/tools/rt_set_streams.jl`. The user-facing selector comes from the
`quadrature_type` field in YAML/TOML and is parsed through `QUAD_MAP` in
`src/IO/Parameters.jl`.

Implemented quadrature choices:

- `GaussQuadHemisphere()`: Gaussian nodes on ``[0,1]`` plus camera and solar
  directions as zero-weight nodes.
- `GaussQuadFullSphere()`: full-sphere Gaussian construction, then uses the
  positive-hemisphere nodes.
- `RadauQuad()`: block-Radau construction that includes the solar direction,
  matching the direct-beam treatment in Sanghavi et al. (2014), Appendix B.

## Elemental Layer

Source file: `src/CoreRT/CoreKernel/elemental.jl`

Main entry point: `elemental!`

Kernels:

- `get_elem_rt!`: reflection/transmission for one thin layer.
- `get_elem_rt_SFI!`: direct-solar source-function integration for the same
  layer.

For a thin elemental layer of optical thickness ``\delta``, the matrix-operator
method starts from first-order reflection and transmission operators
(Sanghavi et al. 2014, Eqs. 19-20):

```math
\mathbf{T}_\delta
= \mathbf{I}
+ \left[-\mathbf{M}^{-1}
+ \frac{\varpi_0}{2}
\mathbf{M}^{-1}\overline{\mathbf{Z}}^{++}\mathbf{C}
\right]\delta ,
```

```math
\mathbf{R}_\delta
= \frac{\varpi_0}{2}
\mathbf{M}^{-1}\overline{\mathbf{Z}}^{-+}\mathbf{C}\delta .
```

In code, `get_elem_rt!` uses finite-thickness single-scattering expressions
rather than only the limiting formula. The important terms are:

```julia
r⁻⁺[i,j,n] = ϖ_λ[n] * Z⁻⁺[i,j,n2] *
             (μ[j] / (μ[i] + μ[j])) * wct[j] *
             -expm1(-dτ_λ[n] * (1 / μ[i] + 1 / μ[j]))

t⁺⁺[i,j,n] = ϖ_λ[n] * Z⁺⁺[i,j,n2] *
             (μ[j] / (μ[i] - μ[j])) * wct[j] *
             expdiff_neg(dτ_λ[n] / μ[i], dτ_λ[n] / μ[j])
```

The `-expm1(-x)` and `expdiff_neg(a, b)` forms are deliberate numerical
stabilizers for optically thin layers, especially in `Float32`.

The elemental optical-thickness bound is chosen before this point by
`doubling_number` in `src/CoreRT/tools/rt_helper_functions.jl`. The result is
`dτ` plus `ndoubl`, so a sufficiently thin layer can be doubled back to the
physical layer thickness.

## Solar SFI

Current code computes solar source vectors inside the elemental kernel. It does
not yet implement the proposed unified offline source-function architecture, and
thermal emission is not currently a parallel offline source in the core.

`get_elem_rt_SFI!` implements the direct-beam single-scattering integrals from
Fell (1997), Eqs. 1.52-1.54. The code first forms phase-matrix times solar
Stokes vector:

```julia
Z⁺⁺_I₀ += Z⁺⁺[i, ii, n2] * F₀[ii-i_start+1, n2]
Z⁻⁺_I₀ += Z⁻⁺[i, ii, n2] * F₀[ii-i_start+1, n2]
```

The upwelling source uses a special limit when the quadrature stream is the
solar stream:

```julia
J₀⁺[i, 1, n] =
    wct02 * ϖ_λ[n] * Z⁺⁺_I₀ *
    (dτ_λ[n] / μ[i]) * exp(-dτ_λ[n] / μ[i])
```

For other streams it uses the stabilized exponential difference:

```julia
J₀⁺[i, 1, n] =
    wct02 * ϖ_λ[n] * Z⁺⁺_I₀ *
    (μ[i_start] / (μ[i] - μ[i_start])) *
    expdiff_neg(dτ_λ[n] / μ[i], dτ_λ[n] / μ[i_start])
```

The downwelling source uses:

```julia
J₀⁻[i, 1, n] =
    wct02 * ϖ_λ[n] * Z⁻⁺_I₀ *
    (μ[i_start] / (μ[i] + μ[i_start])) *
    -expm1(-dτ_λ[n] * (1 / μ[i] + 1 / μ[i_start]))
```

Finally, both source vectors are multiplied by the direct-beam attenuation from
the top of atmosphere to the top of the current layer:

```julia
J₀⁺[i, 1, n] *= exp(-τ_sum[n] / μ[i_start])
J₀⁻[i, 1, n] *= exp(-τ_sum[n] / μ[i_start])
```

## Doubling

Source file: `src/CoreRT/CoreKernel/doubling.jl`

Main entry points: `doubling!`, `doubling_helper!`

Repeated doubling uses the same adding equations with a layer interacting with
an identical copy of itself. The multiple-reflection factor is the matrix
geometric series:

```math
(\mathbf{I} - \mathbf{R}\mathbf{R})^{-1}.
```

In code the expensive batched inverse and its products are grouped in:

```julia
compute_geometric_progression!(temp1, tt⁺⁺_gp_refl,
                               r⁻⁺, t⁺⁺,
                               I_static, temp2,
                               temp1_ptr, temp2_ptr)
```

Then each doubling step updates the source and scattering operators:

```julia
doubling_source_update!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt⁺⁺_gp_refl, expk)
doubling_rt_update!(r⁻⁺, t⁺⁺, tt⁺⁺_gp_refl, expk)
```

After the final doubling, `apply_D_matrix!` and `apply_D_matrix_SFI!` use the
polarization symmetry matrix ``\mathbf{D} = \mathrm{diag}(1,1,-1,-1)`` per
stream to fill the reverse-direction operators. This corresponds to the
D-matrix symmetry relations in Sanghavi et al. (2014), Eqs. 29-32.

## Adding / Interaction

Source file: `src/CoreRT/CoreKernel/interaction.jl`

Main entry point: `interaction!`

The code dispatches on a `ScatteringInterface_*` type chosen by
`get_scattering_interface`:

- `ScatteringInterface_00`: neither accumulated composite nor added layer
  scatters.
- `ScatteringInterface_01`: the added layer scatters and the composite does not.
- `ScatteringInterface_10`: the composite scatters and the added layer does not.
- `ScatteringInterface_11`: both scatter; this is the full matrix-operator
  adding case.

The full case implements Sanghavi et al. (2014), Eqs. 23-28. With the current
code convention, the added layer is the lower layer and the composite layer is
the upper layer. The two geometric factors are:

```julia
temp2 .= I_static .- r⁻⁺ ⊠ R⁺⁻
batch_inv!(temp1, temp2, temp1_ptr, temp2_ptr)
T01_inv = T⁻⁻ ⊠ temp1

temp2 .= I_static .- R⁺⁻ ⊠ r⁻⁺
batch_inv!(temp1, temp2, temp1_ptr, temp2_ptr)
T21_inv = t⁺⁺ ⊠ temp1
```

The composite source and R/T operators are then updated in place:

```julia
J₀⁻ .= J₀⁻ .+ T01_inv ⊠ (r⁻⁺ ⊠ J₀⁺ .+ j₀⁻)
R⁻⁺ .= R⁻⁺ .+ T01_inv ⊠ r⁻⁺ ⊠ T⁺⁺
T⁻⁻ .= T01_inv ⊠ t⁻⁻

J₀⁺ .= j₀⁺ .+ T21_inv ⊠ (J₀⁺ .+ R⁺⁻ ⊠ j₀⁻)
T⁺⁺ .= T21_inv ⊠ T⁺⁺
R⁺⁻ .= r⁺⁻ .+ T21_inv ⊠ R⁺⁻ ⊠ t⁻⁻
```

This is the heart of the adding-doubling solver: elemental layers become full
homogeneous layers by doubling, then `interaction!` stacks those layers into the
inhomogeneous atmosphere.

## Truncation

Truncation reduces the high-order forward peak before the RT solve. vSmartMOM
uses vector delta-m / delta-fit style corrections following Sanghavi and
Stephens (2015), while the Greek-coefficient representation follows Sanghavi
et al. (2014), Appendix A.

Primary source files:

- `src/Scattering/truncate_phase.jl`
- `src/CoreRT/tools/atmo_prof.jl`
- `src/CoreRT/LayerOpticalProperties/`

For a truncation order ``L_\mathrm{tr}``, the scalar truncation factor is
Sanghavi and Stephens (2015), Eq. 26:

```math
f_\mathrm{tr} =
\frac{\beta_{L_\mathrm{tr}}}{2L_\mathrm{tr}+1}.
```

The forward peak is removed from the Greek coefficients and the remaining
coefficients are renormalized. For the coefficient families affected by the
truncation term, Sanghavi and Stephens (2015), Eqs. 27a-27d, give:

```math
\beta_l^* =
\frac{\beta_l - \frac{2l+1}{2L_\mathrm{tr}+1}\beta_{L_\mathrm{tr}}}
{1 - f_\mathrm{tr}},
\qquad
\delta_l^* =
\frac{\delta_l - \frac{2l+1}{2L_\mathrm{tr}+1}\beta_{L_\mathrm{tr}}}
{1 - f_\mathrm{tr}}.
```

```math
\gamma_l^* =
\frac{\gamma_l}{1 - f_\mathrm{tr}},
\qquad
\epsilon_l^* =
\frac{\epsilon_l}{1 - f_\mathrm{tr}}.
```

The layer optical properties are adjusted consistently:

```math
\tau^* = \tau(1 - f_\mathrm{tr}\varpi_0),
\qquad
\varpi_0^* =
\frac{\varpi_0(1 - f_\mathrm{tr})}
{1 - f_\mathrm{tr}\varpi_0}.
```

The resulting truncation factor `fᵗ`, truncated Greek coefficients, optical
depth, and single-scattering albedo enter layer assembly through
`construct_atm_layer`, where the atmospheric optical properties are rescaled
before `rt_kernel!` sees them.

## Linearization

Linearized RT follows Sanghavi et al. (2014), Appendix C. The public workflow is
documented in [Compute Jacobians](../jacobians.md).

The matrix-operator formalism makes the linearization exact at the operator
level: each elemental, doubling, and interaction update has a corresponding
tangent-linear update. This is the source of the "linearization is an exact
method" statement in the original vSmartMOM papers, and it is the reason the
hand-coded path remains useful alongside AD-based optical-property derivatives.

Main files:

- `src/CoreRT/tools/lin_model_from_parameters.jl`: builds the forward model and
  derivative optical-property containers.
- `src/CoreRT/CoreKernel/elemental_lin.jl`: derivatives of elemental `r`, `t`,
  and `j`.
- `src/CoreRT/CoreKernel/doubling_lin.jl`: derivative propagation through
  doubling.
- `src/CoreRT/CoreKernel/interaction_lin.jl`: derivative propagation through
  adding.
- `src/CoreRT/CoreKernel/lin_added_layer_all_params.jl`: assembles derivative
  matrices for all active parameter blocks.
- `src/CoreRT/parameter_layout.jl`: names the final Jacobian parameter ordering.

Use `ParameterLayout` accessors such as `aerosol_range`, `gas_range`, and
`surface_index` instead of manually writing offsets like `7*NAer + NGas + 1`.

## Inelastic Extension

The Raman-aware path extends the elastic solver rather than replacing it. The
optical-property theory follows Sanghavi (2022), and the RT coupling follows
Sanghavi and Frankenberg (2023).

Main files:

- `src/Inelastic/`: molecular constants, Cabannes/Raman phase terms, and Raman
  cross sections.
- `src/CoreRT/CoreKernel/elemental_inelastic.jl` and
  `elemental_inelastic_plus.jl`: inelastic elemental source and coupling terms.
- `src/CoreRT/CoreKernel/doubling_inelastic.jl`: inelastic doubling.
- `src/CoreRT/CoreKernel/interaction_inelastic.jl`: inelastic adding.
- `src/CoreRT/rt_run.jl`: dispatch from `RRS`, `VS_0to1`, `VS_1to0`, and
  `_plus` variants into the matching kernel path.

The current implementation uses a linear-in-inelastic-scattering approximation:
only one inelastic event is included along a photon path. Elastic multiple
scattering still uses the full adding-doubling machinery.

## Surface Interaction and HDRF

Surfaces enter as lower-boundary BRDF operators and are coupled through the same
matrix-operator language used for layer adding. Scalar surface models
(`LambertianSurface*`, RPV, Ross-Li) populate the scalar Stokes block, while
polarized surface models such as Cox-Munk populate Mueller-matrix blocks from
Fresnel reflection over a wind-roughened facet distribution.

The relevant code lives in `src/CoreRT/Surfaces/`,
`src/CoreRT/CoreKernel/interaction_hdrf.jl`, and the postprocessing utilities in
`src/CoreRT/tools/postprocessing_vza*.jl`. Canopy surfaces solve their internal
vegetation/soil stack before presenting an effective lower boundary to the
atmospheric solver. HDRF/BHR postprocessing then converts the accumulated
operator state into view-angle products and hemispheric quantities.

## Other Dispatch Arms

Several specialized paths share the same operator language:

- `rt_kernel_ss.jl`: single-scattering-oriented kernels.
- `rt_kernel_multisensor.jl`: sensor-level output through the atmosphere.
- `interaction_hdrf.jl` and surface kernels: surface coupling and HDRF/BHR
  post-processing.
- `rt_kernel_canopy!`: swaps in canopy directional scattering while keeping the
  elemental/doubling/interaction sequence.

These are extensions around the same core state: `R`, `T`, and source vectors
for each direction.
