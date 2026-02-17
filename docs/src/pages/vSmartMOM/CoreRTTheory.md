# Core RT Theory: Elemental, Doubling, and Adding

This page maps the core radiative transfer operators in `vSmartMOM` to their mathematical sources and implementation files.

## Paper lineage

1. Sanghavi et al. (2013, JQSRT 116): scalar matrix-operator formulation and linearization.
2. Sanghavi et al. (2014, JQSRT 133): vectorized `vSmartMOM` formulation used by this codebase.
3. Sanghavi and Stephens (2015, JQSRT 159): vector truncation (`delta-m`, `delta-fit`, `delta-BGE`).
4. Fell (1997 thesis, German): source-function expressions used in code comments for SFI terms.

## Layer operator form

The core layer system is represented as

```math
\begin{bmatrix}
I_{\Delta}^{+} \\
I_{0}^{-}
\end{bmatrix}
=
\begin{bmatrix}
J_{\Delta}^{+} \\
J_{0}^{-}
\end{bmatrix}
+
\begin{bmatrix}
T_{\Delta 0} & R_{\Delta 0} \\
R_{0 \Delta} & T_{0 \Delta}
\end{bmatrix}
\begin{bmatrix}
I_{0}^{+} \\
I_{\Delta}^{-}
\end{bmatrix},
```

which corresponds to Eq. (16) in Sanghavi et al. (2014), and Eq. (16) in Sanghavi et al. (2013) for the scalar case.

## Symbol map: paper -> code

During adding, the code treats the added layer as `2-1` and the existing composite as `1-0`.

| Paper symbol | Meaning | Code variable |
|---|---|---|
| `R_{21}` | reflection of added layer | `added_layer.r⁻⁺` |
| `R_{12}` | reverse reflection of added layer | `added_layer.r⁺⁻` |
| `T_{21}` | transmission (downward convention in code) | `added_layer.t⁺⁺` |
| `T_{12}` | reverse transmission | `added_layer.t⁻⁻` |
| `J_{21}` | source term (added layer) | `added_layer.j₀⁺` |
| `J_{12}` | source term (added layer, reverse) | `added_layer.j₀⁻` |
| `R_{10}` | reflection of existing composite | `composite_layer.R⁻⁺` |
| `R_{01}` | reverse reflection of existing composite | `composite_layer.R⁺⁻` |
| `T_{10}` | transmission of existing composite | `composite_layer.T⁺⁺` |
| `T_{01}` | reverse transmission of existing composite | `composite_layer.T⁻⁻` |
| `J_{10}` | source term (composite) | `composite_layer.J₀⁺` |
| `J_{01}` | source term (composite, reverse) | `composite_layer.J₀⁻` |

## Elemental layer

Code: `src/CoreRT/CoreKernel/elemental.jl`

- Entry point: `elemental!`
- Kernels: `get_elem_rt!` and `get_elem_rt_SFI!`

The elemental kernel computes finite-thickness single-scattering `r` and `t` terms directly in stream space (the finite-thickness form corresponding to the infinitesimal-layer equations in Sanghavi et al. (2014), Eq. (19)).

The SFI source terms (`j₀⁺`, `j₀⁻`) are implemented in `get_elem_rt_SFI!`, with comments explicitly referencing Fell thesis Eq. 1.52-1.54.

### Fell (1997) mapping for SFI

From the provided Fell section (*1.4.1 Einfachstreuung*), the key single-scattering edge expressions are:

- Eq. (1.52): top-edge singly scattered term `L^{ES}(0;-\vec{s})`
- Eq. (1.53): bottom-edge singly scattered term `L^{ES}(\tau;\vec{s})` for `\mu \neq \mu_0`
- Eq. (1.54): limiting bottom-edge form for `\mu = \mu_0`
- Eq. (1.55)-(1.56): optically thin (`\Delta\tau \ll 1`) linearized forms

Direct implementation mapping in `get_elem_rt_SFI!`:

- `J₀⁻` expression follows Fell Eq. (1.52):
  - code: `J₀⁻[...] = ... (μ0/(μ+μ0)) * (1 - exp(-dτ*(1/μ + 1/μ0)))`
- `J₀⁺` expression for `μ != μ0` follows Fell Eq. (1.53):
  - code: `J₀⁺[...] = ... (μ0/(μ-μ0)) * (exp(-dτ/μ) - exp(-dτ/μ0))`
- `J₀⁺` expression for `μ = μ0` follows Fell Eq. (1.54):
  - code branch when stream index is the solar stream: `... (dτ/μ0) * exp(-dτ/μ0)`

Additional note:

- The code applies `exp(-τ_sum/μ0)` after these terms to account for attenuation by overlying layers before current-layer source injection in the adding workflow.
- The thin-layer formulas (Fell Eq. (1.55)-(1.56)) correspond to first-order small-`dτ` approximations and are conceptually aligned with the optional thin-layer elemental path in `elemental!`.

## Doubling

Code: `src/CoreRT/CoreKernel/doubling.jl`

Entry point: `doubling_helper!`

For homogeneous layers, repeated doubling applies the adding equations recursively, with geometric-series factors of the form

```math
(I - R R)^{-1}.
```

The update structure follows the Eq. (23)-(28) family in Sanghavi et al. (2014) (same structure as the scalar Eq. (22)-(27) family in Sanghavi et al. (2013)).

After doubling, the code applies the polarization symmetry transform with the D-matrix (`D = diag(1,1,-1,-1)` per stream), consistent with Eq. (29)-(32) symmetry relations in Sanghavi et al. (2014):

- `apply_D_matrix!`
- `apply_D_matrix_SFI!`

## Adding / interaction

Code: `src/CoreRT/CoreKernel/interaction.jl`

The full scattering case (`ScatteringInterface_11`) implements the bidirectional adding equations using two inverses:

```math
(I - R_{21}R_{01})^{-1}, \quad (I - R_{01}R_{21})^{-1},
```

and updates `J`, `R`, `T` in both directions. This is the direct matrix-operator adding form from Sanghavi et al. (2014), Eq. (23)-(28).

The other interfaces are algebraic reductions:

- `ScatteringInterface_00`: no scattering in either layer
- `ScatteringInterface_01`: added layer scatters, composite does not
- `ScatteringInterface_10`: composite scatters, added layer does not

## Where workflow is executed

- Per-layer optical properties and interface bookkeeping:
  - `src/CoreRT/tools/atmo_prof.jl`
  - `src/CoreRT/tools/rt_helper_functions.jl` (`doubling_number`, `get_scattering_interface`)
- Core per-layer execution:
  - `src/CoreRT/CoreKernel/rt_kernel.jl`
  - sequence: `elemental! -> doubling! -> interaction!`
- Top-level loop over Fourier moments, layers, and bands:
  - `src/CoreRT/rt_run.jl`

## Truncation coupling

Truncation math is implemented in `src/Scattering/truncate_phase.jl`:

- Eq. (38a)-(38d) from Sanghavi and Stephens (2015) map to updates of `βᵗ, δᵗ, αᵗ, ζᵗ`.
- Eq. (39)-style least-squares systems are used for `γ` and `ϵ` fits.

The resulting truncation factor `fᵗ` enters layer assembly in `construct_atm_layer` (`src/CoreRT/tools/atmo_prof.jl`) via `(1 - fᵗ)` weighting and subsequent optical-property rescaling.
