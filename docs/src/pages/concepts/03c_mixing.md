# 3c · Mixing, Scattering vs Total τ, and δ-M Truncation

> **For:** anyone who needs to understand how the per-layer scatterers + absorption combine, why the doubling step works on `τ_scat` and not `τ_λ`, and how the forward-peak is removed before the RT solve.
>
> **Prev:** [3b · Mie & Rayleigh](03b_scattering.md) · **Next:** [4 · The MOM Solver](04_mom_solver.md)

This page closes the layer-optics build. It explains how the ingredients
mix, why the elemental layer is sized by *scattering* optical depth (a
non-obvious design choice that pays off in [Concepts/04](04_mom_solver.md)
and [Concepts/07](07_architecture.md)), and how the forward-peak of the
aerosol phase function is removed before the RT solver sees it.

## Two operators on `CoreScatteringOpticalProperties`

Two methods of `CoreScatteringOpticalProperties` carry the algebraic
semantics that show up in `compEffectiveLayerProperties.jl`:

```julia
# src/CoreRT/types.jl:1063–1093 — mixing scatterers in a layer
function Base.:+(x::CoreScatteringOpticalProperties, y::CoreScatteringOpticalProperties)
    τ  = x.τ .+ y.τ
    wx = x.τ .* x.ϖ                 # scattering-weight from x
    wy = y.τ .* y.ϖ                 # scattering-weight from y
    w  = wx .+ wy
    ϖ  = w ./ τ
    Z⁺⁺ = (wx .* xZ⁺⁺ .+ wy .* yZ⁺⁺) ./ w
    Z⁻⁺ = (wx .* xZ⁻⁺ .+ wy .* yZ⁻⁺) ./ w
    CoreScatteringOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺)
end
```

Mixing is **τϖ-weighted**. The phase matrix from each scattering species is
weighted by its contribution to the *scattered* intensity, not its total
optical depth. That is the only way to conserve total scattered intensity
per stream pair when species are combined. A pure-absorption species
(``\varpi = 0``) carries no weight in the mix — its presence shows up only
through the ``\tau`` increment, not through ``\mathbf{Z}``.

```julia
# src/CoreRT/types.jl:1096+ — vertical concatenation along the spectral axis
function Base.:*(x::CoreScatteringOpticalProperties, y::CoreScatteringOpticalProperties)
    arr_type = array_type(architecture(x.τ))
    x = expandOpticalProperties(x, arr_type)
    y = expandOpticalProperties(y, arr_type)
    CoreScatteringOpticalProperties(
        [x.τ; y.τ], [x.ϖ; y.ϖ],
        cat(x.Z⁺⁺, y.Z⁺⁺, dims=3),
        cat(x.Z⁻⁺, y.Z⁻⁺, dims=3))
end
```

`*` concatenates along the spectral axis — used when several bands share the
same atmospheric column and want to be processed in one call.

## Scattering vs total optical depth — the load-bearing trick

> **The elemental layer is sized by the *scattering* optical depth, not the
> total optical depth. Absorption is layered in afterwards. This is what lets
> a layer with ``\tau_\mathrm{abs} \approx 50`` (deep gas absorption) coexist
> with ``\tau_\mathrm{scat} \approx 0.05`` (thin aerosol) in the same RT
> call without blowing up the doubling count.**

Walk through `construct_atm_layer` in
[`src/CoreRT/tools/atmo_prof.jl:320–380`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/tools/atmo_prof.jl#L320-L380):

1. Initialize cumulative ``\tau``, ``\varpi`` for *scattering only*, plus an
   accumulator ``A = \tau \cdot \varpi`` for the τϖ-weighted ``\mathbf{Z}``.
2. Add Rayleigh scattering: ``\tau_\mathrm{rayl}``, ``\varpi_\mathrm{Cab}``,
   contribution to ``A \cdot \mathbf{Z}`` from the Rayleigh greek source.
3. Add each aerosol with δ-M rescaling already applied (via `createAero`):
   ``\tau_\mathrm{aer}``, ``\tilde{\varpi}``, contribution to ``A \cdot \mathbf{Z}`` from
   that aerosol's truncated Greek matrix.
4. Normalize ``\mathbf{Z}`` by ``A``: ``\mathbf{Z} \leftarrow A \cdot \mathbf{Z} / A``.
5. **Add gas absorption** to get the total optical depth and rescale the
   single-scattering albedo:

```math
\tau_\lambda = \tau_\mathrm{abs} + \tau_\mathrm{scat},
\qquad
\varpi_\lambda = \frac{\tau_\mathrm{scat}\,\varpi}{\tau_\lambda}.
```

The variable named `τ` *inside* `CoreScatteringOpticalProperties` after step 4
is the *scattering* optical depth ``\tau_\mathrm{scat}``. The variable named
`τ_λ` (with the spectral subscript) is the *total*
``\tau_\lambda = \tau_\mathrm{abs} + \tau_\mathrm{scat}``.

The doubling count ``N_\mathrm{doubl}`` is chosen by `get_dtau_ndoubl`
([`src/CoreRT/CoreKernel/rt_kernel.jl:245–253`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/rt_kernel.jl#L245-L253)) from
``\tau_\mathrm{scat} \cdot \varpi`` — **not** from ``\tau_\lambda``:

```julia
function get_dtau_ndoubl(computed_layer_properties, quad_points)
    (; qp_μ) = quad_points
    (; τ, ϖ) = computed_layer_properties      # τ here is τ_scat
    dτ_max  = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])
    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    dτ = τ ./ 2^ndoubl                        # elemental SCATTERING thickness
    return dτ, ndoubl
end
```

This is **SF2023-II Eqs (8)–(9)**: the *constant-`N_doubl` trick*. Two
consequences:

1. A layer with ``\tau_\mathrm{abs} \gg 1`` and small ``\tau_\mathrm{scat}``
   doesn't force the elemental layer to be impossibly thin — ``N_\mathrm{doubl}``
   stays small.
2. Across the spectral grid, ``\tau_\mathrm{abs}`` can swing over orders of
   magnitude (think O₂ A-band line wings) while ``\tau_\mathrm{scat}`` stays
   essentially constant. So ``N_\mathrm{doubl}`` is *the same* at every
   wavelength in the band — and that's exactly what makes the elemental and
   doubling kernels run as one batched matmul over the spectral axis (see
   [Concepts/07](07_architecture.md)).

::: tip Worked example — O₂ A-band layer
A representative O₂ A-band atmospheric layer at the line center:

| Quantity | Value |
|---|---|
| ``\tau_\mathrm{abs}`` (gas absorption) | ≈ 50 |
| ``\tau_\mathrm{scat}`` (Rayleigh + aerosol) | ≈ 0.05 |
| ``\tau_\lambda = \tau_\mathrm{abs} + \tau_\mathrm{scat}`` | ≈ 50.05 |
| ``\varpi`` (scattering-only) | ≈ 0.95 |
| ``\varpi_\lambda = \tau_\mathrm{scat}\varpi / \tau_\lambda`` | ≈ 9.5 × 10⁻⁴ |
| ``N_\mathrm{doubl}`` (sized by ``\tau_\mathrm{scat}\varpi`` ≈ 0.05) | ~6 |

The elemental kernel in [Concepts/04](04_mom_solver.md#elemental-layer)
sees ``\delta\tau_\mathrm{scat} = 0.05 / 2^6 \approx 7.8\times10^{-4}`` for
the scattering-source term and ``\delta\tau_\lambda = 50.05 / 2^6 \approx 0.78``
for the per-wavelength transmission factors. The two enter different
expressions inside the kernel — see Concepts/04.
:::

## δ-M truncation — removing the forward peak

The aerosol phase function has a strong forward peak that requires very high
Fourier order ``L`` to represent faithfully. δ-M truncation (Wiscombe 1977 for
scalar; Sanghavi 2014 App. A and Sanghavi & Stephens 2015 for vector) replaces
the forward peak with a Dirac delta that gets absorbed into the unscattered
direct beam, leaving a smoother truncated phase matrix that the RT solver
can resolve with a moderate number of streams.

The truncation factor (SS2015 Eq. 26):

```math
f_\mathrm{tr} = \frac{\beta_{L_\mathrm{tr}}}{2L_\mathrm{tr}+1}.
```

The optical-depth and SSA rescaling (SS2015 Eq. 8):

```math
\tau^{\!*} = \tau\,(1 - f_\mathrm{tr}\,\varpi_0),
\qquad
\varpi_0^{\!*} = \frac{\varpi_0\,(1-f_\mathrm{tr})}{1 - f_\mathrm{tr}\,\varpi_0}.
```

The modified Greek coefficients (SS2015 Eqs. 27a–f), with ``L_\mathrm{tr}`` the
truncation order:

```math
\beta_l^{\!*} = \frac{\beta_l - \frac{2l+1}{2L_\mathrm{tr}+1}\beta_{L_\mathrm{tr}}}{1-f_\mathrm{tr}}
```

```math
\delta_l^{\!*} = \frac{\delta_l - \frac{2l+1}{2L_\mathrm{tr}+1}\beta_{L_\mathrm{tr}}}{1-f_\mathrm{tr}},\qquad
\gamma_l^{\!*} = \frac{\gamma_l}{1-f_\mathrm{tr}},\qquad
\epsilon_l^{\!*} = \frac{\epsilon_l}{1-f_\mathrm{tr}}.
```

In the code, the rescaling lives in `createAero` (`compEffectiveLayerProperties.jl:67–72`)
for ``\tau`` and ``\varpi``, and in
[`src/CoreRT/LayerOpticalProperties/delta_m_truncation.jl:44–48`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/LayerOpticalProperties/delta_m_truncation.jl#L44-L48) for the
Greek coefficients:

```julia
function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺)
    (; fᵗ, ω̃) = aerosol_optics
    τ_mod = (1 - fᵗ * ω̃) * τAer
    ϖ_mod = (1 - fᵗ) * ω̃ / (1 - fᵗ * ω̃)
    CoreScatteringOpticalProperties(τ_mod, ϖ_mod, AerZ⁺⁺, AerZ⁻⁺)
end
```

### When δ-M alone is not enough: δBGE-fit

δ-M removes the forward peak but introduces error at view angles near the
*exact backscatter* direction (``\theta_\mathrm{view} \approx -60°`` in the
principal plane for typical solar geometry). For hyperspectral retrievals
that fit O₂ A-band line shapes — XCO₂ from OCO-2/3, CH₄ from GOSAT — that
error is large enough to bias retrievals.

SS2015 §3 introduces a vector adaptation of Hu et al. (2000) δ-fit, called
**δBGE-fit** (β-γ-ε fit), which replaces the single δ-m truncation factor
with an SVD-based least-squares fit on the diagonal phase-matrix elements
``b_1(\mu)``, ``b_2(\mu)`` (Eqs. 38a–d). The result preserves the diagonal
elements through the truncation and practically eliminates the error in
``Q``.

::: tip Design choice — δBGE-fit for hyperspectral retrievals
SS2015 §4.1 demonstrates that δ-m truncation distorts O₂ A-band line shapes
near the exact backscatter direction by enough to bias XCO₂ retrievals. The
distortion in ``I`` is small but the distortion in ``Q`` is severe.
δBGE-fit eliminates the ``Q`` error while matching δ-m on ``I``. **Use
δBGE-fit for OCO-2/3 / GOSAT.**
:::

## Code anchors

| Concept | Source |
|---|---|
| Layer-optics builder | [`src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:11–65`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl#L11-L65) |
| Aerosol δ-M wrapper (`createAero`) | `compEffectiveLayerProperties.jl:67–72` |
| δ-M Greek-coefficient rescaling | [`src/CoreRT/LayerOpticalProperties/delta_m_truncation.jl:44–48`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/LayerOpticalProperties/delta_m_truncation.jl#L44-L48) |
| Layer assembly + absorption combination | `src/CoreRT/tools/atmo_prof.jl::construct_atm_layer:320–380` |
| Scattering-only `N_doubl` sizing | `src/CoreRT/CoreKernel/rt_kernel.jl::get_dtau_ndoubl:245–253` |
| Mix scatterers (`+`) | [`src/CoreRT/types.jl:1063–1093`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/types.jl#L1063-L1093) |
| Stack layers along spectral axis (`*`) | `src/CoreRT/types.jl:1096+` |

## Hands-on tutorials

Runnable examples with Plotly figures:

- [CoreRT walkthrough (full layer-optics build)](../tutorials/Tutorial_CoreRT.md)

## References

- **Sanghavi & Stephens (2015)**, *Adaptation of the delta-m and δ-fit truncation methods to vector RT*, JQSRT **159**:53–68, [doi:10.1016/j.jqsrt.2015.03.007](https://doi.org/10.1016/j.jqsrt.2015.03.007). **Definitive vector truncation reference.** Eqs. (8), (26), (27a–f) for δ-m; Eqs. (38a–d) for δBGE-fit.
- **Sanghavi & Frankenberg (2023)**, *Raman scattering Part II*, JQSRT **311**:108791, [doi:10.1016/j.jqsrt.2023.108791](https://doi.org/10.1016/j.jqsrt.2023.108791). Eqs (8)–(9) — the constant-`N_doubl` trick.
- Sanghavi et al. (2014), JQSRT **133**:412–433, [doi:10.1016/j.jqsrt.2013.09.004](https://doi.org/10.1016/j.jqsrt.2013.09.004), App. A. (Vector δ-m derivation.)
- Wiscombe (1977), *The delta-m method: rapid yet accurate radiative flux calculations…*, J. Atmos. Sci. **34**:1408. (Original scalar δ-m.)
- Hu et al. (2000), *δ-fit*, JQSRT **65**:681. (Original δ-fit; SS2015 vectorizes.)
- Crib sheet: `docs/dev_notes/theory_references.md` §D, §F.
