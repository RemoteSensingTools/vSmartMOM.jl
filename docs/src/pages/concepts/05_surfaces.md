# 5 · Surfaces — the bottom-most AddedLayer

> **For:** users configuring scenes with land, ocean, or canopy surfaces; method developers adding a new BRDF.
>
> **Prev:** [4 · The MOM Solver](04_mom_solver.md) · **Next:** [6 · Linearization](06_linearization.md)

In MOM language, the surface is just another layer. Specifically, the
*bottom-most* `AddedLayer` — built once for the lower boundary of the
atmospheric column, then composed with the accumulated atmosphere using the
same `interaction!` machinery that stacks atmospheric layers in
[Concepts/04](04_mom_solver.md).

## Surfaces in MOM language

Each BRDF type is a Julia struct subtyping `AbstractSurfaceType`, plus a
method of `create_surface_layer!` that fills the `r⁻⁺`, `r⁺⁻`, `t⁺⁺`, `t⁻⁻`,
`j₀⁺`, `j₀⁻` matrices of the bottom `AddedLayer` from the BRDF's parameters.
After the atmosphere has been composed TOA → BOA, the loop calls
`interaction!` one more time with the surface layer underneath the
accumulated atmosphere (as a `ScatteringInterface_11` when any atmospheric
layer scatters, `_01` otherwise; see [Concepts/04](04_mom_solver.md)). Then `interaction_hdrf!` ([`src/CoreRT/CoreKernel/interaction_hdrf.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/interaction_hdrf.jl))
runs to extract hemispherical-directional reflectance products for `m=0`.

For opaque lower-boundary surfaces, `t⁻⁻` is zero: an upward field entering
from below the surface would be a sub-surface source, which is not part of the
atmospheric RT problem. The current surface-as-layer implementation keeps
`t⁺⁺` as an identity pass-through so the downwelling field at the surface
remains available for diagnostics such as `T_SFI`, HDRF, and BHR.

```
   CompositeLayer (after final iz=Nz)        Surface BRDF
        │                                          │
        │                                  create_surface_layer!
        │                                          │
        │                                  AddedLayer (bottom)
        │                                          │
        └──────────────► interaction!  ◄───────────┘
                              │
                              ▼
                       interaction_hdrf!
                       HDRF + BHR products
                              │
                              ▼
                       postprocessing_vza!
                              │
                              ▼
                  R, T, hdr, bhr  per (VZA, Stokes, λ)
```

## Surface gallery

The implemented BRDFs and where to find them:

| BRDF | Parameters | Use case | File |
|---|---|---|---|
| `LambertianSurfaceScalar` | albedo (single) | featureless reference, calibration scenes | `Surfaces/lambertian_surface.jl` |
| `LambertianSurfaceSpectrum` | albedo vector | spectrally-varying surface | same |
| `rpvSurfaceScalar` | ρ₀, ρ_c, k, Θ | semi-arid land with hot-spot | `Surfaces/rpv_surface.jl` |
| `RossLiSurfaceScalar` | f_iso, f_vol, f_geo, `zenith_cap_deg` (80° by default; `nothing` = MODIS/VLIDORT kernels) | MODIS-style kernel-driven BRDF | `Surfaces/rossli_surface.jl` |
| `CoxMunkSurface` | wind speed, n_water | sun glint over ocean (polarized) | `Surfaces/coxmunk_surface.jl` |
| `CanopySurface` | LAI, leaf R/T, soil, n_layers | vegetation canopies | `Surfaces/canopy_surface.jl` |

Linearized variants (`*_lin.jl` files) exist for Lambertian and Cox-Munk;
RPV and Ross-Li are forward-only, and canopy linearization is a work in
progress.

## Polarized surfaces (Cox–Munk)

For ocean glint, vSmartMOM uses the Cox & Munk (1954) wave-slope
distribution for wind-roughened sea surfaces, plus full Fresnel reflection
to populate the polarized 4×4 Mueller matrix on each facet. The polarized
`r⁻⁺` blocks are not zero in `Stokes_IQUV` — they couple ``I, Q, U, V``
through the Fresnel angle dependence.

The slope distribution is an isotropic Gaussian whose *total* mean-square
slope is ``0.003 + 0.00512\,U`` (``U`` the 10-m wind speed in m/s); each of
the two slope components carries half of it (`wind_to_sigma2`). Stokes
vectors are transformed explicitly, by a reflecting change of basis, from
each beam's meridian plane into the facet's plane of reflection and back
(`meridian_facet_stokes` in `coxmunk_surface.jl`), so the Mueller matrix follows the
[Stokes convention](../conventions.md) used everywhere else, and its Fourier
moments are stored like those of the atmospheric phase matrices.

Refinements, both on by default (`include_whitecaps`, `shadowing`):

- **Whitecaps:** wind-speed-dependent Lambertian contribution from breaking
  waves (coverage after Monahan & O'Muircheartaigh 1980; albedo
  `whitecap_albedo = 0.22`, Koepke 1984).
- **Shadowing:** Smith (1967) shadowing function for the incident and
  reflected directions, important at grazing geometries.

**Water refractive index.** With the default `n_water = nothing` the Fresnel
coefficients use the Segelstein (1981) index of liquid water at the
wavelength of every spectral point (`water_refractive_index`; n and k
interpolated linearly in wavelength between table rows). A vector `n_water`
gives one index per spectral point. Because the index varies slowly within a
band, the surface operator is evaluated exactly at a few knots and
interpolated linearly in the index between them (`_n_water_knots`). The
knots are placed so that the Fresnel elements of every point are reproduced
to 1e-7 of its reflectance at 11 angles of incidence, and an index with real
part ≤ 1 is always a knot; for water across 0.25–3.9 µm the operator then
agrees to better than 1e-7 relative (1e-6 in Float32). A single index gives
one operator for the whole band, as before. Whitecaps and shadowing do not
depend on the index.

**Exact glint at TOA.** The specular glint needs far more Fourier moments
than the diffuse field. After the Fourier loop the solver therefore adds, at
each requested view, the exact single-bounce glint minus its Fourier series
(`apply_ss_correction!`), with each spectral point's own index. This
correction is applied to the **TOA radiance only**: upwelling radiances at
interior observer heights carry the Fourier-series glint, in the forward and
the linearized solvers alike, so their radiances and Jacobians stay
consistent.

**Ross–Li zenith cap.** `RossLiSurfaceScalar` caps both kernel zenith angles
at `zenith_cap_deg = 80` by default: the LiSparse kernel turns strongly
negative at grazing quadrature nodes, and beyond the cap the BRDF keeps its
value at the cap. This is a vSmartMOM regularization; `zenith_cap_deg =
nothing` evaluates the MODIS (MCD43) kernels at the true angles, as VLIDORT
does. The choice changes TOA intensities by up to 0.8 % and BOA intensities
by up to 1.1 % in the VLIDORT comparison of `test/vlidort_baseline` (Case D).

## Canopy surfaces

`CanopySurface` is a *pre-solver* — it solves a discrete-ordinate problem
within the vegetation layer (LAI-weighted leaf reflectance/transmittance,
soil at the bottom) before presenting an effective lower-boundary
`AddedLayer` to the atmospheric MOM solver. This keeps the canopy physics
(specular leaf scattering, dark soil background, vertical heterogeneity in
LAI) decoupled from the atmospheric pass while reusing the same
matrix-operator language.

## HDRF and BHR products

Many remote-sensing applications report hemispherical-directional
reflectance factor (HDRF) and bi-hemispherical reflectance (BHR) instead
of (or alongside) per-VZA radiances. `interaction_hdrf!`:

```julia
# src/CoreRT/CoreKernel/interaction_hdrf.jl (interaction_hdrf!)
J₀⁺_diffuse = J₀⁺ .- j₀⁺                       # remove the surface's direct beam
hdr_J₀⁻ = r⁻⁺ ⊠ J₀⁺_diffuse .+ j₀⁻             # HDR (per quadrature node)
if m == 0
    bhr_J₀⁻[i, λ] = sum( hdr_J₀⁻[j, λ] * w[j] * μ[j]  for j in nodes )      # BHR up
    bhr_J₀⁺[i, λ] = sum( J₀⁺_diffuse[j, λ] * w[j] * μ[j]  for j in nodes ) +
                    j₀⁺[iμ₀, λ] * μ₀                                         # BHR down
end
```

After the final surface interaction the composite `J₀⁺` also carries the
surface layer's `j₀⁺` (the attenuated direct beam at the μ₀ rows). Its
reflection is already `j₀⁻`, and its flux enters BHR-down explicitly, so it
is removed first; otherwise it is counted twice whenever μ₀ is a weighted
quadrature node (Gauss with the SZA on a node, every Radau run).

The `hdr` output is the per-VZA upwelling radiance just above the surface,
and `bhr_uw` / `bhr_dw` are the hemispherical upwelling / downwelling fluxes
there (only the `m=0` moment contributes). They are **not** normalized by
the incident solar flux: for a Lambertian surface of albedo `A` under a
vacuum, `hdr = μ₀ F₀ A/π`, `bhr_dw = μ₀ F₀` and `bhr_uw = A μ₀ F₀`. Divide
by `μ₀ F₀/π` (HDRF) or by `bhr_dw` (BHR) to form reflectance factors.

**Canopy surfaces: `bhr_dw` is diffuse only.** A `CanopySurface` folds the
direct beam into its reflected source `j₀⁻` and sets its `j₀⁺ = 0`, so the
direct term `j₀⁺[iμ₀] μ₀` above is zero: for canopies `bhr_dw` is the
diffuse (sky) downwelling flux at the canopy top, without the attenuated
direct beam `μ₀ F₀ e^{−τ/μ₀}`. `hdr`, `bhr_uw` and the TOA radiance include
the reflected direct beam. To form a canopy BHR, add the direct flux to
`bhr_dw` before dividing. With `j₀⁺ = 0` the forward BOA downwelling
radiance likewise has no collimated-beam term at a view on the solar
ordinate, which bare surfaces carry and the linearized solver adds for every
surface.

## Sweeping surfaces over a cached atmosphere

Because the surface enters the solve at exactly one point per Fourier
moment — the final `interaction!` after the layer loop — the expensive
atmosphere accumulation is surface-independent and can be paid once:

```julia
cache = rt_run_atmosphere(model)                      # one m-loop, snapshots per moment
R₁, _ = rt_run_surface(cache, LambertianSurfaceScalar(0.05))
R₂, _ = rt_run_surface(cache, rpvSurfaceScalar(0.02, 0.3, 0.7, 0.1))
# or in one call:
results = rt_run_multi_surface(model, [surf_a, surf_b, surf_c])
```

Each replayed surface costs roughly `1/Nz` of a full `rt_run` and is
**bit-exact** against the monolithic result. Size the cache for the widest
BRDF you intend to replay via `target_brdfs` (a Lambertian-built cache
holds only `m = 0` unless told otherwise). Scope: elastic (`noRS`), single
band, non-canopy surfaces.

For **Lambertian albedos specifically**, skip the replay altogether:
the `m=0` surface reflection is rank-one, so the albedo dependence
collapses to an exact scalar closure (Sherman–Morrison),

```math
R(a) = R_{\mathrm{black}} + w_0\,\frac{2a\,E_{\mathrm{dw}}}{1 - a\,\bar S}\;t_{\mathrm{up}},
```

built once from the cache's `m=0` blocks:

```julia
c = lambertian_closure(cache)
R = c(0.3)                        # O(1) per albedo — no matrix algebra
J = albedo_jacobian(c, 0.3)       # exact ∂R/∂a
a = invert_albedo(c, R_measured)  # closed-form albedo retrieval
```

`S̄` is the atmosphere's spherical albedo viewed from below, `E_dw` the
total (direct + diffuse) downward flux at BOA, and `t_up = T⁻⁻u` the
upward transmission of an isotropic BOA source — the same quantities the
classical two-run "MODTRAN" decomposition fits empirically, here read off
one atmosphere run exactly. Lambertian-only caches are automatically
stored slim (`cache_mode = :auto`), keeping the full matrices at `m = 0`
only.

For batch production (many geometries × surfaces), `ScenarioSweep` +
`run_sweep` cross SZA and BRDF axes, rebuilding only
`ObsGeometry`/`QuadPoints` per SZA (`remake_geometry` — shares optics, so
no HITRAN/Mie recomputation) and amortizing each SZA's atmosphere across
all BRDFs.

## Code anchors

| Concept | Source |
|---|---|
| Surface BRDF kernels | `src/CoreRT/Surfaces/` (10 files) |
| `create_surface_layer!` interface | implemented per BRDF in `Surfaces/*.jl` |
| Surface coupling step | [`src/CoreRT/CoreKernel/interaction_hdrf.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/interaction_hdrf.jl) (`interaction_hdrf!`) |
| Polarized Cox–Munk Mueller blocks | [`src/CoreRT/Surfaces/coxmunk_surface.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/Surfaces/coxmunk_surface.jl) + `fresnel.jl` |
| Canopy pre-solver | [`src/CoreRT/Surfaces/canopy_surface.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/Surfaces/canopy_surface.jl) |
| Surface BRDF map (YAML) | `src/IO/Parameters.jl::BRDF_MAP` |
| Linearized Cox–Munk | [`src/CoreRT/Surfaces/coxmunk_surface_lin.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/Surfaces/coxmunk_surface_lin.jl) |
| Atmosphere/surface split + replay | `src/CoreRT/rt_run_split.jl` |
| Analytic Lambertian albedo closure | `src/CoreRT/tools/lambertian_closure.jl` |
| Scenario sweeps (SZA × BRDF) | `src/CoreRT/tools/scenario_sweep.jl` |

To **add a new surface BRDF**, see the developer guide
[Add a Surface BRDF](../extending/surfaces.md).

## Hands-on tutorials

Runnable examples with Plotly figures:

- [Surfaces & BRDF gallery](../tutorials/Tutorial_Surfaces.md)
- [Canopy surface](../tutorials/Tutorial_Canopy.md)

## References

- Sanghavi et al. (2014), JQSRT **133**:412–433, [doi:10.1016/j.jqsrt.2013.09.004](https://doi.org/10.1016/j.jqsrt.2013.09.004). Eqs. (33)–(37) — surface boundary conditions for MOM.
- Wanner, Li, Strahler (1995), *On the derivation of kernels for kernel-driven models of bidirectional reflectance*, JGR **100**:21077. (Ross-thick / Li-dense.)
- Cox, C. & Munk, W. (1954), *Measurement of the roughness of the sea surface from photographs of the sun's glitter*, J. Opt. Soc. Am. **44**:838.
- Rahman, Pinty, Verstraete (1993), *Coupled surface-atmosphere reflectance (CSAR) model* (RPV).
- Koepke, P. (1984), *Effective reflectance of oceanic whitecaps*, Appl. Opt. **23**:1816.
- Crib sheet: `docs/dev_notes/theory_references.md` §C (S2014 Eqs 33–37).
