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
accumulated atmosphere. Then `interaction_hdrf!` ([`src/CoreRT/CoreKernel/interaction_hdrf.jl:1–42`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/interaction_hdrf.jl#L1-L42))
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
| `RossLiSurfaceScalar` | f_iso, f_vol, f_geo | MODIS-style kernel-driven BRDF | `Surfaces/rossli_surface.jl` |
| `CoxMunkSurface` | wind speed, n_water | sun glint over ocean (polarized) | `Surfaces/coxmunk_surface.jl` |
| `CanopySurface` | LAI, leaf R/T, soil, n_layers | vegetation canopies | `Surfaces/canopy_surface.jl` |

Linearized variants (`*_lin.jl` files) exist for Lambertian, RPV, Ross-Li,
and Cox-Munk; canopy linearization is a work in progress.

## Polarized surfaces (Cox–Munk)

For ocean glint, vSmartMOM uses the Cox & Munk (1954) wave-slope
distribution for wind-roughened sea surfaces, plus full Fresnel reflection
to populate the polarized 4×4 Mueller matrix on each facet. The polarized
`r⁻⁺` blocks are not zero in `Stokes_IQUV` — they couple ``I, Q, U, V``
through the Fresnel angle dependence.

Optional refinements:

- **Whitecaps:** wind-speed-dependent Lambertian contribution from breaking
  waves (Koepke 1984).
- **Shadowing:** Mishchenko-Travis shadowing function for grazing geometries.

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
# src/CoreRT/CoreKernel/interaction_hdrf.jl:1–42
hdr_J₀⁻ = r⁻⁺ ⊠ J₀⁺ .+ j₀⁻                     # HDRF (per quadrature node)
if m == 0
    bhr_J₀⁻[i, λ] = sum( hdr_J₀⁻[j, λ] * w[j] * μ[j]  for j in nodes )  # BHR
end
```

HDRF is the per-VZA (or per-quadrature) reflectance into the upper
hemisphere normalized by incident solar flux; BHR is its hemispherical
integral, only well-defined for `m=0`.

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
| Surface coupling step | [`src/CoreRT/CoreKernel/interaction_hdrf.jl:1–42`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/interaction_hdrf.jl#L1-L42) |
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
