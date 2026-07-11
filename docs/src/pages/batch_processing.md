# Fast Re-runs & Batch Processing

**For:** anyone producing *many* radiative-transfer results that differ in only
a few inputs — surface sweeps, aerosol/trace-gas ensembles, geometry grids,
ML training sets. If you find yourself calling `model_from_parameters` +
`rt_run` inside a loop, read this page first: almost always, most of that work
is recomputing things that did not change.

**Next:** [Quick Start](quickstart.md), [Configure a Scene](IO/Overview.md),
[Surfaces (Concepts)](concepts/05_surfaces.md).

## Where the time goes

A cold `model_from_parameters` + `rt_run` pays for four distinct things:

1. **Line-by-line absorption** — HITRAN parsing + cross-section computation,
   per band × species. Often the single largest setup cost.
2. **Mie scattering** — per aerosol, per band; expensive for wide size
   distributions.
3. **The RT solve** — the Fourier loop over `m`, each moment accumulating
   `Nz` layers of elemental/doubling/interaction over the whole spectral
   batch. This is the per-run cost that dominates once the model is built.
4. **The surface phase** — one surface layer + one final interaction per
   Fourier moment. Cheap relative to (3): roughly `1/Nz` of the solve.

Every tool on this page exists to skip one or more of these when your
scenario axis doesn't touch it.

## Decision table

| What changes between runs | What to call | What is reused (skipped) |
|---|---|---|
| Surface BRDF / albedo only | [`rt_run_atmosphere`](@ref) once, then [`rt_run_surface`](@ref) per surface | (1) + (2) + (3) — only the surface phase re-runs |
| Lambertian albedo, many values or a retrieval | [`lambertian_closure`](@ref) on the cache | (1)–(4): each albedo is a broadcast, exact |
| Aerosol amount / vertical placement | [`update_aerosol_loading!`](@ref) + `rt_run` | (1) + (2) — no Mie recompute |
| Aerosol size distribution / refractive index | [`update_aerosol_microphysics!`](@ref) + `rt_run` | (1) — Mie re-runs for that aerosol only |
| Trace-gas VMRs, T / p / humidity profiles | [`update_model!`](@ref) + `rt_run` | HITRAN parsing + Mie — cross-sections re-evaluate from cached line models |
| SZA or viewing geometry | [`remake_geometry`](@ref) + `rt_run` | (1) + (2) — only quadrature/geometry rebuild |
| Many view angles in one scene | put them all in `vza`/`vaz` (one solve covers them), or [`rt_run_streams`](@ref vSmartMOM.CoreRT.rt_run_streams) for arbitrary post-hoc angles | the whole solve |
| SZA grid × several surfaces | [`ScenarioSweep`](@ref) + [`run_sweep`](@ref) | composes the two rows above |
| Spectral bands, polarization, quadrature type, *number* of species or aerosols | full rebuild (`read_parameters` → `model_from_parameters`) | nothing — these change array shapes everywhere |

## Many surfaces over one atmosphere

The atmosphere accumulation is surface-independent, so build it once and
replay only the surface phase per BRDF — bit-exact against a monolithic
`rt_run`:

```julia
using vSmartMOM

params = read_parameters("config/lambertian_land.yaml")
model  = model_from_parameters(params)

surfaces = [LambertianSurfaceScalar(0.05),
            LambertianSurfaceScalar(0.30),
            rpvSurfaceScalar(0.02, 0.3, 0.7, 0.1)]

# One atmosphere pass; size the cache for the widest BRDF you'll replay:
cache = rt_run_atmosphere(model; target_brdfs = surfaces)

for surf in surfaces
    R_SFI, T_SFI, = rt_run_surface(cache, surf)
    # ... store / process ...
end

# or equivalently, in one call:
results = rt_run_multi_surface(model, surfaces)
```

Each replay costs ~`1/Nz` of a full run. Scope: elastic (`noRS`), single
band, non-canopy surfaces (guarded with `ArgumentError`). A cache built from
a Lambertian-only model holds only the `m = 0` moment — pass non-Lambertian
BRDFs via `target_brdfs` so the cache is sized for them.

## Lambertian albedo sweeps and retrievals

For Lambertian albedos, skip even the replay: the albedo dependence is an
exact scalar closure read off the cache once —

```julia
c = lambertian_closure(cache)

R_03   = c(0.3)                    # O(1) per albedo — a broadcast, no RT
R_spec = c(albedo_vector)          # per-spectral-point albedos
J      = albedo_jacobian(c, 0.3)   # exact ∂R/∂a (no finite differences)
a_hat  = invert_albedo(c, R_meas)  # closed-form albedo retrieval per spectral point
```

This is exact (matches the replay to machine precision), so it is the right
tool for albedo retrievals and dense albedo grids. Solar-beam-only scenes
(a `SurfaceSIF` or per-source emission term makes the closure invalid; the
constructor throws in that case).

## Changing only aerosols

Build a [`BatchContext`](@ref) once; per scenario, update only what changed:

```julia
ctx = BatchContext(params)

# Concentration / vertical placement — cheap, NO Mie recompute:
for τ in (0.05, 0.1, 0.2, 0.5)
    update_aerosol_loading!(ctx, 1; τ_ref = τ)      # aerosol #1
    R, T = rt_run(ctx.model)
end

# Size distribution / refractive index — re-runs Mie for that aerosol only:
new_aer = Aerosol(LogNormal(log(0.2), log(1.8)), 1.4, 1e-6)   # dist, nᵣ, nᵢ
update_aerosol_microphysics!(ctx, 1, new_aer)
R, T = rt_run(ctx.model)
```

(`update_aerosol_loading!` additionally takes `profile_dist` for vertical
placement; `update_aerosol_microphysics!` accepts `τ_ref` alongside the new
`Aerosol` — see the docstrings.)

## Changing trace gases or thermodynamic profiles

`update_model!` re-evaluates absorption cross-sections from *cached* HITRAN
line models — no re-parsing, no Mie:

```julia
ctx = BatchContext(params)

for scene in scenes
    update_model!(ctx;
        T      = scene.T,        # length-Nz temperature [K]
        p_half = scene.p_half,   # length-(Nz+1) half-level pressures [hPa]
        q      = scene.q,        # specific humidity [kg/kg]
        vmr    = scene.vmr)      # Dict, same keys as the original config
    R, T = rt_run(ctx.model)
end
```

Partial updates compose: pass only the keywords that changed; the rest fall
back to the previous scene state. Adding a new *species* (or changing the
band set) needs a fresh `BatchContext`.

## Changing geometry

`remake_geometry` rebuilds only `ObsGeometry` + `QuadPoints` and shares the
optics by reference — bit-exact against a full `model_from_parameters`
rebuild, with zero HITRAN/Mie work:

```julia
for sza in 0:10:70
    m = remake_geometry(model, params; sza = sza)
    R, T = rt_run(m)
end
```

Prefer putting all view angles into one model (`vza`/`vaz` vectors — one
solve covers every angle) and looping only over SZA. For a crossed
SZA × view × surface grid, use the sweep driver, which also amortizes each
SZA's atmosphere across all surfaces:

```julia
sweep  = ScenarioSweep(; sza = [20.0, 40.0, 60.0],
                         vza = [0.0, 30.0], vaz = [0.0, 90.0],
                         view_mode = :cross,
                         brdfs = surfaces)
result = run_sweep(model, params, sweep)   # result.R :: (N_sza, N_view, N_brdf, pol_n, nSpec)
```

## Putting it together: a scenario production loop

The axes compose, but **order matters**: an `AtmosphereRTCache` (and any
`LambertianClosure` built from it) snapshots the optics *at build time* —
every `update_model!` / `update_aerosol_*!` call invalidates it.

The rule: **update the scene → build the cache → sweep surfaces.** Never
reuse a cache across scene updates.

```julia
ctx = BatchContext(params)

for scene in scenes                                   # slowest axis: atmosphere
    update_model!(ctx; T = scene.T, q = scene.q, vmr = scene.vmr)
    update_aerosol_loading!(ctx, 1; τ_ref = scene.aod)

    for sza in scene.szas                             # middle axis: geometry
        m     = remake_geometry(ctx.model, params; sza = sza)
        cache = rt_run_atmosphere(m; target_brdfs = surfaces)   # AFTER all updates

        for surf in surfaces                          # fastest axis: surfaces
            R, = rt_run_surface(cache, surf)
        end
        albedo_grid = lambertian_closure(cache).(0.0:0.05:0.9)  # free extras
    end
end
```

Practical notes:

- **Thread safety:** one `BatchContext` per thread — `update_model!` and
  `rt_run(ctx.model)` mutate shared arrays. Same for a cache: don't share an
  `AtmosphereRTCache` across threads that replay concurrently into it.
- **BLAS threads:** for large spectral batches, cap BLAS
  (`radiative_transfer.numerics.blas_threads: 4` in the config) — oversubscription
  with Julia threads degrades the batched kernels.
- **Cache memory:** a full cache stores 4 matrices of size `NquadN² × nSpec`
  per Fourier moment. Lambertian-only caches are automatically stored slim
  (full matrices at `m = 0` only, ~`(m_max+1)×` smaller); see
  [`AtmosphereRTCache`](@ref).

## What still forces a full rebuild

Changing any of these re-shapes arrays across the whole pipeline — go back to
`read_parameters` → `model_from_parameters` (and a fresh `BatchContext`):

- spectral bands / grids,
- polarization type, quadrature type, `nstreams` / truncation,
- the *number* of aerosols or absorbing species,
- architecture (CPU ↔ GPU) or float type.
