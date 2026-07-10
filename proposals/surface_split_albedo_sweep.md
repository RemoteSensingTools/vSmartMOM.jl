# PR Sketch — Surface sweeps without re-running the atmospheric adding/doubling

**Goal:** evaluate many surfaces (albedos, BRDF parameter sets) against one
atmosphere, paying the m-loop layer accumulation (elemental → doubling →
interaction over `Nz` layers) exactly once.

**Status:** proposal. A working prototype of the general mechanism already
exists on the `gchp-io` branch (`rt_run_atmosphere` / `rt_run_surface` /
`rt_run_multi_surface`); this sketch is (1) the plan to extract and harden it
onto `main`, and (2) an analytic Lambertian closure layered on top that makes
continuous albedo sweeps and albedo retrieval O(1) per evaluation.

**Author context:** surveyed `main` (post-RamanLuT, commit `8afe0aec`) and the
`gchp-io` prototype (bundled inside commit `8df0d67a`, 2026-05-15). Related
student-facing material: `proposals/gchp_aerosol_optics/STUDENT_BRIEF.md`
(the ML-training-set use case is the main consumer of these sweeps).

---

## 0. Why this works: the architecture already decouples the surface

In `rt_run`'s Fourier loop (`src/CoreRT/rt_run.jl`), the surface enters the
solve at exactly one point per moment `m`: after the z-loop finishes, a
surface `AddedLayer` is built (`create_surface_layer!`, rt_run.jl:475-494) and
fused into the composite with **one** final `interaction!` (rt_run.jl:502-508).
Everything before that — the expensive part, `Nz` layers × (elemental +
doubling + per-layer interaction), per moment, per spectral point — is
completely surface-independent.

Cost model per additional surface, once the atmosphere is cached:

| phase | monolithic `rt_run` | cached replay |
|---|---|---|
| layer accumulation | `(m_max+1) · Nz` kernel triples | 0 |
| surface build + final interaction | `(m_max+1) · 1` | `(m_max+1) · 1` |
| postprocessing | per m | per m |

For a typical `Nz ≈ 20–70`, each extra surface costs roughly `1/Nz` of a full
run (plus snapshot-restore copies). For Lambertian-only sweeps the analytic
closure in §3 removes even that.

---

## 1. Prior art on `gchp-io` (the starting point)

What exists (all inside commit `8df0d67a` "Add GCHP TOMAS aerosol
diagnostics" — bundled with ~45 unrelated files, so **extract, don't rebase**):

- **`src/CoreRT/rt_run_split.jl`** (297 lines, new):
  - `rt_run_atmosphere(model; target_brdfs, i_band, sources) -> AtmosphereRTCache`
    — runs `rt_run` with `stop_after_atmosphere = true` and an
    `atm_snapshot_callback` that deep-copies the pre-surface composite blocks
    (`R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻`, per-source slots) once per Fourier
    moment, plus the surface scattering-interface trait and `τ_sum` at BOA.
    Cache `m_max` is sized as `max(model m_max, component_m_max over
    target_brdfs)` so a Lambertian model can still serve Cox-Munk replays.
  - `rt_run_surface(cache, brdf) -> (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw)`
    — per m: `copyto!` the snapshot into a fresh composite, run
    `create_surface_layer!` → `inject_surface_SIF!` →
    `surface_source_contribute!` → `interaction!` → `interaction_hdrf!` →
    `postprocessing_vza!`/`_hdrf!`, mirroring rt_run's surface block
    line-for-line; Cox-Munk TMS `apply_ss_correction!` replayed per BRDF.
  - `rt_run_multi_surface(model, brdfs)` = one `rt_run_atmosphere` + a
    `rt_run_surface` per BRDF.
- **`AtmosphereRTCache`** (gchp-io `src/CoreRT/types.jl:466-491`): fully
  concrete type parameters; snapshot arrays stay device-resident (deepcopy of
  a `CuArray` stays on GPU).
- **Hooks in `rt_run`** (+30/−4 lines): three new kwargs
  (`atm_snapshot_callback`, `stop_after_atmosphere`, `m_max_override`), the
  callback firing per m **after** the z-loop and **before**
  `create_surface_layer!`; `stop_after_atmosphere` skips surface +
  postprocessing via `continue`; `m_max_override` may only raise `m_max`.
- **`src/CoreRT/tools/scenario_sweep.jl`** (206 lines):
  `ScenarioSweep`/`SweepResult`/`run_sweep` — SZA × view-pair × BRDF batching;
  rebuilds the model per SZA (`model_from_parameters`), amortizes the
  atmosphere across BRDFs within each SZA via `rt_run_multi_surface`. Consumed
  by `src/IO/Benchmark/netcdf_writer.jl` (ML training-set generation).

Declared scope on the branch: `noRS` only, single band, Canopy rejected with
`ArgumentError`.

**What's not ideal yet (gaps this PR must close):**

1. **Zero tests.** The docstring's acceptance criterion —
   `rt_run_multi_surface(model, [model.surfaces[1]])[1]` bit-exact vs
   `rt_run(model)` — is never checked anywhere.
2. **73 commits of drift** since the merge-base (2026-05-13): fused LU
   interaction (`a2f08e74`, `598c9fef`), vectorized VZA weighting
   (`4c8d1a2e`), aerosol type-stability (`a332d6d8`), the RamanLuT merge. The
   rt_run.jl hook diff is small (+30/−4) so re-application should be clean,
   but `rt_run_surface`'s mirrored surface block must be re-mirrored against
   *today's* rt_run.jl, not May's.
3. **Bundled commit** — needs extraction into a dedicated, reviewable PR.
4. **Memory footprint** is undocumented:
   `bytes ≈ (m_max+1) · sizeof(FT) · nSpec · NquadN · (4·NquadN + 2)`
   (device-resident). E.g. Float64, NquadN=60, nSpec=10⁴, m_max=20:
   ≈ 24 GB — fine for small nSpec / scalar runs, prohibitive for full spectral
   bands on GPU. §4 has the Lambertian slim-cache mitigation.
5. `Vector{Any}` accumulation + post-hoc narrowing in the capture callback —
   works, but tighten while porting.
6. `scenario_sweep.jl` rebuilds the full model per SZA; should compose with
   `BatchContext` (deferred, noted in the branch itself).

---

## 2. PR 1 (the small PR): port the split to `main`, with tests

**Branch:** `feat/surface-split` off `main`. Contents:

1. **rt_run.jl hooks** — re-apply the 3 kwargs + callback site. Insertion
   point today: immediately before the `if brdf isa CanopySurface` surface
   block (rt_run.jl:474), symmetric to the existing post-surface
   `streams_callback` (Phase H, rt_run.jl:553). Guard the Cox-Munk
   TMS block with `!stop_after_atmosphere` as on the branch.
   Cost when unused: one branch per m, no allocations, bit-exact backwards
   compat (same contract as `streams_callback`).
2. **`AtmosphereRTCache`** into `src/CoreRT/types.jl` (or a new
   `types_surface_split.jl` included from CoreRT.jl), and
   **`src/CoreRT/rt_run_split.jl`** re-mirrored against current rt_run.jl.
   Keep the branch's scope guards (noRS-only, no canopy, single band) as
   explicit `ArgumentError`s.
3. **Exports:** `rt_run_atmosphere`, `rt_run_surface`, `rt_run_multi_surface`,
   `AtmosphereRTCache`.
4. **Tests — `test/test_surface_split.jl`** (the piece the prototype never had):
   - *Identity:* `rt_run_multi_surface(model, [model.surfaces[1]])[1]` vs
     `rt_run(model)` for (a) Lambertian Rayleigh-only, (b) Lambertian +
     aerosol, (c) RPV, (d) Cox-Munk (exercises the replayed TMS correction).
     Same-process, same-platform, same code path ⇒ assert bit-exact (`rtol=0`
     is safe here because both sides are computed in the same CI job — this
     is not a cross-platform golden; see the CI-matrix lesson from the
     integration branch).
   - *Swap correctness:* `rt_run_surface(cache, brdf₂)` vs a fresh monolithic
     `rt_run` with `brdf₂` swapped into the model — bit-exact per platform.
   - *Cache-width guard:* Lambertian-sized cache + Cox-Munk replay throws
     `ArgumentError`; Canopy and RRS throw.
   - *m_max widening:* Lambertian model + `target_brdfs = [CoxMunk]` produces
     a cache whose Cox-Munk replay matches monolithic Cox-Munk.
5. **Docs:** new subsection in `docs/src/pages/concepts/05_surfaces.md`
   ("Sweeping surfaces over a cached atmosphere") + API reference entries +
   memory-footprint formula in the `AtmosphereRTCache` docstring.

**Explicitly out of PR 1:** `scenario_sweep.jl` + `netcdf_writer.jl` (they
drag in the GCHP training-set machinery; follow-up PR), linearized (Jacobian)
variant, Raman, multi-band.

---

## 3. PR 2: analytic Lambertian closure (O(1) per albedo, exact)

The gchp-io prototype replays a full `interaction!` per surface. For
Lambertian surfaces this is unnecessary: the surface reflection operator is
**rank-one**, so the final adding step collapses to a scalar rational
function of albedo — the classical `ρ_TOA(ρ_s) = ρ_atm + T↓ρ_sT↑/(1−ρ_s·S̄)`
closure already used *empirically* (two-run fit) in
`test/benchmarks/modtran_equivalent_fields.jl` and
`plans/MODTRAN_equivalent_equations.md`. Here we get the coefficients exactly,
from one cached atmosphere, with no calibration runs.

### 3.1 Derivation (matches `interaction.jl` and `lambertian_surface.jl` conventions exactly)

From `create_surface_layer!` (lambertian_surface.jl:303), the m=0 surface
reflection is

```
r⁻⁺ = 2a · u vᵀ          u = Stokes-I slot indicator (length NquadN)
                         v = D u,  D = Diagonal(qp_μN .* wt_μN)
```

and the reflected direct beam (lambertian_surface.jl:175-182) is
`j₀⁻ = 2a·μ₀·F₀ᴵ·e^(−τ_surf/μ₀) · u`. The final adding step
(`interaction_helper!`, interaction.jl:234) updates

```
J₀⁻ ← J₀⁻ + T⁻⁻ · (E − r⁻⁺R⁺⁻)⁻¹ · (r⁻⁺J₀⁺ + j₀⁻)
```

With `r⁻⁺ = 2a·uvᵀ`, Sherman–Morrison gives
`(E − 2a·u(vᵀR⁺⁻))⁻¹ u = u / (1 − a·S̄)` where

```
S̄    = 2 · vᵀ R⁺⁻ u                      spherical albedo of the atmosphere
                                          viewed from below (per spectral pt)
E_dw  = vᵀ J₀⁺ + μ₀ F₀ᴵ e^(−τ_surf/μ₀)    total downward flux at BOA,
                                          diffuse + direct (per spectral pt)
t_up  = T⁻⁻ u                             upward transmission of an isotropic
                                          BOA source (NquadN × nSpec, full
                                          Stokes — picks up polarization)
```

so the m=0 TOA source at every output stream is **exactly**

```
J₀⁻(a) = J₀⁻_atm + [ 2a · E_dw / (1 − a·S̄) ] · t_up          (†)
```

All m>0 moments are untouched by a Lambertian surface (`_zero_surface_layer!`
⇒ the final interaction is the identity on `J₀⁻`/`R⁻⁺`). Hence the
post-processed product is

```
R_SFI(a)[iVZA, :, s] = R_black[iVZA, :, s]
                     + w₀ · 2a·E_dw[s]/(1 − a·S̄[s]) · t_up[view slots, s]
```

with `w₀ = 0.5/π` (the m=0 Fourier weight) and `R_black` the a=0 result.
Identical algebra gives the BOA/downwelling analogue via
`(E − R⁺⁻r⁻⁺)⁻¹` and `r_dn = R⁺⁻u`, and the `bhr_uw`/`bhr_dw` diagnostics.

Everything (†) needs is in the **m=0 slice** of the `AtmosphereRTCache`; the
per-albedo evaluation is a broadcast — no matrix algebra at all.

### 3.2 API sketch

```julia
struct LambertianClosure{FT}
    R_black :: Array{FT,3}   # (n_vza, pol_n, nSpec) — a=0 replay, all m summed
    t_up    :: Array{FT,3}   # (n_vza, pol_n, nSpec) — T⁻⁻u at the view slots
    S̄       :: Vector{FT}    # (nSpec,)
    E_dw    :: Vector{FT}    # (nSpec,)
    w₀      :: FT
end

lambertian_closure(cache::AtmosphereRTCache) -> LambertianClosure
    # R_black via one rt_run_surface(cache, LambertianSurfaceScalar(0))
    # (reuses the tested replay path); S̄, E_dw, t_up from the cached m=0
    # R⁺⁻ / J₀⁺ / T⁻⁻ blocks — three matvecs + two dot products per spec pt.

(c::LambertianClosure)(a)          # a scalar or length-nSpec vector
    -> R_SFI  # (n_vza, pol_n, nSpec)

albedo_jacobian(c, a)   # exact ∂R/∂a = w₀·2·E_dw/(1−a·S̄)² · t_up
invert_albedo(c, R_meas; i_vza)    # closed form:
    # x = (R_meas − R_black)/(w₀·t_upᴵ);  a = x / (2·E_dw + S̄·x)
```

The retrieval inversion is exact (per spectral point, Stokes-I row) — no
iteration, no linearized RT run. This is the piece neither the gchp-io replay
nor the `NSurf` Jacobian path provides.

**Tests:** closure vs `rt_run_surface(cache, LambertianSurfaceScalar(a))` for
a ∈ {0, 0.05, 0.3, 0.8} — *not* bit-exact (Sherman–Morrison vs batched LU),
assert `rtol ≈ 50·eps(FT)` scaled; closure with a spectral vector vs
`LambertianSurfaceSpectrum`; `invert_albedo ∘ evaluate ≈ identity`;
`albedo_jacobian` vs the `rt_run_lin` NSurf Jacobian at matching a.

---

## 4. Memory: the Lambertian slim cache

For Lambertian-only targets the full per-m snapshot is overkill:

- m = 0: keep all six blocks (needed by the closure and the a=0 replay).
- m > 0: the surface layer is zero ⇒ replay is the identity on `J₀⁻`;
  `interaction_hdrf!` contributes nothing. Only `J₀⁻` (`NquadN × 1 × nSpec`)
  is needed to re-run postprocessing.

Footprint drops from `(m_max+1)·4·NquadN²·nSpec` to
`4·NquadN²·nSpec + m_max·NquadN·nSpec` — for the NquadN=60, nSpec=10⁴,
m_max=20 example: 24 GB → ~1.2 GB. Implement as
`rt_run_atmosphere(...; cache_mode = :auto)` choosing `:slim` when every
target BRDF has `component_m_max == 0`, `:full` otherwise. (An orthogonal
`device = :cpu` offload flag for `:full` caches on GPU is a cheap follow-up.)

---

## 5. Non-Lambertian BRDFs (RPV, RossLi, Cox-Munk, Canopy)

**The replay (§2) *is* the general mechanism** — it already handles RPV,
RossLi and Cox-Munk on the branch, because every BRDF enters through the same
single per-m `create_surface_layer!` + `interaction!`. Per surface you pay
`(m_max+1)` surface builds + final interactions (`O(NquadN³·nSpec)` each) —
still ~`Nz`× cheaper than a full run. No closure exists in general: a BRDF's
`r⁻⁺` is full-rank at every m ≤ m_max, so `(E − r⁻⁺R⁺⁻)⁻¹` needs a real
solve.

Structure that *can* be exploited, in decreasing payoff order:

1. **Amplitude-scaled families** — RPV's `ρ₀` is a pure amplitude
   (`r⁻⁺(ρ₀) = ρ₀·K_m` at fixed `k, Θ, ρ_c`). Then per m and spectral point,
   one eigendecomposition `K_m R⁺⁻ = WΛW⁻¹` turns every subsequent amplitude
   into `(E − ρ₀·WΛW⁻¹)⁻¹ = W·Diag(1/(1−ρ₀λᵢ))·W⁻¹` — `O(N²)` per ρ₀ instead
   of `O(N³)`. Breakeven: an eigendecomposition costs ~10–25 LU solves, so
   this pays only for sweeps of ≳25 amplitudes per scene. Verdict: document,
   don't build until a consumer (e.g. ρ₀ retrieval over imagers) demands it.
2. **Kernel-linear families** — RossLi is linear in `(fiso, fvol, fgeo)`:
   precompute the three kernel matrices per m once, surface *build* per
   parameter set becomes a 3-term axpy. The final interaction is still a full
   solve (the response is rational, not linear, in the f's). Since
   `create_surface_layer!` is cheap for RossLi anyway, savings are marginal —
   skip. (Its `fiso` term is exactly the rank-one Lambertian part; a mixed
   Woodbury update is possible but not worth the complexity.)
3. **Cox-Munk** — no useful parametric structure (wind speed enters through
   the facet distribution nonlinearly), and its surface *build* (facet
   quadrature per m) plus the TMS correction are the per-surface cost. The
   replay already amortizes what can be amortized. One extra option: cache the
   built surface `AddedLayer` per (BRDF, m) when sweeping *atmospheres* under
   a fixed ocean surface — the mirror image of this PR's sweep; note for the
   BatchContext follow-up.
4. **Canopy** — `create_surface_layer!` re-runs the canopy's own internal
   adding (canopy layers + soil) per call, so the atmosphere split helps but
   the canopy solve repeats per parameter set. Two-level split (cache the
   canopy-above-soil composite, Sherman–Morrison the *soil* albedo underneath
   it — the soil is Lambertian, so §3 applies verbatim one level down) is the
   natural follow-up and would also replace the current finite-difference
   soil Jacobian. Deferred; lift the canopy-preload rejection first, as the
   branch notes.

---

## 6. Follow-ups (explicitly deferred)

- `scenario_sweep.jl` + NetCDF training-set writer port (compose with
  `BatchContext` so the per-SZA rebuild reuses cached HITRAN/Mie).
- Linearized variant: the same snapshot/replay applied to `rt_run_lin`
  (surface Jacobians without atmosphere re-accumulation).
- Raman (`ie*` 4-D snapshot blocks), multi-band caches, canopy two-level split.
- Full-cache GPU→CPU offload flag.

## 7. Suggested PR sequencing

| PR | content | size |
|---|---|---|
| 1 | rt_run hooks + `AtmosphereRTCache` + `rt_run_split.jl` + tests + docs | ~450 lines + tests |
| 2 | `LambertianClosure` (§3) + slim cache (§4) + tests | ~250 lines + tests |
| 3 | scenario sweep / training-set writer port, BatchContext composition | follow-up |
