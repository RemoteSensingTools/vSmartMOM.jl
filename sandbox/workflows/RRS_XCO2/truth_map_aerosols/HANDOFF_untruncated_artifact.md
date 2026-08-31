# Handoff: the "untruncated" angular artifact is a Float32 defect

Last updated: 2026-08-19
Working branch: `suniti_multi_sensor`
Supersedes the "Principal diagnostic finding" section of
[`../WORK_SUMMARY.md`](../WORK_SUMMARY.md).

## Executive summary

The untruncated angular references plotted in
`truncation_vs_untruncated_IQU.png` (111 streams) and
`truncation_vs_untruncated_50stream_IQ.png` (50 streams) are **not**
references. They are the least accurate runs in either figure. Every feature
that distinguishes them from the δBGE-truncated ladder — the local maxima near
signed VZA = ∓20°, the minimum near −45°, the spikes near ∓55°, and the
VZA-to-VZA jaggedness — is numerical error, not physics.

The cause is **Float32 precision, amplified by stream count**. It has nothing
to do with the Legendre truncation, the phase function, or the Fourier order.

The δBGE ladder (`l_trunc` = 8, 12, 16, 24, 32) is the correct answer, and it
already contains the local maximum at signed VZA = −30° that was expected on
physical grounds (exact backscatter, Θ = 180°, at SZA = 30° / Δφ = 180°).

## The decisive measurement

TOA `I` at 757.001655 nm, SZA = 30°, VZA = +20°, Δφ = 0°, scene 9.
Reference is the most converged Float64 run (111 streams).

| run | streams | Legendre | float | I | error |
|---|---:|---:|---|---:|---:|
| `big64_n111_l32` | 111 | δBGE 32 | F64 | 5.475489 | — |
| `big64_n50_l32` | 50 | δBGE 32 | F64 | 5.475486 | −0.0001% |
| `ref64_n17_l32` | 17 | δBGE 32 | F64 | 5.475457 | −0.0006% |
| `unt64_n50_l99` | 50 | **untruncated, l=221** | F64 | 5.475476 | −0.0003% |
| `ref32_n17_l32` | 17 | δBGE 32 | F32 | 5.477041 | +0.028% |
| `big32_n50_l32` | 50 | δBGE 32 | F32 | 5.454085 | **−0.391%** |
| `unt32_n50_l99` | 50 | **untruncated, l=221** | F32 | 5.454073 | −0.391% |
| `big32_n111_l32` | 111 | δBGE 32 | F32 | 5.749065 | **+5.00%** |

Read as a 2×2 at 50 streams:

|  | δBGE(32), l=32 | untruncated, l=221 | spread |
|---|---:|---:|---:|
| **Float64** | 5.475486 | 5.475476 | 2e-6 |
| **Float32** | 5.454085 | 5.454073 | 2e-6 |
| **spread** | 0.39% | 0.39% | |

**Truncation level changes the answer by 2 parts per million. Precision changes
it by 0.39%.** That is the entire finding in one table.

The pre-existing benchmark runs fall exactly on the Float32 row: the 50-stream
untruncated dataset gives 5.4541 (vs 5.454085 here, 3e-6 apart) and the
111-stream untruncated run gives 5.74905 (vs 5.749065 here, 3e-6 apart). Both
are reproduced to six significant figures using a **δBGE(32)-truncated** phase
function, which is why the phase function can be excluded as a factor.

Confirmed independently at VZA = +45° (reference = F64/50 streams = 5.625426):
F64/17str −0.0007%, F32/17str −0.097%, F32/50str +0.580%, and the existing
untruncated runs at +0.581% (50 streams) and −2.925% (111 streams).

## Supporting evidence

**1. The untruncated series adds no physics.** From `untruncated_beta_757nm.dat`,
the mixture Greek β_l:

| l | 16 | 32 | 48 | 99 | 150 | ≥180 |
|---|---:|---:|---:|---:|---:|---:|
| β_l | 5.3e-2 | 3.5e-3 | 4.1e-4 | 2.9e-6 | 2.5e-8 | ~1e-13 |

Past l ≈ 50 there is nothing; past l ≈ 160 the coefficients *are* the Float64
Mie noise floor. These are 0.1–0.2 µm accumulation-mode particles at 757 nm:
the mixture phase function peaks at only f₁₁(0°) = 10.6, with f₁₁(2.5°) = 10.45
and just 0.5% of scattered energy inside 2.5°. There is no forward spike. The
`m_max = 220` ceiling came from the r_max = 10 µm Mie integration limit, i.e.
from a size-distribution tail carrying negligible optical weight.

**2. The δBGE ladder is converged.** `l_trunc` = 8, 12, 16, 24, 32 agree with
one another to ≤0.1% at every VZA. Five independent truncation levels landing
on the same curve is convergence.

**3. The δBGE curve reproduces the single-scattering geometry**, including the
local maximum at signed VZA = −30°, which is exactly Θ = 180°. The mixture
phase function rises toward backscatter (f₁₁ = 0.166 at Θ = 130° → 0.280 at
180°) and Rayleigh peaks there (1.48), so the feature is physical. The
111-stream Float32 error is ~5% at that angle and buries a feature whose true
amplitude is ~1%.

**4. The error is confined to m = 0.** Comparing the two moment-resolved
datasets at identical m ≤ 99 support: `Σ_{m>0}` agrees between 50 and 111
streams to 0.002 absolute (0.04% of I), while m = 0 differs by up to 5.1%. The
Fourier series is fine; the azimuthally averaged component is not. The
near-horizon streams dominate that component, which is where the precision
loss is worst.

## Mechanism

`get_dtau_ndoubl` ([`src/CoreRT/CoreKernel/rt_kernel.jl:266-287`](../../src/CoreRT/CoreKernel/rt_kernel.jl#L266-L287)):

```julia
threshold  = dτ_max_threshold                    # default 0.001
floor_val  = dτ_min_floor                        # default 1024·eps(FT)
μ_min      = minimum(qp_μ[wt_μ .> eps(FT)])      # positive-weight nodes only
dτ_max     = max(floor_val, minimum([maximum(τ .* ϖ), threshold * μ_min]))
_, ndoubl  = doubling_number(dτ_max, maximum(τ .* ϖ))
dτ         = τ ./ 2^ndoubl
```

The accuracy criterion `dτ_max = 0.001·μ_min` keeps `dτ/μ ≤ 0.001` for every
stream. `μ_min` shrinks like 1/N²:

| nstreams | μ_min | dτ_max = 0.001·μ_min | approx. ndoubl for τϖ ≈ 0.1 |
|---:|---:|---:|---:|
| 17 | 4.7e-3 | 4.7e-6 | ~15 |
| 50 | 5.7e-4 | 5.7e-7 | ~18 |
| 111 | 1.2e-4 | 1.2e-7 | ~20 |

So raising the stream count forces exponentially thinner elemental layers and
a longer doubling ladder. The elemental r/t entries scale as `dτ/μ`, so in
Float32 a seed layer of `dτ ≈ 1.2e-7` carries very few significant digits, and
that noise compounds through the doubling ladder. Float64 has ~9 more decimal
digits of headroom and is unaffected.

### Measured: the controlling variable is `dτ_max`, not `nstreams`

`doubling_sweep.dat` varies `dτ_max_threshold` at fixed δBGE(32) physics.
Error is against the converged Float64 answer, 5.475489:

| dτ_max | dτ/μ_min | streams | threshold | I (F32) | error |
|---:|---:|---:|---:|---:|---:|
| 1.16e-8 | 1e-4 | 111 | 1e-4 | 8.047872 | **+46.98%** |
| 1.16e-7 | 1e-3 | 111 | 1e-3 | 5.749065 | +5.00% |
| 5.67e-7 | 1e-3 | 50 | 1e-3 | 5.454085 | −0.39% |
| 1.16e-6 | 1e-2 | 111 | 1e-2 | 5.482080 | +0.12% |
| 4.71e-6 | 1e-3 | 17 | 1e-3 | 5.477041 | +0.028% |
| 5.67e-6 | 1e-2 | 50 | 1e-2 | 5.475245 | −0.0045% |
| 1.16e-5 | 1e-1 | 111 | 1e-1 | 5.475493 | **+0.0001%** |
| 1.22e-4 | 1.05 | 111 | floor | 5.474097 | −0.025% |

Float64 at the same settings: 111 streams +0.00000%, 50 streams −0.00007%,
17 streams −0.00059%.

Three conclusions:

1. **The error is monotonic in `dτ_max` and collapses onto it regardless of
   stream count.** 50 streams at `threshold=1e-2` (dτ_max = 5.7e-6) is *more*
   accurate than 17 streams at the default (dτ_max = 4.7e-6). `nstreams` only
   matters because the criterion ties `dτ_max` to `μ_min`.
2. **There is no "too thick" regime to speak of.** At `dτ/μ_min = 1.05` the
   error is only −0.025%. This vindicates the exact finite-δ elemental
   formulation: layer thickness is not an accuracy problem for it, only
   precision is. The failure is entirely on the thin side.
3. **The floor would have prevented the artifact outright.** Had
   `dτ_min_floor` been the intended `1024·eps(Float32) = 1.22e-4`, the
   111-stream Float32 run would have returned −0.025% instead of +5.00%.

For Float32, keep `dτ_max ≳ 1e-5` (≈ 100·eps(F32)). The existing
`1024·eps(FT)` default satisfies this with margin.

### Correction to an earlier hypothesis

An initial reading of the code suggested that the Float32 floor was *active*
and making the elemental layer too **thick** (`dτ/μ_min` → 1.05 at 111
streams). The direction was backwards. Two probe runs with the floor
explicitly lowered returned bit-for-bit identical results, because the floor
was already inactive at 2.27e-13; and the sweep above shows `dτ/μ_min = 1.05`
is in fact the *most* accurate Float32 setting tested at 111 streams. The
`dτ/μ_min` = 0.026 / 0.215 / 1.05 table from that hypothesis does not describe
what ran — the true ratio was 0.001 throughout.

## Bugs found

**Bug 1 (root cause) — precision-dependent defaults are cast, not recomputed.**
[`_convert_numerics`](../../src/CoreRT/tools/model_from_parameters.jl#L611-L617)
re-types a user `RTNumericalParameters` to the model float type by casting each
stored value. `dτ_min_floor`'s default is `1024·eps(FT)`, so a Float64 YAML
parse stores 2.27e-13 and casting to Float32 keeps 2.27e-13 instead of the
intended 1.22e-4. Any workflow that parses a Float64 config and *then* sets
`float_type = Float32` programmatically loses the guard entirely. That is
exactly what `generate_truth_map.jl:105` (`params.float_type = FT`, `FT =
Float32` at line 27) does, so it affects every script in `RRS_XCO2/scripts/`.
Per the sweep above, restoring the intended floor alone takes the 111-stream
Float32 error from +5.00% to −0.025%, so this is the direct proximate cause of
the artifact, not merely a contributing factor.

Verified directly:

```
YAML float_type          = Float64
params dτ_min_floor      = 2.2737367544323206e-13
after float_type=Float32:
  converted dτ_min_floor = 2.2737368e-13     <-- stale
  intended F32 default   = 0.00012207031
```

**Bug 2 — the `_ss` path never got the `RTNumericalParameters` refactor.**
[`rt_kernel_ss.jl:352`](../../src/CoreRT/CoreKernel/rt_kernel_ss.jl#L352) still
carries the pre-refactor form, flagged in-source with `# SUNITI, check?`:

```julia
dτ_max = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])
```

It ignores both knobs (which the enclosing function accepts and forwards to
other callees) and uses `minimum(qp_μ)` over **all** nodes rather than
positive-weight ones — so a grazing user VZA drives `dτ` down on that path,
the exact failure the main path was fixed to avoid. Not on the truth-map path,
but latent.

**Bug 3 — `gfct` scaling inconsistency in the directional variant.**
[`rt_kernel.jl:307-308`](../../src/CoreRT/CoreKernel/rt_kernel.jl#L307-L308)
computes `dτ_max` from `maximum(gfct * τ .* ϖ)` but passes `maximum(τ .* ϖ)`,
without `gfct`, as the target to `doubling_number`. These disagree whenever
`G[iμ₀] ≠ 1`.

## Consequences for the experiment

- **Retire both untruncated references.** They are ~150× less accurate than the
  runs they were built to validate. The 12.3 CPU-hours (2.899 h + 9.419 h)
  spent generating them measured a precision defect, not angular convergence.
- **The `greek_beta_cutoff = 1e-5` conclusion is moot.** β₉₉ = 2.9e-6 for the
  mixture; the radiance is converged by l ≈ 16. Demanding l = 99 forced
  nstreams ≥ 50, which is what pushed the calculation into the failing regime.
  For this aerosol set a cutoff near 1e-3 (l ≈ 40) is already generous.
- **The production truth map is safe as configured.** `nstreams: 5`,
  `truncation: auto` sits at the low-stream end where Float32 is accurate to
  ~0.03%, and δBGE(8) already matches the converged ladder to 0.1%.
- **Rule of thumb pending a proper guard:** Float32 is reliable up to roughly
  20 streams per hemisphere for this problem; beyond that use Float64.

## Recommended fixes

1. Make `_convert_numerics` recompute precision-dependent defaults instead of
   casting them — or store the floor as a dimensionless multiple of `eps(FT)`
   so re-typing is automatically correct.
2. Warn or error when the resolved `dτ_max` is below a precision-safe bound for
   the working float type, so an unsupportable `nstreams`/`float_type`
   combination announces itself rather than returning a plausible wrong answer.
   Any such guard must be honest about the trade: at very high stream counts a
   Float32 run may have no safe setting at all, in which case the right answer
   is to refuse rather than clamp.
3. Bring `rt_kernel_ss.jl` onto `get_dtau_ndoubl` and fix the `gfct`
   inconsistency.
4. Add a regression test asserting Float32/Float64 agreement across a stream
   sweep at fixed physics — the property that would have caught this.

## Artifacts

New:
- `scripts/diagnose_untruncated_precision.jl` — the 2×2×2 probe
  (streams × Legendre support × precision). Output:
  `truth_map_aerosols/untruncated_precision_probe.dat`.
- `scripts/diagnose_doubling_sweep.jl` — sweeps `dτ_max_threshold` over four
  decades at fixed δBGE(32). Output: `truth_map_aerosols/doubling_sweep.dat`.
  Result: monotonic in `dτ_max`; the failure is entirely on the thin side.
- `scripts/plot_untruncated_artifact_diagnosis.py` →
  `truth_map_aerosols/untruncated_artifact_diagnosis.png`.
  Panel (d) plots the measured Float32 error against `dτ_max` from
  `doubling_sweep.dat`.

Reinterpreted (data still valid, conclusions changed):
- `untruncated_stokes_by_m_nstreams50_mmax99.{jld2,dat}`
- `untruncated_stokes_by_m_nstreams111_mmax99.{jld2,dat}`
- `truncation_iqU.dat`, `untruncated_convergence_cpu*.jld2`

## Open items

1. Fix Bugs 1–3 and add a Float32/Float64 stream-sweep regression test.
2. Decide whether the production truth map moves to Float64. On the evidence
   above `nstreams: 5` in Float32 is accurate to ~0.03%, so this is a cost
   question, not a correctness one.
3. Optional: if an untruncated cross-check is ever wanted again, run it in
   Float64 (`unt64_n50_l99` agreed with δBGE(32) to 3 ppm) — or in Float32
   with `numerics.dτ_max_threshold: 1e-1`, which cost 332 s versus 731 s for
   the Float64 equivalent and matched it to 1e-6.
