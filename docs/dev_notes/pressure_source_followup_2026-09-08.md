# Pressure absorption, source contracts, and quality review — 2026-09-08

Completed the two implementation priorities from `SESSION_HANDOFF.md` in
`vSmartMOM-release`, on `integration/release-candidate`, based on `41f0a053`.
Changes remain local and uncommitted. No push, tag, registration, deployment,
or external message was sent. Scientific formulation and linearization credit
remains with Suniti Sanghavi; integration and validation are supporting work.

## Pressure-dependent absorption

`src/CoreRT/tools/absorption_pressure.jl` supplies the missing cross-section
response at the upstream optical-property boundary. The coordinate is **hPa
at the bottom interface of the final model grid**, with other interfaces,
layer temperature, humidity, and VMR fixed. Thus `dp_full[end]/dp_surf = 1/2`:

```text
dτ/dp_surf = σ * VMR * dN_dry/dp_surf
          + (dσ/dp_full) * N_dry * VMR / 2.
```

The existing molecular-column contribution is retained. Fixed gases, variable
gases and q-driven H2O now accumulate the cross-section term through the same
helper. Disabling H2O abundance Jacobians does not disable its pressure tangent.
CIA/MT_CKD retain their separate density/column pressure factors.

- **AtmosphericAbsorption line-by-line models:** reuse the actual prepared
  lines/windows; seed a dimensionless pressure scale in the collisional
  widths, shifts, velocity-changing frequency and line mixing. ForwardDiff
  differentiates the actual scalar line-profile/CPF implementation inside a
  portable KernelAbstractions kernel. Divide by layer pressure for the hPa
  derivative. Doppler width, strength and correlation remain fixed at fixed
  temperature/composition. The RT operator receives plain arrays.
- **Native ABSCO:** use the exact active pressure-interval slope, evaluating
  both pressure-node spectra through the same public interpolation API. This
  preserves each pressure node's own temperature grid, broadener interpolation
  and spectral resampling. No raw-table difference at mismatched temperatures.
- **Regular AtmosphericAbsorption LUTs:** share the pressure-interval slope;
  a scalar adapter also aligns their linearized path with batched forward use.
- **Legacy BSpline LUTs:** differentiate the actual interpolant with a dual
  pressure coordinate, preserving its out-of-band spectral zeros.
- **Optional CPU Erfcx reference CPF:** a package-owned strategy supplies the
  exact complex scalar chain rule, since SpecialFunctions does not accept
  `Complex{Dual}` directly. No methods were added to foreign numeric types.

The AA preparation/profile dependency is isolated in this file and validated
against installed AtmosphericAbsorption 0.1.2. Future upstream changes to line
preparation should be checked against these finite-difference tests.

At an interior native LUT pressure knot the right interval is selected; at or
outside either clamped pressure endpoint the derivative is zero. Hard line
wing cutoffs and interpolation kinks are not differentiable at their changes.
This pressure column does not differentiate profile regridding, observer-grid
insertion, or a coordinate that moves other pressure interfaces. CIA/MT_CKD
abundance tangents remain incomplete. No new retrieval-bias estimate is claimed.

## Source requests

`Sources/validation.jl` shares checks across forward, linearized, SS, streams,
TOA, routed interior-observer and atmosphere/surface split paths. Cache replay
also checks the new surface against its prepared fluorescence sources.

- An explicit solar SZA must agree with model geometry at model precision.
- Multiple solar beams are rejected, including a `BlackbodySource` combined
  with another beam. Same-geometry irradiances can be summed into one beam.
- Nonzero prescribed or retrievable `SurfaceSIF` requires a supporting surface.
  A prescribed zero remains a harmless placeholder. `supports_surface_sif`
  is an extensible dispatch trait, currently true for the four Lambertian
  variants that implement injection and tangents.
- `rt_run_ss_exact` now honors stored/overridden unpolarized solar irradiance
  and `NoSource`. Its atmosphere-only reference rejects nonzero/retrievable
  SIF, thermal emission, and polarized incident illumination.

The unexported legacy `rt_run_test_ms` helper retains unit-beam transport;
its docstring directs callers to the public `rt_run` wrapper that applies
source validation and normalization. Source composition/preparation utilities
remain separate from the solver's supported physical combinations.

## Validation

All numerical gates ran on a frozen implementation. Tests ran from `test/`.
Julia 1.12.6 used four Julia threads and one BLAS thread; CUDA used the idle
A100 on device 0. Device 1 and other users' workflows were left untouched.
Logs live in `/home/cfranken/vsmartmom-release-evidence/`.

| Check | Final result | Log |
|---|---|---|
| Complete CPU suite | **12,671 passed, 16 skipped, zero failures/errors**, exit 0, 24m34.1s | `cpu-priorities.log` |
| Complete GPU runner | **20,301/20,301**, exit 0, 6m59.5s | `gpu-priorities.log` |
| Focused pressure tests, CPU | **195 passed**, one CUDA skip, exit 0 | `pressure-final.log` |
| Same pressure file with CUDA | **203/203**, including eight CUDA derivative comparisons | `pressure-final-cuda.log` |
| Source contracts | **47/47**, exit 0 | `source-validation-final-v3.log` |
| Real HITRAN/ABSCO optical depth, RT and convolution | **16/16** additional assertions after the synthetic suite | `pressure-real.log` |
| Julia 1.10.12 pressure/source/Aqua checks | **252 passed**, one CUDA skip, exit 0 | `priorities-julia110-final.log` |
| Strict Documenter/VitePress and generated assets | Exit 0, deployment disabled | `docs-priorities.log` |
| Executable documentation contracts | **20/20**, exit 0 | `docs-contracts-priorities.log` |

The complete CPU run supersedes the September 7 record's missing complete
post-fix rerun; that historical record is not rewritten. JET still reports
**285 advisory findings versus historical baseline 232**, unchanged from the
previous candidate audit. Passing the suite does not mean JET is clean.
Hosted platform CI, Metal hardware and publication gates remain outstanding.
The ten-tutorial execution record remains the earlier documentation audit;
this follow-up rebuilt docs and reran the small executable contracts.

The regression begins with independently rebuilt forward models. Before the
fix, Float64 synthetic optical-depth pressure errors were 31.3–32.2% for LBL
and 16.1–16.2% for ABSCO (`pressure-initial.log`). Tests cover Float32/64,
fixed/variable/mixed/H2O cases, regular/native/legacy tables, clamping and
nonconstant pressure slopes, additive helpers, inactive lines, pressure shifts,
line mixing, default/reference CPF, and Hartmann–Tran. Radiance and fixed
convolution checks are separate and cover embedded and external solar modes.

Float32 finite differences use zero-shift synthetic lines: a tiny shift of a
6200 cm⁻¹ line can be smaller than one Float32 ULP of the stored line centre.
Nonzero pressure shifts are tested in Float64. Earlier Float32 shifted-centre
FD failures were retained, not relabeled as passing; CPU/CUDA derivative parity
does not remove this finite-precision issue in the underlying forward line
centre arithmetic.

The [real-spectroscopy probe](pressure_source_followup_2026-09-08/real_absorption_probe.jl)
uses two layers, five CO2 spectral samples, fixed final-grid temperature and
composition, `h=0.02 hPa`, and normalized fixed detector weights. Both fixed and
variable gas construction paths are exercised. It selects the genuine HITRAN
constructor branch with empty LUT storage and reads a narrow native ABSCO slab
from a separately supplied file; no spectroscopy was copied into the repo.

| Real CO2 case | Relative L2 optical-depth derivative error | Radiance derivative error | Convolved derivative error |
|---|---:|---:|---:|
| HITRAN LBL, fixed/variable | 1.485e-8 | 4.287e-8 | 3.986e-8 |
| Native ABSCO, fixed/variable | 3.013e-12 | 2.954e-8 | 2.908e-8 |

These are fixture-specific derivative comparisons, not full retrieval or
instrument-shift/width validation. The existing retrieved-state precision
investigations remain separate.

The first implementation tests caught an invalid fully specified Dual
constructor and a prepared-source dispatch omission; both were corrected.
Test-fixture aliasing and tuple/named-output assumptions were also corrected.
An initial Julia 1.10 attempt reused the Julia 1.12 manifest and failed in
PrecompileTools before loading the package; a separately resolved 1.10 test
environment passed. Earlier failed logs are retained under distinct filenames.

## Independent quality review and next work

The [source-quality audit](source_quality_audit_2026-09-08.md) covers scientific
readability and idiomatic Julia across source, extensions and tests. It ranks
ten bounded work areas and separates reproduced defects, inspected contracts
and unmeasured optimization opportunities. A second review of this session's
pressure/source diff found no additional pressure sign/normalization blocker;
it prompted the regular-LUT adapter and exact-SS source corrections.

The first separate follow-up should be **failed batch-update state integrity**.
A [synthetic reproduction](pressure_source_followup_2026-09-08/batch_update_failure_probe.jl)
confirmed that a late LUT `BoundsError` can clear absorption yet leave a callable
model. With unreduced profiles, input and remembered temperatures alias the
live profile and change too; with a reduced profile, remembered state remains
old while the live model changes. The resulting radiance changes about 7.55%
in that fixture. This pre-existing defect is documented, not fixed here.

Other high-value changes are strict required-field HITRAN parsing, explicit
backend mutation/workspace contracts, discoverable Julia property/comparison
interfaces, measured concrete-type preservation, and shared physics preparation.
The audit recommends retaining the scientific kernel decomposition rather than
a broad rewrite or blanket formatting/type-annotation pass.
