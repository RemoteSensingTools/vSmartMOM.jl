# Priority corrections after independent review

Checkout: `vSmartMOM-release`, branch `integration/release-candidate`.
These changes build on the pressure/source and code-quality work and are
included in implementation commit **36e7d1e0**. They do not change the MOM
equations or publish a release. See `SESSION_HANDOFF.md` for the subsequent
user-authorized commit/push checkpoint.

## Corrections

- `BatchContext` uses `copy_parameters(...; share_luts=true)` once for both
  construction and its cached configuration. Profiles, lists and configuration
  are owned; large LUT payloads remain shared read-only. Regression checks
  assert coefficient-array identity, not a size-dependent memory estimate.
- Solar-angle validation compares both values at their narrower floating-point
  precision. Tests cover both Float32/Float64 directions, real mismatches and
  non-finite input.
- `requires_pressure_jacobians(flavor)` follows the existing upstream-work
  trait pattern. It defaults to true; fixed-pressure retrievals can override it
  or pass `compute_pressure_jacobians=false`. Disabled pressure arrays are
  `nothing`, and neither their storage nor their derivatives are constructed.
  Plan construction and the optical-cache boundary reject requests for missing
  pressure derivatives, including raw full-layout solves. A forward-only
  absorber fixture tests selected gas/surface Jacobians through both local and
  physical propagation bases against a supported absorber.
- Both aerosol updaters stage all bands in existing scratch storage. The Mie
  updater also stages optical objects and derives candidate Fourier bounds
  before committing. Regression tests inject a second-band staging failure,
  check live optics/loading/reference extinction/bounds, then verify a valid
  retry against fresh construction. Existing successful-update tests remain
  part of the validation.
- AtmosphericAbsorption is temporarily constrained to `=0.1.2` because the
  line-pressure implementation still uses its internal prepared-line interface.
  The tested Julia 1.12 manifest records tree
  `3f9b48e459c504dda42304551fe4743200f5ad26`. A pin does not identify modified
  same-version development trees; a public upstream derivative API remains the
  long-term solution. Local manifests are untracked, not package lockfiles.
- Migration notes now describe beam/SIF errors, strict HITRAN parsing and
  isotopologue IDs. Surface extension docs identify the inherited analytic
  builder. The new pressure trait has an API entry and extension guidance.

## Validation

| Check | Result |
|---|---|
| Source requests and transactional updates, Julia 1.12 | 128 passed |
| Upstream availability, copies, selective Jacobians, existing batch updates, Aqua | 207 passed; one existing degenerate-Mie skip |
| Complete A100 GPU runner | **20,323/20,323**, exit 0, 7m05.4s |
| Source requests and transactional updates, clean Julia 1.10 resolution | 128 passed |
| Pressure regression, Julia 1.12 with CUDA | **213/213**, exit 0, 2m57.5s |
| Pressure regression, Julia 1.10 CPU | **205 passed**, one CUDA skip, exit 0, 2m03.4s |
| Strict Documenter/VitePress build | Passed, exit 0 |
| `git diff --check` | Passed |

The grouped Julia runs each found one over-strict new assertion: forward and
linearized absorption column integration differ by ordinary rounding because
their multiplication orders differ. Only that assertion was changed to
`rtol=4eps(Float64)`; exact selected radiance/Jacobian comparisons were retained.
The pressure file passed separate reruns on both versions. The other grouped test
results above remain applicable to the unchanged production source. The full
GPU result is on the corrected production tree, not the earlier tranche.

Logs for this correction pass:

- `/tmp/vsmartmom-review-regressions-final.log`: source, transactions, pressure
  (including CUDA), upstream availability, parameter copies, selective
  Jacobians, existing batch updates and Aqua; includes the superseded strict
  assertion described above.
- `/tmp/vsmartmom-review-pressure-final.log`: corrected pressure-file rerun,
  Julia 1.12 with CUDA.
- `/tmp/vsmartmom-review-gpu.log`: complete `test/local/gpu/runtests.jl` on A100.
- `/tmp/vsmartmom-review-julia110.log`: fresh Julia 1.10 dependency resolution
  and transaction/source/pressure regressions (CUDA hidden); also includes the
  superseded strict assertion.
- `/tmp/vsmartmom-review-pressure-julia110-final.log`: corrected pressure-file
  rerun in that clean Julia 1.10 environment.
- `/tmp/vsmartmom-review-docs-final.log`: strict Documenter/VitePress build.

Tests run from `test/`, with four Julia threads and one BLAS thread. CUDA runs
use device 0. An initial test-fixture delegation error and a missing API-doc
entry were corrected before the final reruns; failed assertions are not counted
as passing evidence. Both Julia environments resolve AtmosphericAbsorption to
the 0.1.2 tree recorded above.

## Remaining release gates

Hosted CI and Metal hardware validation remain outstanding. This pass does not
rerun the entire CPU suite; the earlier 12,743-pass CPU result predates these
corrections. JET's historical 232-to-285 advisory delta remains untriaged.
Surface quadrature allocation improvements and a public absorption derivative
API remain follow-ups. No commits, push, tag, registration or publication were
performed in this correction pass.
