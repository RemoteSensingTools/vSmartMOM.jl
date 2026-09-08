# Source quality and scientific readability audit — 2026-09-08

The largest maintainability gains are at the boundaries between configuration,
prepared optics, mutable workspaces, and the numerical kernels. The existing
scientific decomposition is sound and worth preserving: optical-property
algebra, elemental/doubling/adding stages, polarization and scattering-interface
dispatch, named Jacobian layouts, and package-owned backend extensions already
use Julia effectively. A rewrite of the scientific kernels is not recommended.

This is a review and refactoring backlog, not an implementation change or a
new certification of the science. Suniti Sanghavi's formulations, scientific
linearization and Raman methodology remain the scientific foundation.

## Scope and evidence

Reviewed the release worktree, `integration/release-candidate`, starting from
`41f0a053`, while the pressure-tangent and source-validation work proceeded in
parallel. Read `AGENTS.md`, `CLAUDE.md` and `SESSION_HANDOFF.md`. Numeric anchors
below describe the reviewed snapshot; use the accompanying symbol names when
concurrent changes move lines.

The initial inventory contained 166 Julia files / 55,036 lines under `src`,
4 / 680 under `ext`, and 160 / 43,269 under `test`. These are raw file counts,
including comments, experimental scripts and large reference arrays; they are
not executable-line counts. A source-validation file was added concurrently.

Broad text and declaration scans covered those trees. Representative deeper
reads covered CoreRT model/layer types, drivers, optical mixing and Jacobian
propagation, construction/update/profile routines; CPU/CUDA/Metal algebra;
Scattering types and NAI2 forward/linearized/batched paths; standalone absorption
and HITRAN parsing; Raman types/helpers; IO and aerosol readers; SolarModel;
StandaloneSS orchestration; Lambertian/RPV surfaces; and Aqua/JET/type-stability
and test orchestration. This was not a line-by-line review of every file,
independent paper-equation verification, security audit, or a hardware benchmark.
No large test suite or GPU computation was run by this reviewer.

Evidence labels:

- **Reproduced:** a small CPU Julia 1.12.6 probe established the behavior.
- **Inspected:** the relevant implementation establishes the contract issue;
  the proposed failure/performance scenario was not run.
- **Opportunity:** a bounded improvement requiring measurement or scientific
  acceptance before deciding whether it is worthwhile.

The pressure derivative and source support omissions from the handoff are being
addressed separately in this session. They should not be counted again as new
findings in this report.

## Follow-up implementation status

The bounded v2.2 repairs identified here were subsequently implemented in the
same uncommitted release-candidate worktree:

- `BatchContext.update_model!` now prepares profiles, sources and every optical
  field in context-owned trial buffers and commits only after all bands succeed.
  A two-band late-failure regression verifies complete old-scene preservation,
  safe retry and parity with a fresh build for both unreduced and one-layer
  profiles (**44/44**).
- The standalone HITRAN reader now validates fixed-width records and required
  numbers with file/line/field diagnostics, preserves explicitly optional blank
  fields and decodes the traditional `0`/`A`/`B` isotopologue notation. All
  `test_Absorption.jl` test sets pass.
- Batched solve/inverse ownership, aliasing and return contracts are documented
  across CPU/CUDA implementations. A new CPU mutation-contract set passes
  **7/7**; the complete batched-kernel file passes **46/46**. Actual CUDA
  pivot/info workspace reuse remains a measured implementation opportunity.
- `RTModel.propertynames`, tolerance-aware `GreekCoefs.isapprox`, and the
  obsolete `noRS` Cabannes overload are repaired and covered (**6/6**).
- Generic analytic-surface construction and Fourier integration now live in
  `Surfaces/analytic_surface.jl`, rather than being hidden in the RPV-specific
  file. RPV and Ross-Li smoke checks preserve behavior.
- Standalone atmospheric-profile readers accept an explicit floating type;
  defaults remain `Float64` and focused I/O validation passes **42/42**.
- The public BLAS-thread documentation now states that this model-selected knob
  changes persistent process-wide state and must be consistent across concurrent
  solves. The compatibility behavior itself is unchanged.

Items 4, the larger portion of 5, 7's Raman state redesign, 8, 9's executor API,
and the remaining support boundaries in 10 are intentionally not folded into a
release cleanup. They require profiling, hardware, scientific acceptance or a
separately reviewable API design. See
[the follow-up record](code_quality_followup_2026-09-08.md) for rationale and
validation scope.

## Ranked work

| Rank | Bounded change | Evidence and benefit | Relative scope |
|---|---|---|---|
| 1 | Give `update_model!` a clear failure-state contract | Reproduced: failed update clears absorption but leaves a callable model; aliasing also affects input/remembered state | Medium |
| 2 | Reject malformed required HITRAN numbers | Reproduced: invalid numerical text becomes a physical zero | Small |
| 3 | Specify batched algebra mutation and workspace semantics | Inspected plus CPU probe: mutation varies with backend/batch size; advertised CUDA workspace reuse is absent | Small contract change; medium implementation |
| 4 | Preserve concrete types at hot boundaries | Inspected type erasure; speed benefit unmeasured | Medium, staged |
| 5 | Share prepared physics across construction/update paths | Inspected duplication; improves forward/tangent consistency and scientific review | Medium per stage |
| 6 | Complete standard Julia interfaces | Reproduced/inspected: tolerance-aware `isapprox`, discoverable properties | Small |
| 7 | Make Raman state and legacy paths explicit | Reproduced broken orphan overload; active and obsolete formulas are hard to distinguish | Small cleanup, medium state redesign |
| 8 | Turn quality diagnostics into actionable regression evidence | Inspected: advisory JET counts and float-type tests do not establish hot-path inference | Small/medium |
| 9 | Separate process settings from per-scene configuration | Inspected: BLAS configuration and optional scratch are process-global | Medium |
| 10 | Make units, axes, support limits and current equations visible at boundaries | Inspected naming/docs friction; preserves scientific interpretation | Small per module |

### 1. Failed batch updates can leave an inconsistent context

`update_model!` in `src/CoreRT/tools/update_model.jl:543` overwrites the current
profile, then geometry and source buffers, then Rayleigh/aerosol optics. At
`:665` it clears live absorption arrays before computing replacement values.
CIA and MT_CKD files are still loaded at `:688` and `:698`. The remembered
`ctx.current_T`, `current_q`, `current_p_half`, and `current_vmr` assignments
occur at `:712`, but delaying those assignments does not isolate aliased arrays.

A small two-case CPU probe **reproduced** the failure. A synthetic legacy LUT
covers 200–300 K; updating a valid 270/280 K profile to 350/360 K raises
`BoundsError` after live absorption has been cleared. `rt_run(ctx.model)` still
succeeds with the inconsistent state:

- Unreduced profile: `ctx.current_T === model.profile.T` and
  `params.T === model.profile.T`. All three become `[350,360]` despite the
  failed update. Delayed assignment does not preserve remembered/input state.
- Reduced-to-one-layer profile: the arrays do not alias. Live temperature
  becomes `[356.6151...]` while remembered input stays `[270,280]`.
- Both cases start with nonzero absorption and end with all-zero absorption.
  Relative radiance changes are 7.5534% and 7.5575% in these synthetic fixtures,
  not production bias estimates.

The [confirmed reproduction](pressure_source_followup_2026-09-08/batch_update_failure_probe.jl)
passed its assertions (log outside the repo:
`/home/cfranken/vsmartmom-release-evidence/quality-probes/update_failure_state_confirmed.log`).
The first probe's assumption that remembered state always remained unchanged
failed because of aliasing; that earlier script/log is retained separately.
A caller catching the update error had no reliable indication of what state
remained safe to use. The follow-up status above records the transactional fix;
this section retains the original reproduction and design rationale.

Choose an explicit contract. Prefer preparing new optics/profile/geometry in
reusable trial buffers and copying/swapping them into the model only after
successful preparation. If that memory cost is unacceptable, mark a context
invalid after any post-mutation failure and require a complete rebuild before
another solve. Merely wrapping the existing body in `try` does not restore the
original arrays.

Acceptance: induce a failure in the second absorption band or continuum loader;
verify either complete preservation of the prior scene or a deterministic
invalid-context rejection. Then retry an ordinary partial update and compare
with a fresh model. Preserve existing successful-update parity tests.

### 2. Standalone HITRAN parsing turns invalid data into plausible values

`read_hitran` at `src/Absorption/read_hitran.jl:63` uses
`something(tryparse(...), varTypes[i](0))` for every numeric field. Consequently,
malformed wavenumber, strength, molecular IDs and broadening values are all
silently replaced with zero. With the default lower bound of zero, an invalid
wavenumber can survive filtering. This is a scientific data-integrity issue,
not simply a slow parser.

The same routine builds heterogeneous rows (`rows = []`, `:57`), allocates a
one-element wrapper for each append (`:70`), then transposes all rows into
columns (`:75`). Strongly typed column builders would make the fixed-width
format easier to validate and reduce temporary storage. `HitranTable` itself
already stores concrete typed columns.

Implement a small field parser with filename/line/field diagnostics; distinguish
legitimate blank optional fields from malformed required fields, and validate
record length before slicing. Keep format-specific missing-value rules explicit.
Only then simplify row assembly with typed columns and `push!`.

Acceptance: one valid synthetic fixed-width record, malformed required numeric
fields, allowed blank optional fields, short records, and an empty selection.
Preserve valid-file numerical values. Scope matters: this concerns the package's
standalone `Absorption.read_hitran` API; the current production RT LBL pipeline
uses AtmosphericAbsorption's separate parser and is not established to share
this issue.

### 3. Batched linear algebra needs one documented ownership contract

`batch_inv!` in `src/CoreRT/tools/cpu_batched.jl:32` computes each inverse with
`A[:,:,i] \ I`, preserving `A` (confirmed by the CPU probe). CUDA's general
method at `ext/gpu_batched_cuda.jl:99` calls `getrf_strided_batched!(A, ...)` for
multiple spectral slices, overwriting `A` with LU data. Its singleton branch
uses a nonmutating solve and preserves `A`. `batch_solve!` also mutates `A` on
the CPU through `qr!` (`cpu_batched.jl:25`). The short docstrings only describe
filling `X`; Julia's `!` marks mutation somewhere, not which inputs are scratch.

In addition, the CUDA `batch_inv!(X,A,ws::RTWorkspace)` docstring at `:115`
claims preallocated pivot/info reuse, but the implementation at `:127` allocates
through `getrf_strided_batched!` and never reads `ws`.

Document which operands may be destroyed, legal aliasing, return values,
singleton/broadcast behavior and singular-system policy in the package-owned
API. A practical low-allocation contract may permit destruction of `A` on all
backends without requiring every implementation to destroy it. If preservation
is needed, provide an explicit scratch-taking preserving wrapper. Reuse the
provided workspace or correct its claim; do not infer allocation savings from
an unused argument.

Acceptance: tiny nonsingular and singular matrices, batch sizes 1 and >1,
output/input alias rules, and CPU/CUDA/Metal parity where hardware is available.
Measure allocations specifically for workspace overloads. No new GPU behavior
was measured in this review.

### 4. Recover concrete types where numerical work actually starts

There is a distinction between accepting an abstract argument (usually good
Julia API design) and storing values in abstractly typed fields (which loses
inference information). Several important carriers do the latter:

- `RTModel` in `src/CoreRT/types.jl:1522` has only architecture/float parameters;
  `solver::SolverConfig{FT}`, `atmosphere::Atmosphere{FT}` and `optics::Optics{FT}`
  erase their remaining parameters. `sources::AbstractSource` does the same.
- `AddedLayer` / `CompositeLayer` at `types.jl:338` / `:297` and solar columns
  at `:263` store arrays as `AbstractArray{FT,N}` and source slots as `NamedTuple`.
- `GreekCoefs` at `src/Scattering/types.jl:462` omits concrete storage and rank
  from all six array fields.
- `MInvariantCache` at
  `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:108` stores
  nested `Vector`, `Vector{Vector}` and similar partially specified containers.
- `BatchContext` at `src/CoreRT/tools/update_model.jl:113` stores an unparameterized
  model, unparameterized vectors and heterogeneous absorption-model lists.
- Standalone `HitranModel` / `InterpolationModel` at
  `src/Absorption/types.jl:169` / `:194` similarly erase model/interpolator types.

The `RTModel` field types were confirmed with `fieldtypes`. This establishes
loss of static information; it does **not** establish a particular runtime
penalty, because downstream function boundaries can recover specialization and
large BLAS calls may dominate.

Start with one measured boundary. An internal layer carrier can use
`CompositeLayer{FT,M,J,S}` with concrete operator storage `M`, source storage `J`
and slot tuple `S`. Keep scientific field names `R⁻⁺`, `T⁺⁺`, `J₀⁺`, etc.
Alternatively pass unpacked arrays/configuration into a concrete per-moment
function. Keep public constructors ergonomic with inferred parameters and
preserve existing partially applied `CompositeLayer{FT}` method signatures.
Do not force diffuse matrices, source vectors, views and CUDA pointer metadata
to share one storage type.

Existing good examples are `SolverConfig{FT,PT,QT}`, `SourceSet{S<:Tuple}`,
`JacobianPropagationWorkspace`, and the `_nai2_bulk_optics_lin` function boundary
(`src/Scattering/compute_NAI2_lin.jl:38`). The latter explicitly isolates dynamic
size-distribution setup from scalar arithmetic. Follow that pattern.

Acceptance: `@code_warntype`/targeted JET on one concrete call, warm allocations,
CPU time, GPU launch time if relevant, and compilation latency for representative
Float32/Float64, I/IQU and source configurations. Expand only when gains justify
specialization and constructor complexity. Metadata dictionaries, a locked
spectroscopy cache with a typed return, and temporary heterogeneous construction
buffers are not automatically defects.

### 5. Put shared physics preparation behind small named stages

`model_from_parameters` has a long elastic implementation starting at
`src/CoreRT/tools/model_from_parameters.jl:255`, another VRS constructor at
`:746`, and a separate linearized builder. Absorption-model selection, q-driven
H₂O, CIA/continuum accumulation and Rayleigh preparation recur across these
builders and `update_model!`. For example the absorption blocks at `:350–426`
and `:830–890` repeat nearly the same logic; the batch updater has its own
version. The comments in `update_model.jl:579` already record a numerical parity
bug caused by previously different mean-temperature handling.

The forward `_rt_run_column` (`src/CoreRT/rt_run.jl:335`) also interleaves
validation, global settings, source preparation, buffer allocation, Fourier
iteration, layer composition, surface coupling and output packaging across
roughly 500 lines. These are different scientific questions and different
specialization boundaries.

Refactor one stage at a time: prepare band absorption models; evaluate band
absorption; prepare molecular/Rayleigh state; build one Fourier moment; package
requested observer outputs. Return named records for their physical quantities
and hold scratch separately. Reuse the forward preparation in the tangent path
before attaching its derivatives. The current pressure work belongs at that
upstream boundary; do not replace the hand-differentiated RT operator chain
with full-solver AD.

Move generic surface-layer scaffolding out of `rpv_surface.jl:60` to a visibly
shared surface file; a scientist implementing a new BRDF should not need to
open RPV to find the universal builder. The shared Lambertian
`surface_albedo` dispatch in `lambertian_surface.jl:61` is a good model.

Avoid one universal mega-driver controlled by more booleans, and avoid macros
that hide the adding equations. Preserve product ordering, exact finite-δ
formulas, D-matrix symmetry, scattering-depth doubling and backend precision
policies. Test each extracted stage against the existing forward/rebuild and
Jacobian fixtures; do not combine file moves with scientific formula changes.

### 6. Small standard-interface fixes would immediately help scientists

`RTModel.getproperty` at `src/CoreRT/types.jl:1584` supports `model.τ_abs`,
`model.profile` and other documented aliases, but no matching `propertynames`
method is defined for `RTModel`. Introspection and `hasproperty` therefore do
not expose the same interface. `PlannedRTModelLin.propertynames`
(`src/CoreRT/types_lin.jl:261`) already demonstrates the intended pattern.
Define the alias list once and include it with real fields in `propertynames`.
Check every advertised alias with `hasproperty` and direct access.

`Base.isapprox(::GreekCoefs,::GreekCoefs)` at
`src/Scattering/types.jl:478` does not accept `rtol`, `atol`, or other standard
keywords; custom scientific tolerances fail before comparing coefficients.
Forward keywords to each field comparison, and use `all` with a function or
generator instead of allocating an intermediate Boolean vector. Define whether
tolerances apply per coefficient array (the existing structural semantics) or
to a combined norm, and keep that decision explicit. This extension is
package-owned-type dispatch, so it is not type piracy.

Do not turn this into a wholesale naming pass. Unicode maps naturally to the
scientific papers. Preserve it, including the existing `⊠` algebra. Short
semantic names for scratch roles and shape comments are more valuable than
romanizing symbols.

### 7. Raman's legacy state obscures what is live and supported

The method `compute_ϖ_Cabannes(RS_type::noRS,depol,λ₀)` at
`src/Inelastic/inelastic_helper.jl:71` assigns scalar `1.0` to
`noRS.ϖ_Cabannes::Vector{FT}`. A CPU call reproduced the resulting `MethodError`.
No call to this particular overload was found in production code or tests;
do not interpret this as a demonstrated failure of ordinary elastic RT.
Remove/deprecate it if obsolete, or repair its documented scalar/vector and
mutation semantics with a small regression test. Its non-`!` name also hides
its intended mutation.

More generally `RRS` at `src/Inelastic/types.jl:18` combines molecular constants,
Greek coefficients, mutable per-band selections, solar/SIF buffers and computed
phase arrays. `fscattRayl`, `bandSpecLim` and related fields are untyped; the
runtime role of the object changes during the solve. Large commented-out
formula blocks in `inelastic_helper.jl:76` onward and
`stellar_inelastic_helper.jl:156` onward make it hard to tell which convention
is active. Some VS docstrings still say `struct RRS` (`types.jl:46`, `:75`).

Separate fixed molecular data from a prepared spectral redistribution state and
per-solve scratch; explicitly label Earth N₂/O₂ versus stellar H₂ routing.
Replace historical alternative implementations with links to dated development
notes, retaining the scientific rationale, attribution and version provenance.
This is a readability proposal, not authorization to reconcile differing Raman
formulas without scientific review. Keep deliberate Float64 molecular constants:
their subnormal-range rationale is documented at the top of `types.jl`.

### 8. Quality tests should distinguish scalar type, inference and science

`test/test_type_stability.jl` mostly checks that outputs retain Float32/Float64;
its `@inferred` assertions focus on small helpers. Those are useful checks but
do not establish inference through `RTModel → layer setup → interaction`.
`test/test_jet.jl:23` retains a historical baseline of 232 findings; advisory
mode succeeds regardless of count, and only the first 20 reports are printed.
The handoff's 285 is previous evidence, not a fresh run here. Broad counts can
both vary with Julia/JET and hide a new serious finding when another disappears.

Keep broad JET advisory, but save structured findings grouped by symbol and
reason, triage real undefined paths separately from macro noise, and gate a
small set of concrete representative calls in a pinned environment. Use narrow
inference checks for boundaries modified in item 4. Separate float-preservation,
allocation/performance and physical-conservation tests in names and reports.
Keep genuine operator identities, finite differences and full-rebuild parity;
do not substitute tests that simply mirror the implementation.

`test/runtests.jl:17` still changes process working directory, and many fixtures
rely on shared top-level imports/state. Gradually use `joinpath(@__DIR__,...)`
and local fixture helpers so focused files run independently, then remove the
cwd requirement once all dependencies are converted. Continue obeying the
current `cd test/` rule in the meantime. The runner already restores the cwd
with `finally`, and Aqua already checks ambiguities; preserve those improvements.

### 9. Per-model execution knobs can mutate process-global state

`_rt_run_column` at `src/CoreRT/rt_run.jl:345` calls
`BLAS.set_num_threads(model.numerics.blas_threads)` for a configured model. This
is a process-wide setting that persists after the solve and can affect unrelated
linear algebra and concurrent scenes. It is not a private property of that
model. Moving it into `try/finally` alone would still race between concurrent
calls. Prefer an explicit application/executor-level setting and document its
scope; retain a compatibility path if necessary.

The optional `_CUDA_BATCH_INV_RUNNERS` dictionary at
`ext/gpu_batched_cuda.jl:18` owns shared scratch keyed only by `(FT,n,batch)`.
Its comment already states that identical-size concurrent calls are unsupported;
the key does not encode CUDA device/stream identity. Keep this benchmark path
opt-in and move reusable ownership into per-solve workspaces before considering
it a general production optimization. Existing locked read-only HITRAN caching
(`src/CoreRT/tools/hitran_cache.jl:1`) and explicit Mie workspace ownership checks
(`src/Scattering/compute_NAI2_nodes_batched.jl:407`) illustrate different,
appropriate treatments of immutable data and mutable scratch.

No race was reproduced here. This is an inspected ownership limitation with a
clear refactoring direction, not a claim that normal single-task runs fail.

### 10. Use boundary documentation to resolve axes, units and support

There are at least three important layouts: RT diffuse operators
`(direction×Stokes,direction×Stokes,spectrum)`; model absorption
`(spectrum,layer)`; and StandaloneSS contributor optics `(layer,spectrum)`
(`src/StandaloneSS/solver.jl:103`). Conversion functions should state both sides
explicitly. A scientist should not infer the meaning of the third or fourth
dimension from temporary variable names. For each tangent carrier, document the
parameter axis and link it to `ParameterLayout`/`ParameterKey` rather than
repeating offsets.

The internal constant `nm_per_m` actually equals `1e7`, i.e. nm per cm;
`src/Inelastic/InelasticScattering.jl:38` now documents the misleading name.
Prefer correctly named private conversion helpers and aliases over globally
changing scientific constants. Keep intentional Float64 Mie arithmetic and
explicit host/device conversions; float genericity does not mean eliminating
numerically necessary wider intermediates.

Other bounded support clarifications:

- `_rt_run_column` selects only the first surface for a multi-band call
  (`src/CoreRT/rt_run.jl:374`) and emits an info message. Reject incompatible
  per-band surfaces or require an explicit common-surface mode rather than
  letting a scientifically meaningful distinction depend on reading a log.
- `_to_aa_arch` has CPU/CUDA but no Metal mapping
  (`src/CoreRT/CoreRT.jl:140`). This is a known end-to-end absorption limitation;
  successful Metal Mie or matrix kernels do not establish whole-model support.
  Validate capabilities at preparation and report the unsupported component.
- Legacy `reduce_profile_binavg` uses unweighted VMR means
  (`src/CoreRT/tools/atmo_prof.jl:577`), as its TODO acknowledges. A dry-column
  conserving variant needs a scientific acceptance test, not a cosmetic change.
  The current default reduction is a separate interpolation path at `:471`;
  this finding does not establish a defect in that default.
- `read_atmos_profile_dict` deliberately converts all profiles to Float64
  (`src/IO/AtmosProfile.jl:47`), while other readers accept `FT`. Clarify this
  IO contract or add an explicit `FT` keyword; do not call it a kernel precision
  defect.

Comments should lead with current quantities, equations, assumptions and
mutation ordering. Lengthy completed "Phase A/B/C" narratives can move to
linked development notes. Retain short explanations of nonobvious numerical
choices and provenance. For example, immediately above an adding product, a
comment stating "reuse the pre-update downward transmission" helps a scientist
more than a history of the scratch-allocation optimization. Keep a readable
reference equation next to fused implementations rather than hiding the
operation in generic metaprogramming.

## Small probe record

One CPU process loaded the package and reported:

- `fieldtypes(RTModel{CPU,Float64})` includes `SolverConfig{Float64}`,
  `Atmosphere{Float64}`, `Optics{Float64}`, `Vector{<:AbstractSurfaceType}` and
  `AbstractSource`, confirming the declaration-level type erasure.
- No `propertynames(::RTModel, ...)` specialization; only the planned-linearized
  wrapper has one.
- `compute_ϖ_Cabannes(noRS(),0.03,760.0)` throws a Float64-to-Vector conversion
  `MethodError`.
- The CPU `batch_inv!` result leaves a copied input matrix unchanged.

A second CPU process constructed a synthetic 160-character HITRAN record with
`BROKEN` in its wavenumber field and valid remaining required numeric fields.
`read_hitran` accepted the record and returned `νᵢ == [0.0]`. The same process
confirmed `isapprox(get_greek_rayleigh(0.03), get_greek_rayleigh(0.03);
rtol=1e-6)` throws `MethodError`. These are small contract probes, not RT
validation. No real spectroscopy was downloaded for either probe; package
precompilation occurred because files were changing concurrently.

## Suggested implementation sequence

The selected pressure/source work and numerical checks are now complete; see
[the follow-up record](pressure_source_followup_2026-09-08.md).
The first three contract repairs plus the small Julia interface fixes are now
implemented. Next take one measured concrete-carrier refactor or one shared
absorption/Rayleigh preparation extraction as an independent change. Only after
those results should Raman-state redesign, larger driver separation or broad
module/file reorganization proceed.

Do not combine these into a formatting-only mega-commit. Do not promise speedups
from annotations alone, automatically replace explicit loops with broadcasts,
require every object to be immutable, remove every `Any`, or rewrite all branch
logic as types. Julia's advantage here is preserving scientific structure while
specializing the computation where it matters.
