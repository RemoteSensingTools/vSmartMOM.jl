# Proposal: batched caller-node Mie seam + reusable GPU workspace (v2)

Status: DRAFT v2 (2026-08-16) — supersedes v1 after external review. Verdict
adopted from that review: **proceed with batching, but not via a single
globally-padded batch**. v1's global-max design broke its own bitwise
contract and carried ~8× memory amplification on realistic TOMAS columns.
v2 is built around **exact-`nmax` grouping + immutable prepared geometry +
single-owner workspaces**; a fully ragged one-launch solver is explicitly
deferred until profiling proves the grouped design launch-bound.

Motivation (unchanged from v1):

1. Serial per-work-item GPU submission loses to threaded CPU on the real
   GCHPIO dust column: A100 1.51 s vs 16-thread CPU 1.44 s at 0.765 µm
   (1.30 vs 0.87 s at 1.61 µm) — launch + allocation overhead dominates.
2. Codex review (2026-08-16): "fragmented, per-work-item Mie calls … serial
   host submission on GPU and fresh device allocations."

Interaction with the (landed) band-interpolation feature
(`GCHPIO compute_aerosol_optics_band`, Google-RT `82638c9`): interpolation
reduces the number of λ evaluations per band to ~5, so batching's payoff is
per-λ-node column throughput and global multi-column scenes — benchmark
AFTER interpolation, not instead of it.

## Why v1 was wrong (review findings, kept as design constraints)

1. **Global angular grid breaks bitwise equivalence.** The single-ensemble
   API sizes its angular grid per ensemble: `n_mu = 2·n_max_global − 1`
   (`compute_NAI2_nodes.jl:393`). A batch-max grid changes array lengths,
   quadrature nodes, projection ops, and rounding — it cannot be bitwise
   equal to the single calls. → **Group ensembles by exactly equal
   `n_max_global`**: within a group, the angular grid and reduction order
   are identical to the single-ensemble path by construction.
2. **Kernel 1 is not ensemble-blind.** It takes scalar `m_re, m_im`
   (`gpu_mie_kernels.jl:163`). → Batched kernel 1 needs a
   `node_to_ensemble` index + per-ensemble RI arrays (per-node RI arrays
   are the forward-compatible superset the smooth-RI follow-up will want —
   choose per-ensemble now, but keep the index layout compatible).
3. **Global padding is materially expensive.** Real dust column inventory
   at 0.55 µm: 882 ensembles, 20 distinct `nmax` values (2–81), 84,672
   nodes at nquad=96 / 903,168 at nquad=1024. Padded-to-max Float32
   (aₙ/bₙ + four phase arrays): 0.33 GB vs 0.041 GB ragged (nquad=96);
   3.50 GB vs 0.44 GB (nquad=1024) — excluding Greek tables, reductions,
   scratch. Exact-`nmax` grouping turns 882 submissions into ~20 group
   launches with no 8× memory penalty.

## API

Two-object design separating λ-independent from λ-dependent state (review
finding 4 — avoids repacking/uploading up to ~9·10⁵ nodes per spectral
node):

```julia
# λ-INDEPENDENT, immutable, device-resident after construction:
geometry = prepare_mie_node_geometry(radii, weights, offsets;
                                     architecture)
#  - concatenated radii + weights (weights normalized per ensemble in
#    Float64 — P2 lesson — then narrowed), offsets (length E+1),
#    node_to_ensemble index; uploaded once.

# λ-DEPENDENT per call:
compute_aerosol_optical_properties_nodes_batched(
    geometry, λ, n_real::Vector, n_imag::Vector,   # length-E per-ensemble RI
    pol, trunc;
    l_max = nothing,
    precision_policy = default_mie_precision_policy(architecture, FT),
    workspace = nothing,       # ::Union{Nothing, MieBatchWorkspace}
) -> Vector{AerosolOptics}

compute_aerosol_extinction_nodes_batched(geometry, λ, n_real, n_imag; ...)
```

FT discipline (review addition): the prepared geometry is locked to an
explicit `FT` at construction; per-call λ / RI arrays of any other element
type are REJECTED with `ArgumentError`, never silently promoted — otherwise
the bitwise single-call contract is ambiguous (which FT would the single
call have used?).

Group gather (review addition): a group's ensembles are generally
NONCONTIGUOUS in the concatenated node arrays. Kernels therefore address
nodes through a device-side per-group index (ensemble permutation +
per-ensemble node offsets built at grouping time, i.e. a gather), while
each ensemble's own nodes keep their ORIGINAL intra-ensemble order — the
reduction loop walks offsets[e]:offsets[e+1]-1 exactly as the single call
does, so the summation order (and hence bitwise equality) is unaffected by
where the ensemble sits in the group permutation.

Per-call, internally:
- size parameters + per-ensemble `n_max_global` are computed for THIS λ,
  ensembles are grouped by exact `n_max_global` (λ-dependent grouping —
  group membership may change between λs; that is correct and cheap),
- one kernel-1 launch per group over that group's concatenated nodes
  (per-ensemble RI via `node_to_ensemble`),
- segmented weighted reduction per ensemble, each segment reducing in the
  SAME order as the single-ensemble kernel (one thread per
  (ensemble, output element), serial Neumaier node loop — segmentation
  changes loop bounds only, so bitwise equality falls out),
- one Greek-projection launch per group over (ensemble, l) work items,
- results reassembled in the ORIGINAL ensemble order.

Contract: for every ensemble, batched result `==` the single-ensemble call
with the same policy/backend — now achievable **by construction** because
each group reproduces the single-call angular grid and reduction order
exactly. Truncation (δBGE then `l_max` cap — P1 ordering) stays per-ensemble
on the host.

All three `MiePrecisionPolicy`s supported; `NativeFloat32`'s
zero-Float64-device-array invariant extends to geometry and workspace
arrays; the Metal guard fires identically.

## MieBatchWorkspace

Explicit semantics (review finding 5):

- Keyed on **backend instance/device** (not just backend type), FT, RA,
  polarization ONLY if it changes raw scratch layout (it does not for
  kernel-1/reduction scratch — so no pol key there).
- Capacities tracked separately: grouped-node count, coefficient count
  (Σ nodes·nmax per largest group), node-angle scratch, Greek table size.
  Grow-only per array; request > capacity reallocates that array.
- **Single-owner, NOT thread-safe**: one workspace per concurrent CPU task
  or GPU stream. Constructor documents this; a debug-mode owner check
  (task-id stamp) is cheap and worth it.
- `nothing` ⇒ allocate-and-drop (current behavior; zero API break).

CPU path: **keep the existing ensemble-level `Threads.@threads` path** as
primary — a single sequential "CPU batch" would regress the measured
16-thread performance. The batched API on CPU = threaded loop over
ensembles reusing per-task workspaces; grouping is a GPU concern.

## Acceptance gates

1. Bitwise `==` batched-vs-single per ensemble, all three policies, KA-CPU
   + real CUDA (Metal wired, hardware-gated); full + extinction variants.
2. Batching invariance: splitting one batch arbitrarily (including across
   group boundaries) changes no ensemble's result (bitwise).
3. Workspace-reuse invariance: warm == cold, and reuse across λs (regroup
   between calls) == fresh.
4. Original-order reassembly property-tested (shuffled ensemble sizes).
5. Performance, measured AFTER band-interpolation is in the consumer:
   real GCHPIO column, one batched call per λ-node; go/no-go bar
   **GPU ≥ 2× over 16-thread CPU** on A100; also benchmark consumer CUDA
   and Metal (the NativeFloat32 targets). If grouped batching can't clear
   2×, the seam stays CPU-primary and the workspace still pays on CPU.
6. Full suite + docs gates; julia-reviewer on API + kernels.

## Sequencing (adopted from the review)

1. Polarized validation of NativeFloat32 + band interpolation (phase
   matrices, Stokes I/Q) — IN FLIGHT on both repos; hard prerequisite.
2. Band-API diagnostic/physical-bound fixes — IN FLIGHT (GCHPIO).
3. This proposal (v2) — done.
4. Implement exact-`nmax` grouped batching (geometry object, grouped
   kernel-1 with per-ensemble RI, segmented order-preserving reductions).
5. Device-aware single-owner workspaces.
6. Benchmarks on A100 / consumer CUDA / Metal, after interpolation.
7. Ragged one-launch solver ONLY if profiling shows the grouped version
   launch-bound (~20 launches/λ-node is unlikely to be).

## Consumer wiring (GCHPIO, after merge + repin)

`unified_seam.jl`'s collect phase already builds the work-item list →
`prepare_mie_node_geometry` once per column (λ-independent, like
`ReconstructionPlan` — same staleness rule: geometry contains aerosol node
state only, nothing met-derived); per λ-node one batched call with the
per-ensemble RI vector for that λ. Hybrid exact-τ leg uses the batched
extinction variant. `_LayerBundle` fitting is orthogonal and unchanged.
