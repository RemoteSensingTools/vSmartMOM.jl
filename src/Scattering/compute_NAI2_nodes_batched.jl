#=
=====================================================================
Batched caller-node Mie seam (exact-nmax grouping + reusable GPU workspace).
=====================================================================
Implements proposals/batched_mie_nodes_seam.md (v2). See that file for the
full design rationale (why v1's global-max-padded batch was rejected, the
memory-amplification numbers, the review findings this design encodes).

Two-object API, separating λ-independent from λ-dependent state:
- `prepare_mie_node_geometry(radii, weights, offsets; architecture)` builds an
  immutable `MieNodeGeometry`: concatenated radii/weights (per-ensemble
  Float64-normalized weights, narrowed once), offsets, and a node→ensemble
  index. Device-resident after construction for GPU architectures.
- `compute_aerosol_optical_properties_nodes_batched(geometry, λ, n_real,
  n_imag, pol, trunc; ...)` / `compute_aerosol_extinction_nodes_batched(...)`
  take per-ensemble refractive index vectors for ONE wavelength and return a
  `Vector{AerosolOptics}` / `Vector{k}` in ORIGINAL ensemble order.

GPU grouping strategy (the core correctness argument, restated close to the
code that implements it): ensembles are grouped by EXACT `n_max_global` at
this λ (λ-dependent — group membership may change between λs; that is
expected and cheap). Within a group:
  - Groups address their (generally noncontiguous, in the FULL geometry's
    original layout) ensembles via a DEVICE-SIDE permutation
    (`_group_permutation`, uploaded once per group as a small `Int` array)
    plus per-ensemble LOCAL offsets built at grouping time: `geometry.radii[perm]`
    / `geometry.weights[perm]` (standard `AbstractArray` fancy indexing — a
    proper GPU gather kernel under CUDA.jl/Metal.jl) assemble the group's
    contiguous node arrays WITHOUT a host round-trip for the float node
    payload; each ensemble's own intra-ensemble node order is untouched (the
    permutation is built from `offsets[e]+1:offsets[e+1]` slices, in order).
  - Kernel 1 (Mie coefficients) launches ONCE over the group's concatenated
    nodes, using the `_batched` sibling kernels (`mie_coefficients_kernel_f64_batched!`/
    `_ds_batched!` in gpu_mie_kernels.jl) with a `node_to_ensemble` index into
    per-ensemble (LOCAL, within-group) refractive-index arrays. Every other
    per-node quantity a node's aₙ/bₙ depend on (`x`, `n_max_i`, the
    per-thread `nmx_i` recursion depth) is already indexed per node, so this
    substitution changes NOTHING about any single node's arithmetic — bitwise
    identical to the single-ensemble kernel launch for that node's own
    ensemble. Per-node `n_max_i` is itself computed device-side
    (`nmax_per_node_kernel!`, calling the same scalar `get_n_max`), so `x`
    never needs a host round-trip either.
  - Kernels 2+3 (amplitude/phase) and 4a (cross-sections) are the EXISTING,
    UNMODIFIED single-ensemble kernels, launched over the whole group: valid
    because grouping is by exact `n_max_global` (hence shared `n_mu`), and
    `k_wavenum` is one scalar shared by the whole λ-call regardless of
    ensemble.
  - The weighted reduction and Greek projection use NEW `segmented_*`
    sibling kernels (`segmented_weighted_sum_kernel!`,
    `segmented_size_reduction_kernel!`, `greek_coefficients_kernel_batched!`):
    one thread per (ensemble[, angle/l]), each looping over ONLY that
    ensemble's own node sub-range within the group (`local_offsets[j]+1:local_offsets[j+1]`)
    in the SAME ascending order the single-ensemble kernel would use for its
    own (un-segmented) array — segmentation changes loop BOUNDS only, so the
    sequence of (value, weight) pairs visited, and therefore the
    Neumaier-compensated sums, are bitwise identical to the single-ensemble
    path.
  - Results are reassembled into a `Vector` indexed by ORIGINAL ensemble id
    (not group-local order), via each group's `ens_ids` list.

Every per-call λ/refractive-index array is LOCKED to `geometry`'s own `FT`
(`eltype(geometry.radii_host)`, fixed at `prepare_mie_node_geometry` time):
passing `wavelength`/`n_real`/`n_imag` of any OTHER element type throws
`ArgumentError` (`_check_batched_ft_lock`) rather than silently promoting or
narrowing — a silent conversion here would silently change which kernel
precision actually runs, not just round an output for display.

CPU path: **no grouping at all** — grouping is a GPU concern (a single
"CPU batch" would regress the measured 16-thread-vs-serial-GPU-submission
performance gap this proposal exists to close). The CPU batched functions are
a `Threads.@threads` loop over ensembles, each delegating to the SAME shared
CPU core (`_nai2_bulk_optics` / `_mie_node_ab_and_C_ext!` in compute_NAI2.jl)
the single-ensemble CPU path uses — with geometry's ALREADY-normalized
per-ensemble weight slice passed straight through (no re-normalization),
which is what makes the batched CPU result bit-identical to
`compute_aerosol_optical_properties_nodes`/`compute_aerosol_extinction_nodes`
called directly on that ensemble's raw (radii, weights): both normalize the
SAME raw weights via the SAME Float64 arithmetic exactly ONCE (geometry does
it at construction time instead of at call time), so there is no
double-normalization rounding perturbation to break bitwise equality.
Internal CPU compute always runs in Float64 (`IC = Float64`, matching the
single-ensemble CPU path exactly) regardless of the CALLER's own
radii/weights element type -- but the OUTPUT still narrows to the caller's
own type (`MieNodeGeometry.output_FT`, captured from the caller's ORIGINAL
radii/weights eltype at `prepare_mie_node_geometry` time, kept separate from
any internal-compute forcing): an all-`Float32` CPU call returns `Float32`,
exactly matching what `compute_aerosol_optical_properties_nodes` itself
would return for the same inputs.

**Public bit-compatibility contract, precisely, by backend** (see
`docs/src/pages/gpu.md`'s "Batched caller-node Mie seam" section for the
full writeup and measured numbers): `Architectures.CPU()` and any portable
KA backend (KA-CPU) achieve TRUE bitwise (`==`) batched-vs-single equality,
by construction, as argued above. Real CUDA/Metal hardware do NOT reproduce
this bitwise -- only near-machine-precision -- because the batched and
single-ensemble calls launch kernels over different `ndrange` shapes, and a
GPU compiler is free to generate different (but equally valid) machine code
for logically identical, race-free source compiled for a different launch
configuration. The one invariant that DOES hold exactly (`==`), even on real
CUDA, for all three precision policies, is
`compute_aerosol_extinction_nodes_batched`'s per-ensemble result equalling
`compute_aerosol_optical_properties_nodes_batched`'s `.k` for that SAME
ensemble (the two share the identical Kernel-1 + cross-section +
segmented-weighted-sum path; extinction-only just never continues into the
launch-shape-sensitive amplitude/phase/Greek reduction) -- this is the exact
invariant GCHPIO's hybrid exact-τ mechanism rests on, and it has its own
dedicated strict real-CUDA regression test (not folded into the
tolerance-based comparisons) in `test/test_mie_nodes_batched.jl`.
=====================================================================
=#

# ============================================================================
# MieNodeGeometry
# ============================================================================

"""
    MieNodeGeometry{GFT, ARCH, RV, IV}

Immutable, λ-INDEPENDENT, device-resident (for GPU architectures) caller-node
size-distribution geometry backing [`compute_aerosol_optical_properties_nodes_batched`](@ref)
/ [`compute_aerosol_extinction_nodes_batched`](@ref). Built once via
[`prepare_mie_node_geometry`](@ref) and reused across many per-λ batched
calls (and, in the intended GCHPIO usage, across many spectral bands) without
re-uploading or re-normalizing anything.

# Fields (all internal; use the accessors / pass the geometry object around,
do not construct or mutate directly)
- `architecture`: fixed at construction; determines `GFT` and whether `radii`/
  `weights`/`node_to_ensemble` are device-resident copies (GPU) or identical
  to the host arrays (CPU — no copy is made).
- `backend`: the ACTUAL KA backend `radii`/`weights`/`node_to_ensemble` were
  built on (`nothing` for `Architectures.CPU()`) — may be an explicit
  override passed to [`prepare_mie_node_geometry`](@ref)'s `backend` keyword,
  NOT necessarily `Architectures.ka_backend(architecture)`. The top-level
  batched dispatchers (`compute_aerosol_optical_properties_nodes_batched`/
  `compute_aerosol_extinction_nodes_batched`) launch kernels on THIS stored
  backend rather than re-deriving one from `architecture`, so an overridden
  geometry stays internally consistent no matter which entry point is used.
- `output_FT`: the CALLER's ORIGINAL `radii`/`weights` promoted element type
  (`float(promote_type(eltype(radii), eltype(weights)))`, captured BEFORE any
  internal-compute forcing) — this is what the batched CPU output narrows to
  (together with `n_real`/`n_imag`/`wavelength`'s own eltypes), exactly
  mirroring [`compute_aerosol_optical_properties_nodes`](@ref)'s own output-type
  formula. Distinct from `GFT` for `Architectures.CPU()` geometries, where
  `GFT` is forced to `Float64` for internal compute regardless of the
  caller's actual `radii`/`weights` eltype (e.g. an all-`Float32` CPU call
  must still return `Float32`, not `Float64`) — for GPU architectures
  `output_FT === GFT` trivially, since `GFT` is never forced there.
- `n_ensembles`, `n_nodes`: `E` and total node count.
- `offsets::Vector{Int}`: length `E+1`, ALWAYS host-resident (needed for
  per-λ host-side grouping and the CPU threaded loop regardless of
  architecture); `offsets[1] == 0`, `offsets[end] == n_nodes`.
- `radii_host::Vector{GFT}`: host mirror of the concatenated, narrowed radii
  — needed every λ-call for the host-side per-ensemble `n_max_global`
  computation that drives GPU grouping, even when `radii` itself is
  device-resident.
- `radii::RV`, `weights::RV`: concatenated, narrowed-to-`GFT` radii and
  PRE-NORMALIZED (per ensemble, in Float64, then narrowed — see
  [`prepare_mie_node_geometry`](@ref)) weights. `RV === Vector{GFT}` (same
  object as `radii_host`/the host weight mirror) for `Architectures.CPU()`;
  a device array type (`CuArray`/`MtlArray`) for GPU architectures.
- `node_to_ensemble::IV`: length `n_nodes`, 1-based ensemble id (in ORIGINAL
  ensemble order) per node.

`GFT` (the geometry's own float type) is `Float64` for `CPU()` architectures
(matching the single-ensemble CPU path's `IC = Float64` internal-computation
convention exactly, so the pre-normalized weights need no further narrowing
or widening round-trip before being fed to the shared `_nai2_bulk_optics`
core) and `float(promote_type(eltype(radii), eltype(weights)))` for GPU
architectures (matching the single-ensemble GPU path's kernel type `FT`).
"""
struct MieNodeGeometry{GFT<:AbstractFloat, ARCH<:Architectures.AbstractArchitecture,
                        RV<:AbstractVector{GFT}, IV<:AbstractVector{Int}}
    architecture::ARCH
    backend::Any                 # the ACTUAL KA backend geometry's device arrays were built on
                                  # (`nothing` for Architectures.CPU()); reused by the top-level
                                  # dispatchers so they never re-derive a (possibly DIFFERENT)
                                  # backend from `architecture` -- see prepare_mie_node_geometry's
                                  # `backend` keyword docstring.
    output_FT::Type               # caller's ORIGINAL radii/weights promoted eltype, independent of
                                  # any internal compute-precision forcing (e.g. CPU geometries force
                                  # GFT=Float64 for internal compute, but a Float32 CALLER must still
                                  # get a Float32 OUTPUT -- see _batched_cpu_full/_batched_cpu_extinction)
    n_ensembles::Int
    n_nodes::Int
    offsets::Vector{Int}
    radii_host::Vector{GFT}
    weights_host::Vector{GFT}
    radii::RV
    weights::RV
    node_to_ensemble::IV
end

@doc raw"""
    prepare_mie_node_geometry(radii, weights, offsets;
        architecture = Architectures.CPU(), backend = nothing) -> MieNodeGeometry

Build the λ-independent, immutable [`MieNodeGeometry`](@ref) backing the
batched caller-node Mie API. Concatenates `E` ensembles' nodes (`E =
length(offsets) - 1`) into flat `radii`/`weights` arrays, normalizes each
ensemble's weights to sum to 1 (Float64, then narrowed — the same
"P2 lesson" convention as [`compute_aerosol_optical_properties_nodes_gpu`](@ref):
normalizing in Float64 before narrowing avoids a Float32 sum of
large-but-finite weights overflowing to `Inf`), and (for GPU architectures)
uploads the result once.

# Arguments
- `radii::AbstractVector`: concatenated wet particle radii for ALL
  ensembles, one per node, **[μm]**, ensemble-contiguous (ensemble `e`'s
  nodes occupy `radii[offsets[e]+1:offsets[e+1]]`).
- `weights::AbstractVector`: concatenated non-negative node weights, same
  layout as `radii`. **Need not be normalized per ensemble** — normalized
  internally, exactly as the single-ensemble API's weights are.
- `offsets::AbstractVector{<:Integer}`: length `E+1`, `offsets[1] == 0`,
  non-decreasing, `offsets[end] == length(radii)`. Every ensemble needs at
  least one node (`offsets[e+1] > offsets[e]`).

# Keyword arguments
- `architecture`: `Architectures.CPU()` (default: `GFT = Float64` always,
  matching the CPU-threaded path's `IC = Float64` convention exactly, no
  device copy made) or a GPU architecture with a registered Mie pipeline
  (`Architectures.has_gpu_mie`; `GFT = float(promote_type(eltype(radii),
  eltype(weights)))`, i.e. the caller's own choice — typically `Float32` for
  `NativeFloat32`/`DSEmulated`). On `MetalGPU()`, `radii`/`weights` must be
  `Float32` (or narrower) — `Float64` throws `ArgumentError` immediately
  (Metal has no Float64 hardware support at all; the batched device arrays
  would be illegal), mirroring the single-ensemble
  `_prepare_node_mie_gpu`'s Metal guard.
- `backend`: `nothing` (default) derives the device-array backend from
  `architecture` (`Architectures.ka_backend`) whenever `architecture` is
  NON-CPU, exactly like every other GPU entry point in this module. Pass an
  EXPLICIT KA backend (e.g. `KernelAbstractions.CPU()`) to override this — the
  main use case is testing a GPU-style geometry (`architecture` non-CPU, so
  `GFT` is the caller's Float32/Float64 choice rather than the CPU path's
  forced `Float64`) on a portable backend without real GPU hardware,
  mirroring how [`compute_aerosol_optical_properties_nodes_gpu`](@ref) takes
  an explicit `backend` decoupled from any "architecture" concept. Whichever
  backend is resolved here (the override, or the `architecture`-derived
  default) is STORED on the returned geometry (`MieNodeGeometry.backend`) and
  reused verbatim by every later call, INCLUDING the top-level
  `compute_aerosol_optical_properties_nodes_batched`/
  `compute_aerosol_extinction_nodes_batched` dispatchers — they launch
  kernels on this stored backend rather than re-deriving one from
  `architecture`, so an overridden geometry stays internally consistent no
  matter which entry point (top-level or the explicit-backend `_gpu` sibling)
  is used on it afterward. Ignored (must be left `nothing`) when
  `architecture isa Architectures.CPU`.

# Returns
[`MieNodeGeometry`](@ref), reusable across arbitrarily many
[`compute_aerosol_optical_properties_nodes_batched`](@ref) /
[`compute_aerosol_extinction_nodes_batched`](@ref) calls at different
wavelengths (grouping is recomputed fresh each call — see those functions'
docstrings — but the underlying node data is never re-normalized or
re-uploaded).
"""
function prepare_mie_node_geometry(radii::AbstractVector, weights::AbstractVector,
                                    offsets::AbstractVector{<:Integer};
                                    architecture::Architectures.AbstractArchitecture = Architectures.CPU(),
                                    backend = nothing)
    n_nodes = length(radii)
    @assert length(weights) == n_nodes "radii and weights must have the same length"
    @assert length(offsets) ≥ 2 "offsets must encode at least one ensemble (length ≥ 2)"
    E = length(offsets) - 1
    @assert offsets[1] == 0 "offsets must start at 0"
    @assert offsets[end] == n_nodes "offsets[end] ($(offsets[end])) must equal length(radii) ($n_nodes)"
    @assert issorted(offsets) "offsets must be non-decreasing"
    @assert all(w -> w ≥ 0, weights) "weights must be ≥ 0"
    @assert all(r -> r > 0, radii) "radii must be > 0"
    @assert !(architecture isa Architectures.CPU) || backend === nothing "backend override is only meaningful for a non-CPU architecture (a CPU-architecture geometry never has device arrays)"

    # Caller's ORIGINAL radii/weights promoted type, captured BEFORE any
    # internal-compute forcing below -- this is the type the batched CPU
    # output must narrow to (see MieNodeGeometry's `output_FT` field and
    # _batched_cpu_full/_batched_cpu_extinction). For GPU architectures this
    # equals GFT trivially (GFT is never forced away from the caller's choice
    # there), so `output_FT` costs nothing extra on that path.
    output_FT = float(promote_type(eltype(radii), eltype(weights)))

    GFT = architecture isa Architectures.CPU ? Float64 :
          float(promote_type(eltype(radii), eltype(weights)))

    resolved_backend = if architecture isa Architectures.CPU
        nothing
    elseif backend !== nothing
        backend
    else
        Architectures.ka_backend(architecture)
    end
    backend = resolved_backend
    if backend !== nothing && _is_metal_backend(backend) && GFT === Float64
        throw(ArgumentError(
            "Metal has no Float64 hardware support: geometry radii/weights element type " *
            "$GFT would need a Float64 device array. Build the geometry from Float32 " *
            "radii/weights on MetalGPU()."))
    end

    offsets_v = collect(Int, offsets)
    radii_host = GFT.(radii)
    weights_host = similar(radii_host)
    node_to_ensemble_host = Vector{Int}(undef, n_nodes)

    # Normalize per ensemble in Float64 THEN narrow to GFT -- same recipe,
    # same order, as the single-ensemble GPU/CPU paths' own internal
    # normalization, so a batched result stays bit-identical to a
    # single-ensemble call on the same raw (radii, weights) (see this file's
    # top-of-file docstring for the full argument).
    w64 = Float64.(weights)
    for e in 1:E
        lo = offsets_v[e] + 1
        hi = offsets_v[e + 1]
        @assert hi ≥ lo "ensemble $e has zero nodes (offsets[$e] == offsets[$(e+1)]); every ensemble needs ≥ 1 node"
        seg = view(w64, lo:hi)
        s = sum(seg)
        @assert isfinite(s) && s > 0 "ensemble $e weights must sum to a positive, finite value"
        seg ./= s
        @inbounds for i in lo:hi
            weights_host[i] = GFT(w64[i])
            node_to_ensemble_host[i] = e
        end
    end

    if architecture isa Architectures.CPU
        return MieNodeGeometry(architecture, backend, output_FT, E, n_nodes, offsets_v, radii_host, weights_host,
                                radii_host, weights_host, node_to_ensemble_host)
    else
        AT = KernelAbstractions.allocate
        radii_dev = AT(backend, GFT, n_nodes)
        weights_dev = AT(backend, GFT, n_nodes)
        n2e_dev = AT(backend, Int, n_nodes)
        KernelAbstractions.copyto!(backend, radii_dev, radii_host)
        KernelAbstractions.copyto!(backend, weights_dev, weights_host)
        KernelAbstractions.copyto!(backend, n2e_dev, node_to_ensemble_host)
        return MieNodeGeometry(architecture, backend, output_FT, E, n_nodes, offsets_v, radii_host, weights_host,
                                radii_dev, weights_dev, n2e_dev)
    end
end

# ============================================================================
# MieBatchWorkspace
# ============================================================================

@doc raw"""
    MieBatchWorkspace(backend, FT::Type, RA::Type)

Reusable, GROW-ONLY GPU scratch space for the batched caller-node Mie API,
keyed on the **backend VALUE** (`objectid(backend)`) plus `FT` (kernel
precision) and `RA` (reduction precision, see `_mie_reduction_type`) — NOT
polarization: polarization does not change the raw kernel-1/reduction scratch
layout, so there is no separate workspace per `AbstractPolarizationType`.

**`objectid(backend)` is backend-VALUE identity, not device/stream identity**:
KA backend structs (`KernelAbstractions.CPU()`, `CUDA.CUDABackend()`,
Metal's backend) are `isbits` values carrying only backend-TYPE/config
fields (e.g. `CUDABackend`'s `prefer_blocks`/`always_inline` flags) — NOT
which physical GPU or CUDA stream is currently active. Two separately
constructed `CUDA.CUDABackend()` values are `objectid`-EQUAL regardless of
`CUDA.device()`/`CUDA.stream()` (verified: same `objectid` across 2 distinct
GPUs on a multi-GPU host). This check therefore catches genuinely different
backend TYPES/configs (a KA-CPU workspace handed a CUDA launch, or vice
versa) but does **NOT** catch a workspace built while one CUDA device/stream
was active being reused while a DIFFERENT device/stream is active — that is
a real gap, not yet closed: true device/stream-aware keying (e.g.
incorporating `CUDA.device()`/`CUDA.stream()` into the key) is future work.
Single-device, single-stream usage (the common case) is unaffected.

Capacities are tracked independently PER SCRATCH ARRAY (grouped-node count,
coefficient count, node-angle scratch, Greek-table size): a request larger
than an array's current capacity reallocates JUST that array (grow-only —
never shrinks, so reuse across calls with varying group sizes never
thrashes). `nothing` (the default `workspace` keyword everywhere in this
module) means allocate-and-drop, i.e. the pre-existing zero-workspace
behavior — passing an explicit `MieBatchWorkspace` is a pure optimization,
never a correctness requirement (see the workspace-reuse-invariance test in
`test/test_mie_nodes_batched.jl`: warm == cold, and reuse across two
different wavelengths — where group membership genuinely changes — still
reproduces the fresh/cold result bit-for-bit).

**NOT thread-safe. Single-owner only**: one workspace per concurrent GPU
stream / CPU task. The FIRST use from a given task stamps `owner_task`; a
LATER use from a DIFFERENT task throws an `ArgumentError` immediately rather
than silently corrupting scratch shared across two concurrent callers. This
check is UNCONDITIONAL (every `_ws_array!` call, not gated behind
bounds-checking) — it is one field compare (plus, on first use, one field
write), cheap enough that making it a debug-only check would only trade a
negligible cost saving for the ability to silently disable safety via
`@inbounds`. Every use also verifies `objectid(backend) === backend_id` (the
backend VALUE the workspace was constructed for): this catches a workspace
built for one backend TYPE/config being handed to a launch against a
DIFFERENT one (e.g. KA-CPU vs CUDA vs Metal, or a differently-configured
`CUDABackend`). It does **NOT** distinguish CUDA devices or streams (see the
caveat above) — single-device, single-stream usage only; true device/stream
keying is future work.
"""
mutable struct MieBatchWorkspace
    backend_id::UInt64
    FT::Type
    RA::Type
    owner_task::Union{Nothing,Task}
    buffers::Dict{Symbol,Any}
end

function MieBatchWorkspace(backend, FT::Type, RA::Type)
    return MieBatchWorkspace(objectid(backend), FT, RA, nothing, Dict{Symbol,Any}())
end

@inline function _ws_owner_check!(ws::MieBatchWorkspace, backend)
    # UNCONDITIONAL (not @boundscheck-gated -- review finding 4): this is one
    # field compare plus, on first use, one field write, which is cheap enough
    # that gating it behind bounds-checking (silently disabled by an
    # `@inbounds` caller) is not worth the risk of (a) a workspace's cached
    # device arrays being reused across a DIFFERENT backend VALUE (backend
    # TYPE/config -- KA-CPU vs CUDA vs Metal, or a differently-configured
    # CUDABackend -- NOT device/stream identity: KA backend structs are
    # isbits values that do not encode which physical GPU/stream is active,
    # so this does NOT catch a workspace reused across two CUDA devices/
    # streams with an otherwise-identical CUDABackend() value; see
    # MieBatchWorkspace's docstring), or (b) two concurrent tasks/streams
    # silently corrupting scratch shared by grow-only reuse (one task's
    # buffer can be resized out from under another's in-flight kernel).
    if objectid(backend) !== ws.backend_id
        throw(ArgumentError(
            "MieBatchWorkspace was built for a different backend VALUE (objectid " *
            "mismatch): its cached device arrays cannot be reused across a different " *
            "backend type/config (e.g. KA-CPU vs CUDA vs Metal). Construct a new " *
            "MieBatchWorkspace for this backend. Note: this does NOT detect a workspace " *
            "reused across different CUDA devices/streams sharing the same CUDABackend() " *
            "value -- see MieBatchWorkspace's docstring."))
    end
    t = current_task()
    if ws.owner_task === nothing
        ws.owner_task = t
    elseif ws.owner_task !== t
        throw(ArgumentError(
            "MieBatchWorkspace is single-owner (not thread-safe): this workspace was " *
            "first used from a different task/stream. Construct one MieBatchWorkspace " *
            "per concurrent CPU task or GPU stream."))
    end
    return nothing
end

"""
    _ws_array!(workspace, backend, key::Symbol, ::Type{T}, dims) -> array

Grow-only scratch-array accessor shared by all batched GPU compute functions.
`workspace === nothing` ⇒ allocate-and-drop (fresh `KernelAbstractions.allocate`
every call, the zero-workspace default). Otherwise: reuse (reshaping, never
copying) the buffer cached under `key` if its linear capacity already covers
`prod(dims)` and its eltype matches `T`; else allocate a fresh, larger buffer
under that key (grow-only — the old, smaller buffer for `key` is dropped, not
resized in place, since arrays are immutable-length in Julia).
"""
function _ws_array!(workspace::Nothing, backend, key::Symbol, ::Type{T}, dims) where {T}
    return KernelAbstractions.allocate(backend, T, dims)
end
function _ws_array!(workspace::MieBatchWorkspace, backend, key::Symbol, ::Type{T}, dims) where {T}
    _ws_owner_check!(workspace, backend)
    needed = prod(dims)
    cur = get(workspace.buffers, key, nothing)
    if cur === nothing || eltype(cur) !== T || length(cur) < needed
        cur = KernelAbstractions.allocate(backend, T, needed)
        workspace.buffers[key] = cur
    end
    return reshape(view(cur, 1:needed), dims)::AbstractArray{T}
end

# ============================================================================
# Grouping (host-side, per λ-call)
# ============================================================================

"""
    _group_by_nmax(offsets, radii_host, k_wavenum) -> (n_max_per_ensemble, groups)

Host-side, per-λ grouping: computes each ensemble's `n_max_global` at this
λ (`k_wavenum = 2π/λ`) from `radii_host` (always host-resident, see
[`MieNodeGeometry`](@ref)), then groups ensemble indices (in ORIGINAL order
within each group) by EXACT `n_max_global`. Returns
`(n_max_per_ensemble::Vector{Int}, groups::Vector{Pair{Int,Vector{Int}}})`
— `groups` sorted by `n_max` ascending purely for deterministic iteration
order (does not affect correctness: each group is processed independently).
"""
function _group_by_nmax(offsets::Vector{Int}, radii_host::Vector{<:AbstractFloat}, k_wavenum)
    E = length(offsets) - 1
    n_max_per_ensemble = Vector{Int}(undef, E)
    @inbounds for e in 1:E
        lo = offsets[e] + 1
        hi = offsets[e + 1]
        xmax = k_wavenum * maximum(view(radii_host, lo:hi))
        n_max_per_ensemble[e] = get_n_max(xmax)
    end
    groups_dict = Dict{Int,Vector{Int}}()
    @inbounds for e in 1:E
        push!(get!(() -> Int[], groups_dict, n_max_per_ensemble[e]), e)
    end
    groups = sort(collect(groups_dict); by = first)
    return n_max_per_ensemble, groups
end

# ============================================================================
# CPU batched path (Threads.@threads over ensembles -- NO grouping)
# ============================================================================

function _batched_cpu_full(geometry::MieNodeGeometry, λ, n_real::AbstractVector, n_imag::AbstractVector,
                            truncation::AbstractTruncationType, l_max::Union{Nothing,Integer})
    E = geometry.n_ensembles
    offsets = geometry.offsets
    r_all = geometry.radii_host   # IC = Float64 for a CPU geometry
    w_all = geometry.weights_host # pre-normalized, IC = Float64

    # FT_out is the CALLER's output type: promote_type of geometry's ORIGINAL
    # (pre-forcing) radii/weights eltype (`output_FT`, NOT `eltype(r_all)` --
    # `r_all` is geometry's INTERNAL Float64 compute mirror, forced to Float64
    # for every CPU geometry regardless of the caller's own radii/weights
    # eltype; using it here would silently return Float64 even for an
    # all-Float32 call -- the exact bug this comment guards against, see
    # MieNodeGeometry's `output_FT` field and this module's top-of-file
    # docstring) together with n_real/n_imag/λ, EXACTLY mirroring
    # `compute_aerosol_optical_properties_nodes`'s own `FT = float(promote_type(
    # eltype(radii), eltype(weights), typeof(n_real), typeof(n_imag),
    # typeof(wavelength)))` formula. Loop-invariant (not e-dependent) -- hoisted
    # out of the Threads.@threads loop both to avoid recomputing it E times and
    # so `results` can be a concretely typed `Vector{AerosolOptics{FT_out}}`
    # (an untyped `Vector{AerosolOptics}` would box every element, since
    # `AerosolOptics{FT}` is parametric).
    FT_out = float(promote_type(geometry.output_FT, eltype(n_real), eltype(n_imag), typeof(λ)))
    results = Vector{AerosolOptics{FT_out}}(undef, E)
    Threads.@threads for e in 1:E
        lo = offsets[e] + 1
        hi = offsets[e + 1]
        r_e = view(r_all, lo:hi)
        w_e = view(w_all, lo:hi)

        IC = Float64
        k_wavenum = IC(2π) / IC(λ)
        x_size_param = k_wavenum .* r_e
        n_max_e = get_n_max(maximum(x_size_param))
        n_mu_e = 2n_max_e - 1

        core = _nai2_bulk_optics(r_e, w_e, IC(n_real[e]), IC(n_imag[e]), IC(λ), n_max_e, n_mu_e, IC)

        raw = if FT_out <: AbstractFloat
            greek_coefs = GreekCoefs(convert.(FT_out, core.α), convert.(FT_out, core.β),
                                      convert.(FT_out, core.γ), convert.(FT_out, core.δ),
                                      convert.(FT_out, core.ϵ), convert.(FT_out, core.ζ))
            AerosolOptics(greek_coefs=greek_coefs, ω̃=FT_out(core.bulk_C_sca / core.bulk_C_ext),
                          k=FT_out(core.bulk_C_ext), fᵗ=FT_out(1))
        else
            greek_coefs = GreekCoefs(core.α, core.β, core.γ, core.δ, core.ϵ, core.ζ)
            AerosolOptics(greek_coefs=greek_coefs, ω̃=(core.bulk_C_sca / core.bulk_C_ext),
                          k=(core.bulk_C_ext), fᵗ=one(eltype(core.α)))
        end

        out = _apply_requested_truncation(raw, truncation)
        results[e] = l_max === nothing ? out : _slice_greek(out, l_max)
    end
    return results
end

function _batched_cpu_extinction(geometry::MieNodeGeometry, λ, n_real::AbstractVector, n_imag::AbstractVector)
    E = geometry.n_ensembles
    offsets = geometry.offsets
    r_all = geometry.radii_host
    w_all = geometry.weights_host

    # See _batched_cpu_full for why FT_out uses geometry.output_FT (the
    # caller's ORIGINAL radii/weights eltype) rather than eltype(r_all)
    # (geometry's internal Float64 compute mirror) -- using the latter would
    # silently return Float64 for an all-Float32 call. Hoisted (loop-invariant)
    # so `results` can be a concretely typed Vector{FT_out} rather than boxing
    # through a Vector{Any} + identity.() narrowing pass.
    FT_out = float(promote_type(geometry.output_FT, eltype(n_real), eltype(n_imag), typeof(λ)))
    results = Vector{FT_out}(undef, E)
    Threads.@threads for e in 1:E
        lo = offsets[e] + 1
        hi = offsets[e + 1]
        r_e = view(r_all, lo:hi)
        w_e = view(w_all, lo:hi)

        IC = Float64
        k_wavenum = IC(2π) / IC(λ)
        x_size_param = k_wavenum .* r_e
        n_max_e = get_n_max(maximum(x_size_param))

        m_ref = IC(n_real[e]) - IC(n_imag[e]) * im
        y_max = maximum(x_size_param) * abs(m_ref)
        nmx_max = round(Int, max(n_max_e, y_max) + 51)
        an = zeros(Complex{IC}, n_max_e)
        bn = zeros(Complex{IC}, n_max_e)
        Dₙ = zeros(Complex{IC}, nmx_max)
        n_ = IC.(2 .* collect(1:n_max_e) .+ 1)

        n_e = hi - lo + 1
        C_ext = zeros(IC, n_e)
        for i in 1:n_e
            n_max_i = get_n_max(x_size_param[i])
            _, _, C_ext_i = _mie_node_ab_and_C_ext!(an, bn, Dₙ, n_, x_size_param[i], n_max_i, m_ref, k_wavenum)
            @inbounds C_ext[i] = C_ext_i
        end
        bulk_C_ext = sum(w_e .* C_ext)

        results[e] = FT_out <: AbstractFloat ? FT_out(bulk_C_ext) : bulk_C_ext
    end
    return results
end

# ============================================================================
# GPU batched path (exact-nmax grouping)
# ============================================================================

"""
    _local_offsets(offsets, ens_ids) -> Vector{Int}

Group-local offsets (length `length(ens_ids)+1`, `local_offsets[1] == 0`),
one entry per group ensemble, derived from the FULL geometry's `offsets` for
just the ensembles in `ens_ids` (in the order given).
"""
function _local_offsets(offsets::Vector{Int}, ens_ids::Vector{Int})
    E_g = length(ens_ids)
    local_offsets = Vector{Int}(undef, E_g + 1)
    local_offsets[1] = 0
    @inbounds for j in 1:E_g
        e = ens_ids[j]
        local_offsets[j + 1] = local_offsets[j] + (offsets[e + 1] - offsets[e])
    end
    return local_offsets
end

"""
    _group_permutation(offsets, ens_ids, local_offsets) -> Vector{Int}

Host-built (cheap: pure `Int` bookkeeping, no float payload) permutation
index array of length `local_offsets[end]`: entry `local_offsets[j]+1:local_offsets[j+1]`
holds the FULL-geometry node indices `offsets[e]+1:offsets[e+1]` for
`e = ens_ids[j]` -- i.e. "which original node goes at this position in the
group's contiguous layout", preserving each ensemble's own intra-ensemble
node order untouched. Used as a GATHER index directly against the geometry's
DEVICE-resident `radii`/`weights` (`geometry.radii[perm]`/`geometry.weights[perm]`
-- standard `AbstractArray` fancy-indexing, a proper GPU gather kernel under
CUDA.jl/Metal.jl, not a host round-trip for the float payload), so groups
address noncontiguous ensembles via this permutation rather than by copying
radii/weights through the host.
"""
function _group_permutation(offsets::Vector{Int}, ens_ids::Vector{Int}, local_offsets::Vector{Int})
    n_group_nodes = local_offsets[end]
    perm = Vector{Int}(undef, n_group_nodes)
    @inbounds for (j, e) in enumerate(ens_ids)
        lo = offsets[e] + 1
        hi = offsets[e + 1]
        perm[(local_offsets[j] + 1):local_offsets[j + 1]] .= lo:hi
    end
    return perm
end

"""
    _check_batched_ft_lock(geometry, λ, n_real, n_imag) -> FT

`geometry` locks the batched GPU API to ONE explicit float type `FT =
eltype(geometry.radii_host)` (the type it was built with) -- a per-call
`wavelength`/`n_real`/`n_imag` of any OTHER element type throws
`ArgumentError` immediately rather than silently promoting or narrowing (a
silent promotion here would be the wrong kind of surprise for a GPU kernel
precision knob: it would change what device arrays/kernels actually run,
not just an output-rounding convenience). Build a new geometry via
[`prepare_mie_node_geometry`](@ref) for a different `FT`.
"""
function _check_batched_ft_lock(geometry::MieNodeGeometry, λ, n_real, n_imag)
    FT = eltype(geometry.radii_host)
    if typeof(λ) !== FT
        throw(ArgumentError(
            "geometry is locked to FT=$FT (from prepare_mie_node_geometry); wavelength " *
            "has type $(typeof(λ)) -- convert explicitly (e.g. $(FT)(λ)) rather than " *
            "relying on silent promotion."))
    end
    if eltype(n_real) !== FT || eltype(n_imag) !== FT
        throw(ArgumentError(
            "geometry is locked to FT=$FT (from prepare_mie_node_geometry); n_real/n_imag " *
            "have eltype $(eltype(n_real))/$(eltype(n_imag)) -- convert explicitly " *
            "(e.g. $(FT).(n_real)) rather than relying on silent promotion."))
    end
    return FT
end

function _batched_gpu_full(geometry::MieNodeGeometry, λ, n_real::AbstractVector, n_imag::AbstractVector,
                            truncation::AbstractTruncationType, l_max::Union{Nothing,Integer},
                            backend, precision_policy::MiePrecisionPolicy, workspace)
    FT = _check_batched_ft_lock(geometry, λ, n_real, n_imag)
    E = geometry.n_ensembles

    _check_policy_ft(precision_policy, FT)
    RA = _mie_reduction_type(precision_policy, FT)
    if _is_metal_backend(backend) && (FT === Float64 || RA === Float64)
        throw(ArgumentError(
            "Metal has no Float64 hardware support: $(nameof(typeof(precision_policy))) " *
            "with geometry element type $FT would need a Float64 device array. Use a " *
            "Float32 geometry with NativeFloat32() on MetalGPU()."))
    end

    k_wavenum = FT(2π / λ)
    offsets = geometry.offsets
    r_host = geometry.radii_host

    _, groups = _group_by_nmax(offsets, r_host, k_wavenum)

    m_re_all = FT.(n_real)
    m_im_all = FT.(.-n_imag)  # convention: n = n_real - i*n_imag

    # Concretely typed (AerosolOptics{FT}, not the abstract AerosolOptics) --
    # every ensemble's output is narrowed/passed through to the SAME `FT`
    # (geometry's locked float type), so this never boxes.
    results = Vector{AerosolOptics{FT}}(undef, E)

    for (n_max_g, ens_ids) in groups
        n_mu_g = 2n_max_g - 1
        E_g = length(ens_ids)
        local_offsets = _local_offsets(offsets, ens_ids)
        n_group_nodes = local_offsets[end]
        perm_host = _group_permutation(offsets, ens_ids, local_offsets)

        n2e_group_host = Vector{Int}(undef, n_group_nodes)
        @inbounds for j in 1:E_g
            n2e_group_host[(local_offsets[j] + 1):local_offsets[j + 1]] .= j
        end
        m_re_g = FT[m_re_all[e] for e in ens_ids]
        m_im_g = FT[m_im_all[e] for e in ens_ids]

        # --- Device-side gather: permutation index uploaded once (small,
        # pure Int bookkeeping); radii/weights GATHERED on-device via fancy
        # indexing against geometry's own device-resident arrays -- no host
        # round-trip for the float node payload. ---
        perm_dev = _ws_array!(workspace, backend, :bx_perm, Int, n_group_nodes)
        KernelAbstractions.copyto!(backend, perm_dev, perm_host)
        r_group_dev = _ws_array!(workspace, backend, :bx_rgrp, FT, n_group_nodes)
        w_group_dev = _ws_array!(workspace, backend, :bx_wgrp, FT, n_group_nodes)
        r_group_dev .= view(geometry.radii, perm_dev)
        w_group_dev .= view(geometry.weights, perm_dev)

        nmax_dev = _ws_array!(workspace, backend, :bx_nmax, Int, n_group_nodes)
        n2e_dev  = _ws_array!(workspace, backend, :bx_n2e, Int, n_group_nodes)
        mre_dev  = _ws_array!(workspace, backend, :bx_mre, FT, E_g)
        mim_dev  = _ws_array!(workspace, backend, :bx_mim, FT, E_g)
        an_dev   = _ws_array!(workspace, backend, :bx_an, Complex{FT}, (n_group_nodes, n_max_g))
        bn_dev   = _ws_array!(workspace, backend, :bx_bn, Complex{FT}, (n_group_nodes, n_max_g))

        x_dev = _ws_array!(workspace, backend, :bx_x, FT, n_group_nodes)
        x_dev .= k_wavenum .* r_group_dev   # on-device broadcast, no host round-trip

        nmaxk = nmax_per_node_kernel!(backend)
        nmaxk(nmax_dev, x_dev; ndrange=n_group_nodes)
        KernelAbstractions.synchronize(backend)

        KernelAbstractions.copyto!(backend, n2e_dev, n2e_group_host)
        KernelAbstractions.copyto!(backend, mre_dev, m_re_g)
        KernelAbstractions.copyto!(backend, mim_dev, m_im_g)
        fill!(an_dev, zero(Complex{FT}))
        fill!(bn_dev, zero(Complex{FT}))

        # nmx_max: UNUSED inside Kernel 1's body (see mie_coefficients_kernel_f64!/_ds!'s
        # docstrings -- the per-node recursion depth is recomputed per-thread from
        # that node's own x/m_re/m_im), so any sane upper bound is fine; this is a
        # cheap host-side estimate from the small per-ensemble RI arrays plus n_max_g
        # (no need for the (potentially large) per-node x values at all).
        nmx_max_g = round(Int, n_max_g * (1 + maximum(sqrt.(m_re_g .^ 2 .+ m_im_g .^ 2))) + 51)

        kernel1 = _mie_kernel1_batched(precision_policy, backend)
        kernel1(an_dev, bn_dev, x_dev, n2e_dev, mre_dev, mim_dev, nmax_dev, nmx_max_g; ndrange=n_group_nodes)
        KernelAbstractions.synchronize(backend)

        μ, w_μ = gausslegendre(n_mu_g)
        leg_π, leg_τ = compute_mie_π_τ(FT.(μ), n_max_g)

        f11_dev     = _ws_array!(workspace, backend, :bx_f11, FT, (n_mu_g, n_group_nodes))
        f33_dev     = _ws_array!(workspace, backend, :bx_f33, FT, (n_mu_g, n_group_nodes))
        f12_dev     = _ws_array!(workspace, backend, :bx_f12, FT, (n_mu_g, n_group_nodes))
        f34_dev     = _ws_array!(workspace, backend, :bx_f34, FT, (n_mu_g, n_group_nodes))
        C_sca_dev   = _ws_array!(workspace, backend, :bx_Csca, FT, n_group_nodes)
        C_ext_dev   = _ws_array!(workspace, backend, :bx_Cext, FT, n_group_nodes)
        leg_pi_dev  = _ws_array!(workspace, backend, :bx_legpi, FT, (n_mu_g, n_max_g))
        leg_tau_dev = _ws_array!(workspace, backend, :bx_legtau, FT, (n_mu_g, n_max_g))
        KernelAbstractions.copyto!(backend, leg_pi_dev, FT.(leg_π))
        KernelAbstractions.copyto!(backend, leg_tau_dev, FT.(leg_τ))

        kernel23 = amplitude_phase_kernel!(backend)
        kernel23(f11_dev, f33_dev, f12_dev, f34_dev, an_dev, bn_dev, leg_pi_dev, leg_tau_dev,
                 x_dev, nmax_dev; ndrange=(n_mu_g, n_group_nodes))
        KernelAbstractions.synchronize(backend)

        kernel4a = cross_sections_kernel!(backend)
        kernel4a(C_sca_dev, C_ext_dev, an_dev, bn_dev, k_wavenum, nmax_dev; ndrange=n_group_nodes)
        KernelAbstractions.synchronize(backend)

        w_dev  = _ws_array!(workspace, backend, :bx_w, RA, n_group_nodes)
        wr_dev = _ws_array!(workspace, backend, :bx_wr, RA, n_group_nodes)
        w_dev  .= RA.(w_group_dev)
        # NOT `RA.(r_group_dev) .^ 2` (pre-widening before squaring): the
        # single-ensemble reference computes `4π .* r.^2 .* w` with `r`/`w`
        # STILL in native FT, so `r[i]^2` happens in FT precision FIRST and
        # only gets promoted to Float64 by the (untyped, hence Float64) `4π`
        # literal in the SAME broadcast fusion -- squaring pre-widened-to-RA
        # values instead is a materially different (if more "accurate")
        # computation, not bitwise equal to the single-ensemble path.
        wr_dev .= 4π .* r_group_dev .^ 2 .* w_group_dev
        loff_dev = _ws_array!(workspace, backend, :bx_loff, Int, E_g + 1)
        KernelAbstractions.copyto!(backend, loff_dev, local_offsets)

        bulk_Csca_dev = _ws_array!(workspace, backend, :bx_bCsca, RA, E_g)
        bulk_Cext_dev = _ws_array!(workspace, backend, :bx_bCext, RA, E_g)
        wsum = segmented_weighted_sum_kernel!(backend)
        wsum(bulk_Csca_dev, C_sca_dev, w_dev, loff_dev; ndrange=E_g)
        wsum(bulk_Cext_dev, C_ext_dev, w_dev, loff_dev; ndrange=E_g)
        KernelAbstractions.synchronize(backend)

        bulk_f11_dev = _ws_array!(workspace, backend, :bx_bf11, RA, (n_mu_g, E_g))
        bulk_f33_dev = _ws_array!(workspace, backend, :bx_bf33, RA, (n_mu_g, E_g))
        bulk_f12_dev = _ws_array!(workspace, backend, :bx_bf12, RA, (n_mu_g, E_g))
        bulk_f34_dev = _ws_array!(workspace, backend, :bx_bf34, RA, (n_mu_g, E_g))
        sred = segmented_size_reduction_kernel!(backend)
        sred(bulk_f11_dev, f11_dev, wr_dev, loff_dev; ndrange=(n_mu_g, E_g))
        sred(bulk_f33_dev, f33_dev, wr_dev, loff_dev; ndrange=(n_mu_g, E_g))
        sred(bulk_f12_dev, f12_dev, wr_dev, loff_dev; ndrange=(n_mu_g, E_g))
        sred(bulk_f34_dev, f34_dev, wr_dev, loff_dev; ndrange=(n_mu_g, E_g))
        KernelAbstractions.synchronize(backend)

        bulk_C_sca_host = Array(bulk_Csca_dev)
        bulk_C_ext_host = Array(bulk_Cext_dev)
        inv_bulk_C_sca_dev = _ws_array!(workspace, backend, :bx_invCsca, RA, E_g)
        KernelAbstractions.copyto!(backend, inv_bulk_C_sca_dev, one(RA) ./ bulk_C_sca_host)
        inv_row = reshape(inv_bulk_C_sca_dev, 1, E_g)
        bulk_f11_dev .*= inv_row
        bulk_f33_dev .*= inv_row
        bulk_f12_dev .*= inv_row
        bulk_f34_dev .*= inv_row

        P, P², R², T² = compute_legendre_poly(RA.(μ), n_mu_g)
        P_dev   = _ws_array!(workspace, backend, :bx_P, RA, (n_mu_g, n_mu_g))
        P2_dev  = _ws_array!(workspace, backend, :bx_P2, RA, (n_mu_g, n_mu_g))
        R2_dev  = _ws_array!(workspace, backend, :bx_R2, RA, (n_mu_g, n_mu_g))
        T2_dev  = _ws_array!(workspace, backend, :bx_T2, RA, (n_mu_g, n_mu_g))
        wmu_dev = _ws_array!(workspace, backend, :bx_wmu, RA, n_mu_g)
        KernelAbstractions.copyto!(backend, P_dev,  RA.(P))
        KernelAbstractions.copyto!(backend, P2_dev, RA.(P²))
        KernelAbstractions.copyto!(backend, R2_dev, RA.(R²))
        KernelAbstractions.copyto!(backend, T2_dev, RA.(T²))
        KernelAbstractions.copyto!(backend, wmu_dev, RA.(w_μ))

        α_dev = _ws_array!(workspace, backend, :bx_alpha, RA, (n_mu_g, E_g))
        β_dev = _ws_array!(workspace, backend, :bx_beta,  RA, (n_mu_g, E_g))
        γ_dev = _ws_array!(workspace, backend, :bx_gamma, RA, (n_mu_g, E_g))
        δ_dev = _ws_array!(workspace, backend, :bx_delta, RA, (n_mu_g, E_g))
        ϵ_dev = _ws_array!(workspace, backend, :bx_eps,   RA, (n_mu_g, E_g))
        ζ_dev = _ws_array!(workspace, backend, :bx_zeta,  RA, (n_mu_g, E_g))

        greek_kernel = greek_coefficients_kernel_batched!(backend)
        greek_kernel(α_dev, β_dev, γ_dev, δ_dev, ϵ_dev, ζ_dev,
                     bulk_f11_dev, bulk_f33_dev, bulk_f12_dev, bulk_f34_dev,
                     P_dev, P2_dev, R2_dev, T2_dev, wmu_dev, n_mu_g; ndrange=(n_mu_g, E_g))
        KernelAbstractions.synchronize(backend)

        α_host = Array(α_dev); β_host = Array(β_dev); γ_host = Array(γ_dev)
        δ_host = Array(δ_dev); ϵ_host = Array(ϵ_dev); ζ_host = Array(ζ_dev)

        for j in 1:E_g
            e = ens_ids[j]
            raw = if RA !== FT
                greek_coefs = GreekCoefs(convert.(FT, α_host[:, j]), convert.(FT, β_host[:, j]),
                                          convert.(FT, γ_host[:, j]), convert.(FT, δ_host[:, j]),
                                          convert.(FT, ϵ_host[:, j]), convert.(FT, ζ_host[:, j]))
                AerosolOptics(greek_coefs=greek_coefs,
                              ω̃=FT(bulk_C_sca_host[j] / bulk_C_ext_host[j]),
                              k=FT(bulk_C_ext_host[j]), fᵗ=FT(1))
            else
                greek_coefs = GreekCoefs(α_host[:, j], β_host[:, j], γ_host[:, j],
                                          δ_host[:, j], ϵ_host[:, j], ζ_host[:, j])
                AerosolOptics(greek_coefs=greek_coefs,
                              ω̃=(bulk_C_sca_host[j] / bulk_C_ext_host[j]),
                              k=bulk_C_ext_host[j], fᵗ=one(FT))
            end
            out = _apply_requested_truncation(raw, truncation)
            results[e] = l_max === nothing ? out : _slice_greek(out, l_max)
        end
    end

    return results
end

function _batched_gpu_extinction(geometry::MieNodeGeometry, λ, n_real::AbstractVector, n_imag::AbstractVector,
                                  backend, precision_policy::MiePrecisionPolicy, workspace)
    FT = _check_batched_ft_lock(geometry, λ, n_real, n_imag)
    E = geometry.n_ensembles

    _check_policy_ft(precision_policy, FT)
    RA = _mie_reduction_type(precision_policy, FT)
    if _is_metal_backend(backend) && (FT === Float64 || RA === Float64)
        throw(ArgumentError(
            "Metal has no Float64 hardware support: $(nameof(typeof(precision_policy))) " *
            "with geometry element type $FT would need a Float64 device array. Use a " *
            "Float32 geometry with NativeFloat32() on MetalGPU()."))
    end

    k_wavenum = FT(2π / λ)
    offsets = geometry.offsets
    r_host = geometry.radii_host

    _, groups = _group_by_nmax(offsets, r_host, k_wavenum)

    m_re_all = FT.(n_real)
    m_im_all = FT.(.-n_imag)

    results = Vector{FT}(undef, E)

    for (n_max_g, ens_ids) in groups
        E_g = length(ens_ids)
        local_offsets = _local_offsets(offsets, ens_ids)
        n_group_nodes = local_offsets[end]
        perm_host = _group_permutation(offsets, ens_ids, local_offsets)

        n2e_group_host = Vector{Int}(undef, n_group_nodes)
        @inbounds for j in 1:E_g
            n2e_group_host[(local_offsets[j] + 1):local_offsets[j + 1]] .= j
        end
        m_re_g = FT[m_re_all[e] for e in ens_ids]
        m_im_g = FT[m_im_all[e] for e in ens_ids]

        # Shared key namespace with _batched_gpu_full (review finding 5): every
        # buffer below has the IDENTICAL shape/dtype/role in both functions
        # (the Kernel-1 preamble is common to full and extinction-only), so a
        # workspace reused across both APIs (the documented GCHPIO hybrid
        # exact-τ pattern) holds ONE copy of each -- not two -- including
        # an/bn, the largest allocation. Buffers unique to the full path
        # (amplitude/phase/Greek scratch) keep their own `:bx_*` keys below.
        perm_dev = _ws_array!(workspace, backend, :bx_perm, Int, n_group_nodes)
        KernelAbstractions.copyto!(backend, perm_dev, perm_host)
        r_group_dev = _ws_array!(workspace, backend, :bx_rgrp, FT, n_group_nodes)
        w_group_dev = _ws_array!(workspace, backend, :bx_wgrp, FT, n_group_nodes)
        r_group_dev .= view(geometry.radii, perm_dev)
        w_group_dev .= view(geometry.weights, perm_dev)

        nmax_dev = _ws_array!(workspace, backend, :bx_nmax, Int, n_group_nodes)
        n2e_dev  = _ws_array!(workspace, backend, :bx_n2e, Int, n_group_nodes)
        mre_dev  = _ws_array!(workspace, backend, :bx_mre, FT, E_g)
        mim_dev  = _ws_array!(workspace, backend, :bx_mim, FT, E_g)
        an_dev   = _ws_array!(workspace, backend, :bx_an, Complex{FT}, (n_group_nodes, n_max_g))
        bn_dev   = _ws_array!(workspace, backend, :bx_bn, Complex{FT}, (n_group_nodes, n_max_g))

        x_dev = _ws_array!(workspace, backend, :bx_x, FT, n_group_nodes)
        x_dev .= k_wavenum .* r_group_dev

        nmaxk = nmax_per_node_kernel!(backend)
        nmaxk(nmax_dev, x_dev; ndrange=n_group_nodes)
        KernelAbstractions.synchronize(backend)

        KernelAbstractions.copyto!(backend, n2e_dev, n2e_group_host)
        KernelAbstractions.copyto!(backend, mre_dev, m_re_g)
        KernelAbstractions.copyto!(backend, mim_dev, m_im_g)
        fill!(an_dev, zero(Complex{FT}))
        fill!(bn_dev, zero(Complex{FT}))

        nmx_max_g = round(Int, n_max_g * (1 + maximum(sqrt.(m_re_g .^ 2 .+ m_im_g .^ 2))) + 51)

        kernel1 = _mie_kernel1_batched(precision_policy, backend)
        kernel1(an_dev, bn_dev, x_dev, n2e_dev, mre_dev, mim_dev, nmax_dev, nmx_max_g; ndrange=n_group_nodes)
        KernelAbstractions.synchronize(backend)

        C_sca_dev = _ws_array!(workspace, backend, :bx_Csca, FT, n_group_nodes)
        C_ext_dev = _ws_array!(workspace, backend, :bx_Cext, FT, n_group_nodes)
        kernel4a = cross_sections_kernel!(backend)
        kernel4a(C_sca_dev, C_ext_dev, an_dev, bn_dev, k_wavenum, nmax_dev; ndrange=n_group_nodes)
        KernelAbstractions.synchronize(backend)

        w_dev = _ws_array!(workspace, backend, :bx_w, RA, n_group_nodes)
        w_dev .= RA.(w_group_dev)
        loff_dev = _ws_array!(workspace, backend, :bx_loff, Int, E_g + 1)
        KernelAbstractions.copyto!(backend, loff_dev, local_offsets)

        bulk_Cext_dev = _ws_array!(workspace, backend, :bx_bCext, RA, E_g)
        wsum = segmented_weighted_sum_kernel!(backend)
        wsum(bulk_Cext_dev, C_ext_dev, w_dev, loff_dev; ndrange=E_g)
        KernelAbstractions.synchronize(backend)

        bulk_C_ext_host = Array(bulk_Cext_dev)
        for j in 1:E_g
            e = ens_ids[j]
            results[e] = RA !== FT ? FT(bulk_C_ext_host[j]) : bulk_C_ext_host[j]
        end
    end

    return results
end

# ============================================================================
# Public entry points
# ============================================================================

const _WARNED_NO_GPU_MIE_NODES_BATCHED = Ref(false)
const _WARNED_NO_GPU_MIE_EXTINCTION_NODES_BATCHED = Ref(false)

"""
    _validate_batched_node_inputs(geometry::MieNodeGeometry, n_real, n_imag)

Shared per-call input validation for all four batched public entry points
(full/extinction × top-level/`_gpu`): `n_real`/`n_imag` must have length
`geometry.n_ensembles`, and every `n_imag` must be `≥ 0` (convention
`n = n_real - i·n_imag`) -- mirrors this module's sibling
`_validate_node_inputs` (`compute_NAI2_nodes.jl`) convention of a single
`@assert`-based helper rather than repeating the same two checks at every
call site.
"""
@inline function _validate_batched_node_inputs(geometry::MieNodeGeometry, n_real, n_imag)
    E = geometry.n_ensembles
    @assert length(n_real) == E && length(n_imag) == E "n_real/n_imag must have length == geometry.n_ensembles ($E)"
    @assert all(x -> x ≥ 0, n_imag) "Imaginary part of the refractive index must be ≥ 0 for every ensemble"
    return nothing
end

@doc raw"""
    compute_aerosol_optical_properties_nodes_batched(geometry::MieNodeGeometry,
        wavelength, n_real::AbstractVector, n_imag::AbstractVector,
        polarization, truncation;
        l_max = nothing,
        precision_policy = nothing,
        workspace = nothing) -> Vector{AerosolOptics}

Batched sibling of [`compute_aerosol_optical_properties_nodes`](@ref): given a
[`MieNodeGeometry`](@ref) (`E` ensembles, prepared once via
[`prepare_mie_node_geometry`](@ref)) and ONE wavelength, returns a
`Vector{AerosolOptics}` of length `E` — element `e` is what
`compute_aerosol_optical_properties_nodes(radii_e, weights_e, n_real[e],
n_imag[e], wavelength, polarization, truncation; l_max, precision_policy)`
would return for ensemble `e`'s own (raw) nodes, in ORIGINAL ensemble order.
**Output type**: `AerosolOptics{FT_out}` where `FT_out` is derived from
`geometry`'s ORIGINAL `radii`/`weights` eltype (`MieNodeGeometry.output_FT`,
captured before any internal-compute-precision forcing) promoted with
`n_real`/`n_imag`/`wavelength`'s eltypes — an all-`Float32` CPU call returns
`Float32`, matching `compute_aerosol_optical_properties_nodes` exactly (this
is a load-bearing contract: GCHPIO's consumer code narrows to its own
working precision and must be able to rely on this without a manual
post-hoc narrowing pass).

# Bit-compatibility contract (precise, by backend)
- **`Architectures.CPU()` and KA-CPU** (portable KernelAbstractions backends
  that aren't real CUDA/Metal hardware): bit-for-bit identical (`==`) to the
  single-ensemble path, by construction (same arithmetic order, same
  per-ensemble normalize-once semantics — see this file's top-of-file module
  docstring for the full argument).
- **Real CUDA hardware**: near-machine-precision, NOT bit-for-bit (different
  `ndrange` launch shapes between the batched and single-ensemble kernel
  invocations can produce different, but equally valid, PTX/SASS for
  logically identical code — see `docs/src/pages/gpu.md`'s "Batched
  caller-node Mie seam" section for the measured magnitudes).
- **Exact even on real CUDA**: `compute_aerosol_extinction_nodes_batched`'s
  result for ensemble `e` equals THIS function's `.k` for that same ensemble
  — the invariant GCHPIO's hybrid exact-τ mechanism depends on (see
  [`compute_aerosol_extinction_nodes_batched`](@ref)).

See `test/test_mie_nodes_batched.jl` for the acceptance-gate tests (bitwise
batched-vs-single on CPU/KA-CPU, batching-split invariance, workspace-reuse
invariance, original-order reassembly, and the dedicated strict real-CUDA
`.k`-equality regression test).

# GPU: exact-`n_max_global` grouping
For a GPU `geometry.architecture`, ensembles are grouped by EXACT
`n_max_global` at THIS wavelength (λ-dependent — group membership may differ
between wavelengths; recomputed fresh every call), and one kernel-1 launch
covers each group's concatenated nodes (per-ensemble refractive index via a
`node_to_ensemble` index — see `gpu_mie_kernels.jl`'s `_batched` sibling
kernels).

# CPU: no grouping
For `Architectures.CPU()`, this is a `Threads.@threads` loop over ensembles
(grouping is a GPU-only concern — see the module docstring for why a single
sequential "CPU batch" would regress the measured 16-thread performance this
proposal is built around).

# Arguments
- `geometry`: prepared once via [`prepare_mie_node_geometry`](@ref); reused
  across as many calls (different `wavelength`s, `n_real`/`n_imag`, etc.) as
  needed without re-normalizing or re-uploading anything.
- `wavelength`: ONE wavelength for this call (μm, matching `geometry`'s radii
  units).
- `n_real::AbstractVector`, `n_imag::AbstractVector`: length-`E` per-ensemble
  refractive index (convention `n = n_real - i·n_imag`, `n_imag ≥ 0` for
  every ensemble), one complex value per ensemble (shared by all of that
  ensemble's own nodes — the same "one RI per ensemble" convention as the
  single-ensemble API).
- `polarization`, `truncation`: same meaning as
  [`compute_aerosol_optical_properties_nodes`](@ref); `truncation` (δBGE then
  `l_max` cap ordering) is applied PER ENSEMBLE on the host after the raw
  batched compute.

# Keyword arguments
- `l_max`: same meaning as the single-ensemble API, applied per ensemble.
- `precision_policy`: GPU-only, see [`MiePrecisionPolicy`](@ref); `nothing`
  auto-selects via `Architectures.default_mie_precision_policy(architecture,
  eltype(geometry.radii_host))`. Ignored on the CPU path.
- `workspace`: `nothing` (default, allocate-and-drop) or a
  [`MieBatchWorkspace`](@ref) to reuse device scratch across calls.
"""
function compute_aerosol_optical_properties_nodes_batched(
    geometry::MieNodeGeometry, wavelength::Real,
    n_real::AbstractVector, n_imag::AbstractVector,
    polarization::AbstractPolarizationType,
    truncation::AbstractTruncationType;
    l_max::Union{Nothing,Integer} = nothing,
    precision_policy::Union{Nothing,MiePrecisionPolicy} = nothing,
    workspace = nothing,
)
    _validate_batched_node_inputs(geometry, n_real, n_imag)

    architecture = geometry.architecture
    if architecture isa Architectures.CPU
        return _batched_cpu_full(geometry, wavelength, n_real, n_imag, truncation, l_max)
    elseif Architectures.has_gpu_mie(architecture)
        FT = eltype(geometry.radii_host)
        # geometry.backend -- the ACTUAL backend `geometry`'s device arrays
        # were built on -- NOT a fresh Architectures.ka_backend(architecture)
        # re-derivation: `prepare_mie_node_geometry`'s `backend` keyword can
        # override the natural backend (see its docstring), and re-deriving
        # here would silently ignore that override, launching kernels against
        # a DIFFERENT backend than the one geometry's arrays actually live on.
        backend = geometry.backend
        policy = precision_policy === nothing ?
                 Architectures.default_mie_precision_policy(architecture, FT) : precision_policy
        return _batched_gpu_full(geometry, wavelength, n_real, n_imag, truncation, l_max, backend, policy, workspace)
    else
        if !_WARNED_NO_GPU_MIE_NODES_BATCHED[]
            @warn "no GPU Mie pipeline for $(architecture); computing batched caller-node Mie on CPU"
            _WARNED_NO_GPU_MIE_NODES_BATCHED[] = true
        end
        return _batched_cpu_full(geometry, wavelength, n_real, n_imag, truncation, l_max)
    end
end

@doc raw"""
    compute_aerosol_optical_properties_nodes_batched_gpu(geometry::MieNodeGeometry,
        wavelength, n_real::AbstractVector, n_imag::AbstractVector,
        polarization, truncation, backend;
        l_max = nothing,
        precision_policy = NativeFloat64(),
        workspace = nothing) -> Vector{AerosolOptics}

GPU/KernelAbstractions implementation backing
[`compute_aerosol_optical_properties_nodes_batched`](@ref) — takes an
EXPLICIT KA `backend` (e.g. `KernelAbstractions.CPU()` or
`CUDA.CUDABackend()`), mirroring
[`compute_aerosol_optical_properties_nodes_gpu`](@ref)'s directly-callable
convention so it can be exercised on a CPU-backed KA backend in tests without
a real GPU device, INDEPENDENTLY of `geometry.architecture` (which only
determines what the top-level [`compute_aerosol_optical_properties_nodes_batched`](@ref)
dispatcher picks automatically). This is the entry point the acceptance-gate
tests use to compare, e.g., the KA-CPU-backed batched path against the
KA-CPU-backed single-ensemble path bit-for-bit.
"""
function compute_aerosol_optical_properties_nodes_batched_gpu(
    geometry::MieNodeGeometry, wavelength::Real,
    n_real::AbstractVector, n_imag::AbstractVector,
    polarization::AbstractPolarizationType,
    truncation::AbstractTruncationType, backend;
    l_max::Union{Nothing,Integer} = nothing,
    precision_policy::MiePrecisionPolicy = NativeFloat64(),
    workspace = nothing,
)
    _validate_batched_node_inputs(geometry, n_real, n_imag)
    return _batched_gpu_full(geometry, wavelength, n_real, n_imag, truncation, l_max, backend, precision_policy, workspace)
end

@doc raw"""
    compute_aerosol_extinction_nodes_batched(geometry::MieNodeGeometry,
        wavelength, n_real::AbstractVector, n_imag::AbstractVector;
        precision_policy = nothing, workspace = nothing) -> Vector

Batched, extinction-only sibling of [`compute_aerosol_extinction_nodes`](@ref)
— same `geometry`/grouping/workspace story as
[`compute_aerosol_optical_properties_nodes_batched`](@ref), but returns just
`Vector{k}` (length `E`, one number-mean bulk extinction cross-section per
ensemble, ORIGINAL order, output type from `geometry.output_FT` promoted
with `n_real`/`n_imag`/`wavelength` -- same load-bearing Float32-output
contract as the full API, see its docstring), computing none of the
amplitude matrices, phase-matrix elements, or Greek coefficients — see
[`compute_aerosol_extinction_nodes`](@ref)'s own docstring for the
motivation (GCHPIO's hybrid exact-τ leg is the intended consumer of this
batched variant).

# Exact-even-on-CUDA invariant
`compute_aerosol_extinction_nodes_batched`'s result for ensemble `e` equals
[`compute_aerosol_optical_properties_nodes_batched`](@ref)'s `.k` for that
SAME ensemble EXACTLY (`==`), on EVERY backend including real CUDA, for all
three precision policies — see that function's "Bit-compatibility contract"
section for why (shared Kernel-1 + cross-section + segmented-weighted-sum
code path). GCHPIO's hybrid exact-τ mechanism rests on this; see the
dedicated strict real-CUDA regression test in `test/test_mie_nodes_batched.jl`.
"""
function compute_aerosol_extinction_nodes_batched(
    geometry::MieNodeGeometry, wavelength::Real,
    n_real::AbstractVector, n_imag::AbstractVector;
    precision_policy::Union{Nothing,MiePrecisionPolicy} = nothing,
    workspace = nothing,
)
    _validate_batched_node_inputs(geometry, n_real, n_imag)

    architecture = geometry.architecture
    if architecture isa Architectures.CPU
        return _batched_cpu_extinction(geometry, wavelength, n_real, n_imag)
    elseif Architectures.has_gpu_mie(architecture)
        FT = eltype(geometry.radii_host)
        # geometry.backend -- the ACTUAL backend `geometry`'s device arrays
        # were built on -- NOT a fresh Architectures.ka_backend(architecture)
        # re-derivation: `prepare_mie_node_geometry`'s `backend` keyword can
        # override the natural backend (see its docstring), and re-deriving
        # here would silently ignore that override, launching kernels against
        # a DIFFERENT backend than the one geometry's arrays actually live on.
        backend = geometry.backend
        policy = precision_policy === nothing ?
                 Architectures.default_mie_precision_policy(architecture, FT) : precision_policy
        return _batched_gpu_extinction(geometry, wavelength, n_real, n_imag, backend, policy, workspace)
    else
        if !_WARNED_NO_GPU_MIE_EXTINCTION_NODES_BATCHED[]
            @warn "no GPU Mie pipeline for $(architecture); computing batched caller-node Mie extinction on CPU"
            _WARNED_NO_GPU_MIE_EXTINCTION_NODES_BATCHED[] = true
        end
        return _batched_cpu_extinction(geometry, wavelength, n_real, n_imag)
    end
end

@doc raw"""
    compute_aerosol_extinction_nodes_batched_gpu(geometry::MieNodeGeometry,
        wavelength, n_real::AbstractVector, n_imag::AbstractVector, backend;
        precision_policy = NativeFloat64(), workspace = nothing) -> Vector

GPU/KernelAbstractions implementation backing
[`compute_aerosol_extinction_nodes_batched`](@ref) — explicit-`backend`
directly-callable counterpart, mirroring
[`compute_aerosol_extinction_nodes_gpu`](@ref)'s convention, independent of
`geometry.architecture` (see
[`compute_aerosol_optical_properties_nodes_batched_gpu`](@ref) for why this
entry point exists — it is what the acceptance-gate tests use for
same-backend batched-vs-single comparisons).
"""
function compute_aerosol_extinction_nodes_batched_gpu(
    geometry::MieNodeGeometry, wavelength::Real,
    n_real::AbstractVector, n_imag::AbstractVector, backend;
    precision_policy::MiePrecisionPolicy = NativeFloat64(),
    workspace = nothing,
)
    _validate_batched_node_inputs(geometry, n_real, n_imag)
    return _batched_gpu_extinction(geometry, wavelength, n_real, n_imag, backend, precision_policy, workspace)
end
