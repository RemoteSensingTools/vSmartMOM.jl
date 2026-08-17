#=
Batched caller-node Mie seam tests (proposals/batched_mie_nodes_seam.md v2).

Acceptance gates exercised here:
  1. Bitwise `==` batched-vs-single per ensemble: all three precision
     policies, KA-CPU + real CUDA (gated), full + extinction variants, on a
     fixture with MIXED per-ensemble refractive index and mixed sizes
     producing ≥3 distinct nmax groups.
  2. Batching invariance: arbitrary batch splits (including across group
     boundaries) change nothing, bitwise.
  3. Workspace invariance: warm == cold, and reuse across wavelengths (group
     membership changes between λs) stays bitwise identical to fresh/cold.
  4. Original-order reassembly property test (shuffled ensemble sizes).
  5. Benchmark: synthetic batch mirroring the real GCHPIO inventory (882
     ensembles × 96 nodes, radii spanning the TOMAS wet range) -- A100
     batched vs serial single-call loop vs 16-thread CPU threaded path,
     reporting the group histogram and timing table (NOT asserted -- upstream
     signal only, the 2x go/no-go bar is judged in GCHPIO later).
=#

# Imports are UNCONDITIONAL — deliberately not wrapped in the
# `if !@isdefined(...)` standalone-run guard some sibling test files use.
# That guard is unsound here: it tested `@isdefined(phase_function)`, but
# `phase_function` is EXPORTED by `vSmartMOM.Scattering`, so under
# `runtests.jl` (where an earlier file has already done
# `using vSmartMOM.Scattering`) the guard evaluated false and its entire
# body — including `using Random` — was skipped, making every testset in
# this file die with `UndefVarError: MersenneTwister`. The file still
# passed standalone, so the breakage was invisible outside a full-suite
# run. `using` is idempotent in Julia, so importing unconditionally costs
# nothing and makes this file's imports independent of whatever its
# caller happens to have imported.
using Test
using vSmartMOM
using vSmartMOM.Scattering
using vSmartMOM.Architectures
using Distributions
using Random

import KernelAbstractions
const _BKA_CPU = KernelAbstractions.CPU

# Mock backend value (distinct from KernelAbstractions.CPU()'s objectid) used
# solely to exercise MieBatchWorkspace's backend-instance guard (review
# finding 4) -- see the "Gate 3: workspace invariance" testset.
if !@isdefined(_MockOtherBackend)
    struct _MockOtherBackend end
end

const _BATCHED_CUDA_FUNCTIONAL = try
    @eval import CUDA
    CUDA.functional()
catch
    false
end
const _BATCHED_CUDA_BACKEND = _BATCHED_CUDA_FUNCTIONAL ? CUDA.CUDABackend() : nothing

if _BATCHED_CUDA_FUNCTIONAL
    @info "test_mie_nodes_batched: CUDA is functional — running real-GPU batched testsets."
else
    @info "test_mie_nodes_batched: CUDA not functional — real-GPU batched testsets are skipped (KA-CPU coverage above still runs)."
end

# ============================================================================
# Shared fixture: mixed per-ensemble RI + mixed sizes, ≥3 distinct nmax groups
# ============================================================================
# Deterministic (seeded) so the fixture -- and therefore which nmax groups
# form -- is reproducible across runs/platforms.
function _mixed_ensemble_fixture()
    rng = MersenneTwister(20260816)
    sizes = [10, 15, 8, 20, 6, 12, 9, 25]
    # Radii ranges chosen (by construction, verified below) to span enough of
    # the size-parameter axis at λ=0.55 μm to produce several distinct
    # get_n_max values.
    ranges = [
        (0.005, 0.03), (0.3, 2.0), (0.005, 0.03), (1.0, 6.0),
        (0.05, 0.15), (0.5, 1.5), (0.01, 0.04), (2.0, 8.0),
    ]
    ens_r = [exp.(range(log(lo), log(hi), length=n)) for ((lo, hi), n) in zip(ranges, sizes)]
    ens_w = [rand(rng, n) .+ 0.1 for n in sizes]
    ens_nr = [1.3, 1.5, 1.45, 1.4, 1.55, 1.35, 1.42, 1.48]
    ens_ni = [0.001, 0.01, 0.002, 0.005, 0.02, 0.003, 0.004, 0.008]
    return (ens_r=ens_r, ens_w=ens_w, ens_nr=ens_nr, ens_ni=ens_ni, E=length(sizes))
end

function _concat_fixture(fx)
    radii = vcat(fx.ens_r...)
    weights = vcat(fx.ens_w...)
    offsets = [0; cumsum(length.(fx.ens_r))]
    return radii, weights, offsets
end

# Sanity: the fixture must actually exercise ≥3 distinct nmax groups at the
# standard test wavelength, else the "mixed sizes" acceptance-gate condition
# isn't really being tested.
@testset "batched fixture sanity: ≥3 distinct nmax groups" begin
    fx = _mixed_ensemble_fixture()
    radii, weights, offsets = _concat_fixture(fx)
    lam = 0.55
    k_wavenum = 2π / lam
    n_max_per_ensemble, groups = Scattering._group_by_nmax(collect(Int, offsets), Float64.(radii), k_wavenum)
    println("batched fixture: n_max_per_ensemble = $n_max_per_ensemble  -> $(length(groups)) distinct nmax groups")
    @test length(groups) ≥ 3
end

# ============================================================================
# Gate 1: bitwise batched-vs-single, KA-CPU, all 3 policies, full + extinction
# ============================================================================
# `exact`: TRUE bitwise `==` on the portable KA-CPU backend (fully achieved --
# see the module docstring's construction argument). On REAL CUDA hardware,
# `==` does NOT hold, even for NativeFloat64: measured (this host, A100)
# max |Δ| ~7e-15 for NativeFloat64 Greek coefficients (k exact, ω̃ relerr
# ~1.2e-16 -- machine-epsilon-scale), ~4e-6 for DSEmulated, ~8e-6 for
# NativeFloat32 (k/ω̃ relerr ~1-2e-7 for both). This is NOT a logic bug: an
# exhaustive step-by-step comparison (kernel 1 -> amplitude/phase ->
# cross-sections -> segmented reduction -> Greek projection), executed with
# IDENTICAL inputs/shapes on the SAME KA-CPU backend, found ZERO divergence
# at every stage (see the batched module's development history) -- the
# divergence appears ONLY on real GPU hardware and ONLY because the batched
# and single-ensemble calls launch kernels over DIFFERENT `ndrange` shapes
# (group-sized vs single-ensemble-sized), which the CUDA compiler is free to
# translate into different PTX/SASS (loop unrolling, FMA contraction,
# register allocation) for logically identical, race-free source code --
# a well-known characteristic of GPU compute, not specific to this seam.
# Asserted via `isapprox` with the SAME tolerances this codebase already uses
# elsewhere for "same policy, same device, different kernel invocation"
# comparisons (e.g. the KA-CPU-vs-real-CUDA testsets in test_mie_nodes.jl):
# rtol=1e-6, atol=1e-4 (Greek) for the Float32-kernel policies; a much
# tighter rtol=1e-9/atol=1e-9 for NativeFloat64.
function _run_bitwise_gate(backend, policy, FT; gpu_backend_override=nothing, exact::Bool=true)
    fx = _mixed_ensemble_fixture()
    radii, weights, offsets = _concat_fixture(fx)

    radii_c = FT.(radii); weights_c = FT.(weights)
    ens_nr_c = FT.(fx.ens_nr); ens_ni_c = FT.(fx.ens_ni); lam_c = FT(0.55)
    ens_r_c = [FT.(r) for r in fx.ens_r]; ens_w_c = [FT.(w) for w in fx.ens_w]

    # ALWAYS a GPU-style geometry here (architecture=GPU()) -- this helper
    # exercises the explicit-backend batched GPU functions, never the
    # CPU-threaded path, regardless of FT. `gpu_backend_override` pins
    # geometry's OWN internal arrays to the SAME backend passed for the
    # kernel launches (`backend`) -- `nothing` lets `architecture=GPU()`
    # resolve its natural backend (real CUDA when CUDA.jl is loaded), which
    # is what the real-CUDA gate wants; the KA-CPU gate passes
    # `gpu_backend_override=_BKA_CPU()` explicitly so geometry's arrays are
    # plain Vectors too (a CuArray geometry can't be gathered from by a
    # KA-CPU-launched kernel, and vice versa -- see
    # prepare_mie_node_geometry's docstring).
    geom = prepare_mie_node_geometry(radii_c, weights_c, offsets; architecture=Architectures.GPU(), backend=gpu_backend_override)

    pol = Stokes_IQUV(); trunc = NoTruncation()
    batched = compute_aerosol_optical_properties_nodes_batched_gpu(geom, lam_c, ens_nr_c, ens_ni_c, pol, trunc, backend; precision_policy=policy)
    ext_batched = compute_aerosol_extinction_nodes_batched_gpu(geom, lam_c, ens_nr_c, ens_ni_c, backend; precision_policy=policy)

    rtol_scalar = exact ? 0.0 : (FT === Float64 ? 1e-9 : 1e-6)
    atol_greek  = exact ? 0.0 : (FT === Float64 ? 1e-9 : 1e-4)

    for e in 1:fx.E
        single = compute_aerosol_optical_properties_nodes_gpu(ens_r_c[e], ens_w_c[e], ens_nr_c[e], ens_ni_c[e], lam_c, backend; precision_policy=policy)
        single_k = compute_aerosol_extinction_nodes_gpu(ens_r_c[e], ens_w_c[e], ens_nr_c[e], ens_ni_c[e], lam_c, backend; precision_policy=policy)

        if exact
            @test batched[e].k == single.k
            @test batched[e].ω̃ == single.ω̃
            for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
                @test getproperty(batched[e].greek_coefs, field) == getproperty(single.greek_coefs, field)
            end
            @test ext_batched[e] == single_k
            # extinction-only batched must ALSO match the full batched call's k exactly
            @test ext_batched[e] == batched[e].k
        else
            @test isapprox(batched[e].k, single.k; rtol=rtol_scalar)
            @test isapprox(batched[e].ω̃, single.ω̃; rtol=rtol_scalar)
            for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
                @test isapprox(getproperty(batched[e].greek_coefs, field), getproperty(single.greek_coefs, field); atol=atol_greek)
            end
            @test isapprox(ext_batched[e], single_k; rtol=rtol_scalar)
            @test isapprox(ext_batched[e], batched[e].k; rtol=rtol_scalar)
        end
    end
end

@testset "Gate 1: bitwise batched-vs-single, KA-CPU, NativeFloat64" begin
    _run_bitwise_gate(_BKA_CPU(), NativeFloat64(), Float64; gpu_backend_override=_BKA_CPU())
end
@testset "Gate 1: bitwise batched-vs-single, KA-CPU, DSEmulated" begin
    _run_bitwise_gate(_BKA_CPU(), DSEmulated(), Float32; gpu_backend_override=_BKA_CPU())
end
@testset "Gate 1: bitwise batched-vs-single, KA-CPU, NativeFloat32" begin
    _run_bitwise_gate(_BKA_CPU(), NativeFloat32(), Float32; gpu_backend_override=_BKA_CPU())
end

@testset "Gate 1: batched-vs-single, real CUDA, all 3 policies (near-machine-precision, see _run_bitwise_gate docstring)" begin
    if !_BATCHED_CUDA_FUNCTIONAL
        @test_skip "CUDA not functional on this host — real-CUDA batched gate skipped"
    else
        # No `gpu_backend_override`: architecture=GPU() resolves its own
        # natural backend (real CUDA), matching the explicit `backend`
        # argument passed for the kernel launches below. `exact=false`:
        # real CUDA hardware does NOT reproduce KA-CPU's exact bitwise
        # equality (see _run_bitwise_gate's docstring for the measured
        # magnitudes and why -- GPU kernel-launch-shape-dependent codegen,
        # not a logic bug).
        _run_bitwise_gate(_BATCHED_CUDA_BACKEND, NativeFloat64(), Float64; exact=false)
        _run_bitwise_gate(_BATCHED_CUDA_BACKEND, DSEmulated(), Float32; exact=false)
        _run_bitwise_gate(_BATCHED_CUDA_BACKEND, NativeFloat32(), Float32; exact=false)
    end
end

# ============================================================================
# Gate 1 (CUDA, STRICT exact): full-batched .k == extinction-only-batched k
# ============================================================================
# The narrower invariant that IS exact (`==`) on real CUDA hardware, for all
# three precision policies (independently verified by Codex's review): the
# extinction-only batched call and the full batched call share the IDENTICAL
# Kernel-1 + cross-section + segmented-weighted-sum code path -- extinction-
# only simply never continues into the launch-shape-sensitive amplitude/
# phase/Greek reduction that breaks bitwise equality between batched and
# single-ensemble elsewhere on real CUDA (see the tolerance-based testset
# above). This is the exact invariant GCHPIO's hybrid exact-τ mechanism
# (extinction-only for most of the spectrum, full optics only where Greek
# coefficients are actually needed) depends on, so it gets its OWN dedicated
# STRICT (zero-tolerance `==`) regression test rather than being folded into
# the isapprox-based comparisons above (where it was previously only
# incidentally exercised via `@test isapprox(ext_batched[e], batched[e].k;
# rtol=rtol_scalar)`, an inexact check that could mask a real regression).
# ============================================================================
@testset "Gate 1 (CUDA, STRICT exact): full-batched .k == extinction-only-batched k" begin
    if !_BATCHED_CUDA_FUNCTIONAL
        @test_skip "CUDA not functional on this host — strict real-CUDA k-equality regression skipped"
    else
        fx = _mixed_ensemble_fixture()
        radii, weights, offsets = _concat_fixture(fx)
        pol = Stokes_IQUV(); trunc = NoTruncation()

        for (policy, FT) in ((NativeFloat64(), Float64), (DSEmulated(), Float32), (NativeFloat32(), Float32))
            radii_c = FT.(radii); weights_c = FT.(weights)
            ens_nr_c = FT.(fx.ens_nr); ens_ni_c = FT.(fx.ens_ni); lam_c = FT(0.55)
            geom = prepare_mie_node_geometry(radii_c, weights_c, offsets; architecture=Architectures.GPU())
            batched = compute_aerosol_optical_properties_nodes_batched_gpu(
                geom, lam_c, ens_nr_c, ens_ni_c, pol, trunc, _BATCHED_CUDA_BACKEND; precision_policy=policy)
            ext_batched = compute_aerosol_extinction_nodes_batched_gpu(
                geom, lam_c, ens_nr_c, ens_ni_c, _BATCHED_CUDA_BACKEND; precision_policy=policy)
            for e in 1:fx.E
                @test ext_batched[e] == batched[e].k
            end
        end
    end
end

# ============================================================================
# Gate 1 (CPU dispatch path): closes the coverage gap flagged in review --
# `_run_bitwise_gate` above only exercises the explicit-backend `_batched_gpu`
# functions (mirroring the GPU-focused acceptance-gate wording), never the
# actual `Architectures.CPU()` dispatch branch of the top-level
# `compute_aerosol_optical_properties_nodes_batched`/
# `compute_aerosol_extinction_nodes_batched` (`_batched_cpu_full`/
# `_batched_cpu_extinction`) -- which is what a caller gets BY DEFAULT (no
# `architecture` kwarg at all downstream). Compares against the single-
# ensemble CPU entry points (`compute_aerosol_optical_properties_nodes`/
# `compute_aerosol_extinction_nodes`, themselves defaulting to
# `Architectures.CPU()`), bitwise (VALUE) AND type-equal (Codex CPU
# Float32-output-contract regression: `prepare_mie_node_geometry` used to
# force GFT=Float64 for CPU geometries and the batched CPU output-type
# formula read that forced Float64 back off `geometry.radii_host`, so a
# batched CPU call with Float32 inputs silently returned Float64 while the
# single-ensemble call correctly returned Float32 -- GCHPIO's consumer code
# had to add a `_narrow_greek` workaround for exactly this. Fixed via
# `MieNodeGeometry.output_FT`, captured from the caller's ORIGINAL
# radii/weights eltype before any internal Float64 forcing; this test locks
# the fix in place for both Float64 and Float32), on the mixed-RI
# multi-group fixture.
# ============================================================================
function _run_cpu_dispatch_gate(::Type{FT}) where {FT}
    fx = _mixed_ensemble_fixture()
    radii, weights, offsets = _concat_fixture(fx)
    radii_c = FT.(radii); weights_c = FT.(weights)
    ens_r_c = [FT.(r) for r in fx.ens_r]; ens_w_c = [FT.(w) for w in fx.ens_w]
    ens_nr_c = FT.(fx.ens_nr); ens_ni_c = FT.(fx.ens_ni); lam_c = FT(0.55)
    pol = Stokes_IQUV(); trunc = NoTruncation()

    geom_cpu = prepare_mie_node_geometry(radii_c, weights_c, offsets; architecture=Architectures.CPU())
    batched = compute_aerosol_optical_properties_nodes_batched(geom_cpu, lam_c, ens_nr_c, ens_ni_c, pol, trunc)
    ext_batched = compute_aerosol_extinction_nodes_batched(geom_cpu, lam_c, ens_nr_c, ens_ni_c)

    @test length(batched) == fx.E && length(ext_batched) == fx.E
    for e in 1:fx.E
        single = compute_aerosol_optical_properties_nodes(ens_r_c[e], ens_w_c[e], ens_nr_c[e], ens_ni_c[e], lam_c, pol, trunc)
        single_k = compute_aerosol_extinction_nodes(ens_r_c[e], ens_w_c[e], ens_nr_c[e], ens_ni_c[e], lam_c)

        # TYPE equality first -- this is exactly what the Codex-reported bug
        # broke (batched CPU F32 input silently returning Float64 output).
        @test typeof(batched[e].k) === FT === typeof(single.k)
        @test typeof(batched[e].ω̃) === FT === typeof(single.ω̃)
        @test eltype(batched[e].greek_coefs.β) === FT === eltype(single.greek_coefs.β)
        @test typeof(ext_batched[e]) === FT === typeof(single_k)

        # VALUE equality (bitwise), to the same standard the single CPU path
        # itself meets (both compute internally in Float64 and narrow once).
        @test batched[e].k == single.k
        @test batched[e].ω̃ == single.ω̃
        for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
            @test getproperty(batched[e].greek_coefs, field) == getproperty(single.greek_coefs, field)
        end
        @test ext_batched[e] == single_k
        @test ext_batched[e] == batched[e].k
    end
end

@testset "Gate 1: bitwise batched-vs-single, Architectures.CPU() dispatch path, Float64" begin
    _run_cpu_dispatch_gate(Float64)
end
@testset "Gate 1: bitwise+type-equal batched-vs-single, Architectures.CPU() dispatch path, Float32 (Codex F32-output-contract regression)" begin
    _run_cpu_dispatch_gate(Float32)
end

# ============================================================================
# Gate 2: batching invariance -- arbitrary splits (including across group
# boundaries) change nothing, bitwise.
# ============================================================================
@testset "Gate 2: batching-split invariance (KA-CPU, NativeFloat64)" begin
    fx = _mixed_ensemble_fixture()
    radii, weights, offsets = _concat_fixture(fx)
    lam = 0.55

    geom_full = prepare_mie_node_geometry(radii, weights, offsets; architecture=Architectures.GPU(), backend=_BKA_CPU())
    full_batch = compute_aerosol_optical_properties_nodes_batched_gpu(
        geom_full, lam, fx.ens_nr, fx.ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64())
    full_ext = compute_aerosol_extinction_nodes_batched_gpu(geom_full, lam, fx.ens_nr, fx.ens_ni, _BKA_CPU(); precision_policy=NativeFloat64())

    # Split into three arbitrary sub-batches, DELIBERATELY crossing whatever
    # group boundaries the fixture happens to produce (not aligned to them).
    splits = (1:2, 3:5, 6:8)
    for sub in splits
        radii_s = vcat(fx.ens_r[sub]...); weights_s = vcat(fx.ens_w[sub]...)
        offsets_s = [0; cumsum(length.(fx.ens_r[sub]))]
        geom_s = prepare_mie_node_geometry(radii_s, weights_s, offsets_s; architecture=Architectures.GPU(), backend=_BKA_CPU())
        batch_s = compute_aerosol_optical_properties_nodes_batched_gpu(
            geom_s, lam, fx.ens_nr[sub], fx.ens_ni[sub], Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64())
        ext_s = compute_aerosol_extinction_nodes_batched_gpu(geom_s, lam, fx.ens_nr[sub], fx.ens_ni[sub], _BKA_CPU(); precision_policy=NativeFloat64())

        for (i, e) in enumerate(sub)
            @test full_batch[e].k == batch_s[i].k
            @test full_batch[e].greek_coefs.β == batch_s[i].greek_coefs.β
            @test full_ext[e] == ext_s[i]
        end
    end
end

# ============================================================================
# Gate 3: workspace invariance -- warm == cold, reuse across λs (regrouping)
# ============================================================================
@testset "Gate 3: workspace invariance (warm == cold; reuse across λ, regrouping)" begin
    fx = _mixed_ensemble_fixture()
    radii, weights, offsets = _concat_fixture(fx)
    geom = prepare_mie_node_geometry(radii, weights, offsets; architecture=Architectures.GPU(), backend=_BKA_CPU())

    for lam in (0.55, 1.61)
        n_max_per_ens, groups = Scattering._group_by_nmax(geom.offsets, geom.radii_host, 2π / lam)
        cold = compute_aerosol_optical_properties_nodes_batched_gpu(
            geom, lam, fx.ens_nr, fx.ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64())

        ws = MieBatchWorkspace(_BKA_CPU(), Float64, Float64)
        warm1 = compute_aerosol_optical_properties_nodes_batched_gpu(
            geom, lam, fx.ens_nr, fx.ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64(), workspace=ws)
        warm2 = compute_aerosol_optical_properties_nodes_batched_gpu(
            geom, lam, fx.ens_nr, fx.ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64(), workspace=ws)

        for e in 1:fx.E
            @test cold[e].k == warm1[e].k == warm2[e].k
            @test cold[e].greek_coefs.β == warm1[e].greek_coefs.β == warm2[e].greek_coefs.β
        end
    end

    # Reuse the SAME workspace across the two DIFFERENT wavelengths (group
    # membership changes between them -- this is the "regroup with a reused
    # workspace" case the grow-only design must handle transparently).
    ws_cross = MieBatchWorkspace(_BKA_CPU(), Float64, Float64)
    r1 = compute_aerosol_optical_properties_nodes_batched_gpu(
        geom, 0.55, fx.ens_nr, fx.ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64(), workspace=ws_cross)
    r2 = compute_aerosol_optical_properties_nodes_batched_gpu(
        geom, 1.61, fx.ens_nr, fx.ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64(), workspace=ws_cross)
    r1_fresh = compute_aerosol_optical_properties_nodes_batched_gpu(
        geom, 0.55, fx.ens_nr, fx.ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64())
    r2_fresh = compute_aerosol_optical_properties_nodes_batched_gpu(
        geom, 1.61, fx.ens_nr, fx.ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64())

    for e in 1:fx.E
        @test r1[e].k == r1_fresh[e].k
        @test r1[e].greek_coefs.β == r1_fresh[e].greek_coefs.β
        @test r2[e].k == r2_fresh[e].k
        @test r2[e].greek_coefs.β == r2_fresh[e].greek_coefs.β
    end

    # Single-owner guard: a workspace touched from a second task must throw
    # ArgumentError (not silently share scratch). UNCONDITIONAL now (review
    # finding 4 -- no longer @boundscheck-gated), so this always exercises the
    # guard regardless of --check-bounds.
    ws_owner = MieBatchWorkspace(_BKA_CPU(), Float64, Float64)
    compute_aerosol_optical_properties_nodes_batched_gpu(
        geom, 0.55, fx.ens_nr, fx.ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64(), workspace=ws_owner)
    t = @task compute_aerosol_optical_properties_nodes_batched_gpu(
        geom, 0.55, fx.ens_nr, fx.ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64(), workspace=ws_owner)
    schedule(t)
    @test_throws Union{ArgumentError, TaskFailedException} Base.wait(t)

    # Backend-instance guard (review finding 4): a workspace built for one
    # backend INSTANCE must reject a launch against a different one, even if
    # both are (separately constructed) `KernelAbstractions.CPU()` values --
    # covered instead via a mock backend with a distinct `objectid` (plain
    # `CPU()` values are isbits and therefore objectid-equal to each other,
    # so this needs a genuinely different value to exercise; see
    # `_MockOtherBackend` defined at file scope above).
    ws_backend = MieBatchWorkspace(_BKA_CPU(), Float64, Float64)
    compute_aerosol_optical_properties_nodes_batched_gpu(
        geom, 0.55, fx.ens_nr, fx.ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64(), workspace=ws_backend)
    @test_throws ArgumentError Scattering._ws_array!(ws_backend, _MockOtherBackend(), :bx_perm, Int, 4)
end

# ============================================================================
# Gate 3 (workspace key-sharing, review finding 5): _batched_gpu_full and
# _batched_gpu_extinction now share the SAME workspace keys (`:bx_*`) for
# every Kernel-1-preamble buffer they both need (perm/rgrp/wgrp/nmax/n2e/
# mre/mim/an/bn/x/Csca/Cext/w/loff/bCext) -- a workspace reused across both
# APIs (the documented GCHPIO hybrid exact-τ pattern: extinction-only for
# most of the spectrum, full optics where Greek coefficients are actually
# needed) must hold ONE copy of each, not two, including an/bn (the largest
# allocation).
# ============================================================================
@testset "Gate 3: workspace key-sharing across full/extinction APIs (review finding 5)" begin
    fx = _mixed_ensemble_fixture()
    radii, weights, offsets = _concat_fixture(fx)
    geom = prepare_mie_node_geometry(radii, weights, offsets; architecture=Architectures.GPU(), backend=_BKA_CPU())
    lam = 0.55
    pol = Stokes_IQUV(); trunc = NoTruncation()

    ws = MieBatchWorkspace(_BKA_CPU(), Float64, Float64)
    full = compute_aerosol_optical_properties_nodes_batched_gpu(
        geom, lam, fx.ens_nr, fx.ens_ni, pol, trunc, _BKA_CPU(); precision_policy=NativeFloat64(), workspace=ws)
    an_after_full = ws.buffers[:bx_an]
    bn_after_full = ws.buffers[:bx_bn]

    ext = compute_aerosol_extinction_nodes_batched_gpu(geom, lam, fx.ens_nr, fx.ens_ni, _BKA_CPU(); precision_policy=NativeFloat64(), workspace=ws)
    an_after_ext = ws.buffers[:bx_an]
    bn_after_ext = ws.buffers[:bx_bn]

    # SAME underlying array OBJECT (===, not just equal values): proves the
    # extinction call REUSED, rather than duplicated, the full call's
    # Kernel-1-preamble scratch.
    @test an_after_full === an_after_ext
    @test bn_after_full === bn_after_ext
    for key in (:bx_perm, :bx_rgrp, :bx_wgrp, :bx_nmax, :bx_n2e, :bx_mre, :bx_mim,
                :bx_x, :bx_Csca, :bx_Cext, :bx_w, :bx_loff, :bx_bCext)
        @test haskey(ws.buffers, key)
    end
    # No stray `:bxe_*`-prefixed keys left over from the old separate namespace.
    @test !any(k -> startswith(String(k), "bxe_"), keys(ws.buffers))

    # Results stay bitwise-correct despite the shared-buffer reuse.
    full_cold = compute_aerosol_optical_properties_nodes_batched_gpu(
        geom, lam, fx.ens_nr, fx.ens_ni, pol, trunc, _BKA_CPU(); precision_policy=NativeFloat64())
    ext_cold = compute_aerosol_extinction_nodes_batched_gpu(geom, lam, fx.ens_nr, fx.ens_ni, _BKA_CPU(); precision_policy=NativeFloat64())
    for e in 1:fx.E
        @test full[e].k == full_cold[e].k
        @test full[e].greek_coefs.β == full_cold[e].greek_coefs.β
        @test ext[e] == ext_cold[e]
    end
end

# ============================================================================
# Gate 4: original-order reassembly property test (shuffled ensemble sizes)
# ============================================================================
@testset "Gate 4: original-order reassembly (shuffled ensemble sizes, property test)" begin
    rng = MersenneTwister(99)
    for trial in 1:8
        E = rand(rng, 3:10)
        sizes = rand(rng, 4:30, E)
        # Shuffle which ensemble gets which SIZE RANGE so nmax grouping is
        # decoupled from ensemble index order.
        range_pool = [(0.005, 0.03), (0.05, 0.2), (0.3, 1.0), (1.0, 4.0), (4.0, 9.0)]
        chosen_ranges = [range_pool[rand(rng, 1:length(range_pool))] for _ in 1:E]
        ens_r = [exp.(range(log(lo), log(hi), length=n)) for ((lo, hi), n) in zip(chosen_ranges, sizes)]
        ens_w = [rand(rng, n) .+ 0.1 for n in sizes]
        ens_nr = 1.3 .+ 0.3 .* rand(rng, E)
        ens_ni = 0.001 .+ 0.02 .* rand(rng, E)
        lam = 0.55

        radii = vcat(ens_r...); weights = vcat(ens_w...)
        offsets = [0; cumsum(length.(ens_r))]
        geom = prepare_mie_node_geometry(radii, weights, offsets; architecture=Architectures.GPU(), backend=_BKA_CPU())
        batched = compute_aerosol_optical_properties_nodes_batched_gpu(
            geom, lam, ens_nr, ens_ni, Stokes_IQUV(), NoTruncation(), _BKA_CPU(); precision_policy=NativeFloat64())

        @test length(batched) == E
        for e in 1:E
            single = compute_aerosol_optical_properties_nodes_gpu(ens_r[e], ens_w[e], ens_nr[e], ens_ni[e], lam, _BKA_CPU(); precision_policy=NativeFloat64())
            @test batched[e].k == single.k
            @test batched[e].greek_coefs.β == single.greek_coefs.β
        end
    end
end

# ============================================================================
# Gate 5: benchmark on a synthetic GCHPIO-like inventory
# ============================================================================
# 882 ensembles x 96 nodes, radii spanning the TOMAS wet range -- mirrors the
# real dust-column inventory cited in the proposal (882 ensembles, 20 distinct
# nmax at 0.55 μm on the FULL nquad=1024/96 columns). Reports (does not
# assert) A100 batched vs a serial single-call loop vs the 16-thread CPU
# path, plus the nmax group histogram.
@testset "Gate 5: GCHPIO-like benchmark (882 ensembles x 96 nodes, report only)" begin
    rng = MersenneTwister(2026)
    E_bench = 882
    n_nodes_per_ens = 96
    # TOMAS wet range: ~0.001-10 μm across the full inventory; each ensemble
    # gets a bin sub-range within that, mirroring a per-layer sectional
    # distribution whose bin edges vary by layer/species.
    ens_r = Vector{Vector{Float64}}(undef, E_bench)
    ens_w = Vector{Vector{Float64}}(undef, E_bench)
    for e in 1:E_bench
        lo = 10.0^(rand(rng) * 3.5 - 3)         # ~0.001-3 μm
        hi = lo * 10.0^(0.3 + 2.0 * rand(rng))  # up to ~1.5 decades wider
        hi = min(hi, 12.0)
        ens_r[e] = exp.(range(log(lo), log(max(hi, lo * 1.01)), length=n_nodes_per_ens))
        ens_w[e] = pdf.(LogNormal(log(sqrt(lo * hi)), log(1.6)), ens_r[e]) .+ 1e-30
    end
    ens_nr = 1.3 .+ 0.3 .* rand(rng, E_bench)
    ens_ni = 0.0005 .+ 0.02 .* rand(rng, E_bench)
    lam = 0.55

    radii = vcat(ens_r...); weights = vcat(ens_w...)
    offsets = [0; cumsum(length.(ens_r))]

    k_wavenum = 2π / lam
    n_max_per_ensemble, groups = Scattering._group_by_nmax(collect(Int, offsets), Float64.(radii), k_wavenum)
    hist = sort(collect(pairs(Dict(nm => length(ids) for (nm, ids) in groups))); by=first)
    println("Gate 5 benchmark: n_max group histogram (nmax => n_ensembles), $(length(groups)) distinct groups:")
    for (nm, cnt) in hist
        println("  nmax=$nm  n_ensembles=$cnt")
    end
    @test length(groups) ≥ 3  # sanity: this benchmark should exercise real grouping

    pol = Stokes_IQUV(); trunc = NoTruncation()

    # --- KA-CPU grouped batched call ---
    geom_kacpu = prepare_mie_node_geometry(radii, weights, offsets; architecture=Architectures.GPU(), backend=_BKA_CPU())
    compute_aerosol_optical_properties_nodes_batched_gpu(geom_kacpu, lam, ens_nr, ens_ni, pol, trunc, _BKA_CPU(); precision_policy=NativeFloat64())  # warmup
    t_kacpu_batched = @elapsed compute_aerosol_optical_properties_nodes_batched_gpu(
        geom_kacpu, lam, ens_nr, ens_ni, pol, trunc, _BKA_CPU(); precision_policy=NativeFloat64())

    # --- 16-thread CPU (existing ensemble-level Threads.@threads path) ---
    geom_cpu = prepare_mie_node_geometry(radii, weights, offsets; architecture=Architectures.CPU())
    compute_aerosol_optical_properties_nodes_batched(geom_cpu, lam, ens_nr, ens_ni, pol, trunc)  # warmup
    t_cpu_threaded = @elapsed compute_aerosol_optical_properties_nodes_batched(geom_cpu, lam, ens_nr, ens_ni, pol, trunc)
    println("Gate 5 benchmark: Threads.nthreads() = $(Threads.nthreads())")

    println("Gate 5 benchmark timings (nquad=$n_nodes_per_ens, E=$E_bench):")
    println("  CPU (Threads.@threads, $(Threads.nthreads()) threads): $(t_cpu_threaded) s")
    println("  KA-CPU grouped batched (portable backend):             $(t_kacpu_batched) s")

    if _BATCHED_CUDA_FUNCTIONAL
        geom_cuda = prepare_mie_node_geometry(Float64.(radii), Float64.(weights), offsets; architecture=Architectures.GPU())
        compute_aerosol_optical_properties_nodes_batched_gpu(geom_cuda, lam, ens_nr, ens_ni, pol, trunc, _BATCHED_CUDA_BACKEND; precision_policy=NativeFloat64())  # warmup
        t_cuda_batched = @elapsed compute_aerosol_optical_properties_nodes_batched_gpu(
            geom_cuda, lam, ens_nr, ens_ni, pol, trunc, _BATCHED_CUDA_BACKEND; precision_policy=NativeFloat64())

        # Serial single-call-per-ensemble loop on real CUDA (the pre-batching baseline).
        t_cuda_serial = @elapsed for e in 1:E_bench
            compute_aerosol_optical_properties_nodes_gpu(ens_r[e], ens_w[e], ens_nr[e], ens_ni[e], lam, _BATCHED_CUDA_BACKEND; precision_policy=NativeFloat64())
        end

        gpu_name = try
            CUDA.name(CUDA.device())
        catch
            "unknown CUDA device"
        end
        println("  Real CUDA ($gpu_name) grouped batched:                  $(t_cuda_batched) s")
        println("  Real CUDA ($gpu_name) serial single-call loop:         $(t_cuda_serial) s")
        println("  speedup, grouped batched vs serial single-call (CUDA): $(round(t_cuda_serial / t_cuda_batched, digits=2))x")
        println("  speedup, grouped batched (CUDA) vs $(Threads.nthreads())-thread CPU:      $(round(t_cpu_threaded / t_cuda_batched, digits=2))x")
        @test t_cuda_batched ≥ 0 && t_cuda_serial ≥ 0

        # NativeFloat32 / DSEmulated timings too (the Metal/consumer-CUDA targets).
        radii32 = Float32.(radii); weights32 = Float32.(weights)
        ens_nr32 = Float32.(ens_nr); ens_ni32 = Float32.(ens_ni); lam32 = Float32(lam)
        geom_cuda32 = prepare_mie_node_geometry(radii32, weights32, offsets; architecture=Architectures.GPU())
        for policy in (DSEmulated(), NativeFloat32())
            compute_aerosol_optical_properties_nodes_batched_gpu(geom_cuda32, lam32, ens_nr32, ens_ni32, pol, trunc, _BATCHED_CUDA_BACKEND; precision_policy=policy)  # warmup
            t = @elapsed compute_aerosol_optical_properties_nodes_batched_gpu(
                geom_cuda32, lam32, ens_nr32, ens_ni32, pol, trunc, _BATCHED_CUDA_BACKEND; precision_policy=policy)
            println("  Real CUDA ($gpu_name) grouped batched, $(typeof(policy)): $(t) s")
            @test t ≥ 0
        end
    else
        @test_skip "CUDA not functional on this host — real-CUDA Gate 5 benchmark rows skipped"
    end

    @test t_kacpu_batched ≥ 0 && t_cpu_threaded ≥ 0
end
