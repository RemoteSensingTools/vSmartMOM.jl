#!/usr/bin/env julia
#
# Benchmark: full interaction_11 with algebraic reduction + bundled solve kernel
#
# This benchmark compares three versions of the forward ScatteringInterface_11
# update from src/CoreRT/CoreKernel/interaction.jl:
#
#   1. Current-style explicit-inverse path:
#        two factorizations/inverses, then chained batched multiplies
#
#   2. Algebraically reduced cuBLAS path:
#        one factorization/inverse of M = I - R⁺⁻r⁻⁺, using the identities
#          (I - AB)^(-1) A = A (I - BA)^(-1)
#          (I - AB)^(-1) = I + A (I - BA)^(-1) B
#
#      With A = r⁻⁺ and B = R⁺⁻, the full interaction_11 update becomes:
#
#        X_T = M \ T⁺⁺
#        X_R = M \ (R⁺⁻ t⁻⁻)
#        x_J = M \ (J₀⁺ + R⁺⁻ j₀⁻)
#
#        T⁺⁺_out = t⁺⁺ X_T
#        R⁺⁻_out = r⁺⁻ + t⁺⁺ X_R
#        J₀⁺_out = j₀⁺ + t⁺⁺ x_J
#
#        P = T⁻⁻ r⁻⁺
#        T⁻⁻_out = T⁻⁻ t⁻⁻ + P X_R
#        R⁻⁺_out = R⁻⁺ + P X_T
#        J₀⁻_out = J₀⁻ + T⁻⁻ j₀⁻ + P x_J
#
#   3. Bundled KA kernel:
#        same reduced algebra as (2), but factorization + three solves + all
#        output updates happen inside one kernel launch per batch element.
#
# Usage:
#   julia --project=. test/benchmarks/interaction11_bundled_benchmark.jl
#

using CUDA
using KernelAbstractions
using LinearAlgebra
using Printf

if !CUDA.functional()
    println("CUDA is not functional on this machine; skipping interaction_11 benchmark.")
    exit()
end

const FORCE_PEDANTIC = get(ENV, "VSMARTMOM_ALLOW_TF32", "0") != "1"
const MEMORY_FRACTION_LIMIT = parse(Float64, get(ENV, "VSMARTMOM_BENCH_MEMORY_FRACTION", "0.80"))

# ============================================================================
# Bundled KA kernel using the algebraically reduced formulation
# ============================================================================
@kernel function interaction11_reduced_bundled_kernel!(
    R_dn_out, T_dn_out, R_up_out, T_up_out, J_dn_out, J_up_out,
    pref,
    @Const(R_dn), @Const(R_up), @Const(T_up), @Const(T_dn), @Const(J_up), @Const(J_dn),
    @Const(r_up), @Const(r_dn), @Const(t_up), @Const(t_dn), @Const(j_up), @Const(j_dn),
    ::Val{N}
) where {N}
    k = @index(Group, Linear)
    tid = @index(Local, Linear)

    LU   = @localmem eltype(R_dn_out) (N, N)
    work = @localmem eltype(R_dn_out) (N, N)
    piv  = @localmem Int32 (N,)
    vec  = @localmem eltype(R_dn_out) (N,)

    # Precompute:
    #   pref = T⁻⁻ r⁻⁺
    #   T_dn_out = T⁻⁻ t⁻⁻
    #   J_dn_out = J₀⁻ + T⁻⁻ j₀⁻
    @inbounds begin
        for j in 1:N
            s_pref = zero(eltype(R_dn_out))
            s_tdn = zero(eltype(R_dn_out))
            for l in 1:N
                s_pref += T_dn[tid, l, k] * r_dn[l, j, k]
                s_tdn += T_dn[tid, l, k] * t_dn[l, j, k]
            end
            pref[tid, j, k] = s_pref
            T_dn_out[tid, j, k] = s_tdn
        end

        s_jdn = J_dn[tid, 1, k]
        for l in 1:N
            s_jdn += T_dn[tid, l, k] * j_dn[l, 1, k]
        end
        J_dn_out[tid, 1, k] = s_jdn
    end
    @synchronize()

    # Factor M = I - R⁺⁻ r⁻⁺
    @inbounds for j in 1:N
        s = zero(eltype(R_dn_out))
        for l in 1:N
            s += R_up[tid, l, k] * r_dn[l, j, k]
        end
        LU[tid, j] = ((tid == j) ? one(eltype(R_dn_out)) : zero(eltype(R_dn_out))) - s
    end
    @inbounds piv[tid] = Int32(tid)
    @synchronize()

    @inbounds for p in 1:N
        if tid == 1
            max_val = abs(LU[p, p])
            max_row = p
            for r in (p + 1):N
                v = abs(LU[r, p])
                if v > max_val
                    max_val = v
                    max_row = r
                end
            end
            if max_row != p
                for j in 1:N
                    tmp = LU[p, j]
                    LU[p, j] = LU[max_row, j]
                    LU[max_row, j] = tmp
                end
                tmp_p = piv[p]
                piv[p] = piv[max_row]
                piv[max_row] = tmp_p
            end
        end
        @synchronize()

        if tid > p
            LU[tid, p] /= LU[p, p]
        end
        @synchronize()

        if tid > p
            factor = LU[tid, p]
            for j in (p + 1):N
                LU[tid, j] -= factor * LU[p, j]
            end
        end
        @synchronize()
    end

    # -------------------------------------------------------------------------
    # Solve X_T = M \ T⁺⁺
    # -------------------------------------------------------------------------
    @inbounds for i in 1:N
        row = piv[i]
        work[i, tid] = T_up[row, tid, k]
    end
    @synchronize()

    @inbounds for i in 2:N
        s = zero(eltype(R_dn_out))
        for j in 1:(i - 1)
            s += LU[i, j] * work[j, tid]
        end
        work[i, tid] -= s
    end

    @inbounds for i in N:-1:1
        s = zero(eltype(R_dn_out))
        for j in (i + 1):N
            s += LU[i, j] * work[j, tid]
        end
        work[i, tid] = (work[i, tid] - s) / LU[i, i]
    end
    @synchronize()

    @inbounds for j in 1:N
        s_up = zero(eltype(R_dn_out))
        s_dn = zero(eltype(R_dn_out))
        for l in 1:N
            s_up += t_up[tid, l, k] * work[l, j]
            s_dn += pref[tid, l, k] * work[l, j]
        end
        T_up_out[tid, j, k] = s_up
        R_dn_out[tid, j, k] = R_dn[tid, j, k] + s_dn
    end
    @synchronize()

    # -------------------------------------------------------------------------
    # Solve X_R = M \ (R⁺⁻ t⁻⁻)
    # -------------------------------------------------------------------------
    @inbounds for i in 1:N
        row = piv[i]
        s = zero(eltype(R_dn_out))
        for l in 1:N
            s += R_up[row, l, k] * t_dn[l, tid, k]
        end
        work[i, tid] = s
    end
    @synchronize()

    @inbounds for i in 2:N
        s = zero(eltype(R_dn_out))
        for j in 1:(i - 1)
            s += LU[i, j] * work[j, tid]
        end
        work[i, tid] -= s
    end

    @inbounds for i in N:-1:1
        s = zero(eltype(R_dn_out))
        for j in (i + 1):N
            s += LU[i, j] * work[j, tid]
        end
        work[i, tid] = (work[i, tid] - s) / LU[i, i]
    end
    @synchronize()

    @inbounds for j in 1:N
        s_up = zero(eltype(R_dn_out))
        s_dn = zero(eltype(R_dn_out))
        for l in 1:N
            s_up += t_up[tid, l, k] * work[l, j]
            s_dn += pref[tid, l, k] * work[l, j]
        end
        R_up_out[tid, j, k] = r_up[tid, j, k] + s_up
        T_dn_out[tid, j, k] += s_dn
    end
    @synchronize()

    # -------------------------------------------------------------------------
    # Solve x_J = M \ (J⁺ + R⁺⁻ j⁻)
    # -------------------------------------------------------------------------
    if tid == 1
        for i in 1:N
            row = piv[i]
            s = J_up[row, 1, k]
            for l in 1:N
                s += R_up[row, l, k] * j_dn[l, 1, k]
            end
            vec[i] = s
        end

        for i in 2:N
            s = zero(eltype(R_dn_out))
            for j in 1:(i - 1)
                s += LU[i, j] * vec[j]
            end
            vec[i] -= s
        end

        for i in N:-1:1
            s = zero(eltype(R_dn_out))
            for j in (i + 1):N
                s += LU[i, j] * vec[j]
            end
            vec[i] = (vec[i] - s) / LU[i, i]
        end
    end
    @synchronize()

    @inbounds begin
        s_up = zero(eltype(R_dn_out))
        s_dn = zero(eltype(R_dn_out))
        for l in 1:N
            s_up += t_up[tid, l, k] * vec[l]
            s_dn += pref[tid, l, k] * vec[l]
        end
        J_up_out[tid, 1, k] = j_up[tid, 1, k] + s_up
        J_dn_out[tid, 1, k] += s_dn
    end
end

function ka_interaction11_bundled!(outs, ins, scratch)
    N = size(ins.R_dn, 1)
    batch = size(ins.R_dn, 3)
    backend = CUDABackend()

    interaction11_reduced_bundled_kernel!(backend, N)(
        outs.R_dn, outs.T_dn, outs.R_up, outs.T_up, outs.J_dn, outs.J_up,
        scratch.pref,
        ins.R_dn, ins.R_up, ins.T_up, ins.T_dn, ins.J_up, ins.J_dn,
        ins.r_up, ins.r_dn, ins.t_up, ins.t_dn, ins.j_up, ins.j_dn,
        Val(N);
        ndrange=(N * batch,),
    )
    KernelAbstractions.synchronize(backend)
    return outs
end

# ============================================================================
# cuBLAS reference paths
# ============================================================================
function cublas_batched_mul!(C::CuArray{FT,3}, A::CuArray{FT,3}, B::CuArray{FT,3}) where {FT}
    CUDA.CUBLAS.gemm_strided_batched!('N', 'N', one(FT), A, B, zero(FT), C)
    return C
end

function cublas_batch_inv!(X::CuArray{FT,3}, A::CuArray{FT,3}) where {FT}
    pivot, info = CUDA.CUBLAS.getrf_strided_batched!(A, true)
    CUDA.synchronize()
    CUDA.CUBLAS.getri_strided_batched!(A, X, pivot)
    CUDA.synchronize()
    return X
end

function cublas_current_interaction11!(outs, ins, scratch)
    # Direction 1: explicit inverse of I - r⁻⁺R⁺⁻
    cublas_batched_mul!(scratch.Mbuf, ins.r_dn, ins.R_up)
    scratch.Mbuf .= ins.I_static .- scratch.Mbuf
    cublas_batch_inv!(scratch.invM, scratch.Mbuf)
    cublas_batched_mul!(scratch.solve_prefac, ins.T_dn, scratch.invM)

    cublas_batched_mul!(scratch.matbuf, ins.r_dn, ins.T_up)
    cublas_batched_mul!(outs.R_dn, scratch.solve_prefac, scratch.matbuf)
    outs.R_dn .+= ins.R_dn

    cublas_batched_mul!(outs.T_dn, scratch.solve_prefac, ins.t_dn)

    cublas_batched_mul!(scratch.vecbuf, ins.r_dn, ins.J_up)
    scratch.vecbuf .+= ins.j_dn
    cublas_batched_mul!(outs.J_dn, scratch.solve_prefac, scratch.vecbuf)
    outs.J_dn .+= ins.J_dn

    # Direction 2: explicit inverse of I - R⁺⁻r⁻⁺
    cublas_batched_mul!(scratch.Mbuf, ins.R_up, ins.r_dn)
    scratch.Mbuf .= ins.I_static .- scratch.Mbuf
    cublas_batch_inv!(scratch.invM, scratch.Mbuf)
    cublas_batched_mul!(scratch.solve_prefac, ins.t_up, scratch.invM)

    cublas_batched_mul!(outs.T_up, scratch.solve_prefac, ins.T_up)

    cublas_batched_mul!(scratch.matbuf, ins.R_up, ins.t_dn)
    cublas_batched_mul!(outs.R_up, scratch.solve_prefac, scratch.matbuf)
    outs.R_up .+= ins.r_up

    cublas_batched_mul!(scratch.vecbuf, ins.R_up, ins.j_dn)
    scratch.vecbuf .+= ins.J_up
    cublas_batched_mul!(outs.J_up, scratch.solve_prefac, scratch.vecbuf)
    outs.J_up .+= ins.j_up

    return outs
end

function cublas_reduced_interaction11!(outs, ins, scratch)
    # pref = T⁻⁻ r⁻⁺
    cublas_batched_mul!(scratch.pref, ins.T_dn, ins.r_dn)

    # Base terms for reduced lower-direction outputs
    cublas_batched_mul!(outs.T_dn, ins.T_dn, ins.t_dn)
    cublas_batched_mul!(outs.J_dn, ins.T_dn, ins.j_dn)
    outs.J_dn .+= ins.J_dn

    # One inverse only: M = I - R⁺⁻r⁻⁺
    cublas_batched_mul!(scratch.Mbuf, ins.R_up, ins.r_dn)
    scratch.Mbuf .= ins.I_static .- scratch.Mbuf
    cublas_batch_inv!(scratch.invM, scratch.Mbuf)

    # X_T = inv(M) * T⁺⁺
    cublas_batched_mul!(scratch.solve_mat, scratch.invM, ins.T_up)
    cublas_batched_mul!(outs.T_up, ins.t_up, scratch.solve_mat)
    cublas_batched_mul!(outs.R_dn, scratch.pref, scratch.solve_mat)
    outs.R_dn .+= ins.R_dn

    # X_R = inv(M) * (R⁺⁻ t⁻⁻)
    cublas_batched_mul!(scratch.matbuf, ins.R_up, ins.t_dn)
    cublas_batched_mul!(scratch.solve_mat, scratch.invM, scratch.matbuf)
    cublas_batched_mul!(outs.R_up, ins.t_up, scratch.solve_mat)
    outs.R_up .+= ins.r_up
    cublas_batched_mul!(scratch.accum_mat, scratch.pref, scratch.solve_mat)
    outs.T_dn .+= scratch.accum_mat

    # x_J = inv(M) * (J⁺ + R⁺⁻ j⁻)
    cublas_batched_mul!(scratch.vecbuf, ins.R_up, ins.j_dn)
    scratch.vecbuf .+= ins.J_up
    cublas_batched_mul!(scratch.solve_vec, scratch.invM, scratch.vecbuf)
    cublas_batched_mul!(outs.J_up, ins.t_up, scratch.solve_vec)
    outs.J_up .+= ins.j_up
    cublas_batched_mul!(scratch.accum_vec, scratch.pref, scratch.solve_vec)
    outs.J_dn .+= scratch.accum_vec

    return outs
end

# ============================================================================
# Benchmark utilities
# ============================================================================
function make_identity_batch(FT, N, batch)
    I_cpu = repeat(reshape(Matrix{FT}(I, N, N), N, N, 1), 1, 1, batch)
    return CuArray(I_cpu)
end

function make_inputs(FT, N, batch)
    I_static = make_identity_batch(FT, N, batch)

    return (
        # Composite layer
        R_dn = FT(0.08) .* CUDA.rand(FT, N, N, batch),
        R_up = FT(0.08) .* CUDA.rand(FT, N, N, batch),
        T_up = I_static .+ FT(0.04) .* CUDA.rand(FT, N, N, batch),
        T_dn = I_static .+ FT(0.04) .* CUDA.rand(FT, N, N, batch),
        J_up = FT(0.1) .* CUDA.rand(FT, N, 1, batch),
        J_dn = FT(0.1) .* CUDA.rand(FT, N, 1, batch),
        # Added layer
        r_up = FT(0.08) .* CUDA.rand(FT, N, N, batch),
        r_dn = FT(0.08) .* CUDA.rand(FT, N, N, batch),
        t_up = I_static .+ FT(0.04) .* CUDA.rand(FT, N, N, batch),
        t_dn = I_static .+ FT(0.04) .* CUDA.rand(FT, N, N, batch),
        j_up = FT(0.1) .* CUDA.rand(FT, N, 1, batch),
        j_dn = FT(0.1) .* CUDA.rand(FT, N, 1, batch),
        I_static = I_static,
    )
end

function make_outputs(ins)
    return (
        R_dn = similar(ins.R_dn),
        T_dn = similar(ins.T_dn),
        R_up = similar(ins.R_up),
        T_up = similar(ins.T_up),
        J_dn = similar(ins.J_dn),
        J_up = similar(ins.J_up),
    )
end

function make_current_scratch(ins)
    return (
        Mbuf = similar(ins.R_dn),
        invM = similar(ins.R_dn),
        solve_prefac = similar(ins.R_dn),
        matbuf = similar(ins.R_dn),
        vecbuf = similar(ins.J_dn),
    )
end

function make_reduced_scratch(ins)
    return (
        pref = similar(ins.R_dn),
        Mbuf = similar(ins.R_dn),
        invM = similar(ins.R_dn),
        matbuf = similar(ins.R_dn),
        solve_mat = similar(ins.R_dn),
        accum_mat = similar(ins.R_dn),
        vecbuf = similar(ins.J_dn),
        solve_vec = similar(ins.J_dn),
        accum_vec = similar(ins.J_dn),
    )
end

function make_bundled_scratch(ins)
    return (pref = similar(ins.R_dn),)
end

function time_interaction(f!, outs, ins, scratch; nwarmup=3, nruns=20)
    for _ in 1:nwarmup
        f!(outs, ins, scratch)
    end
    CUDA.synchronize()

    times = Float64[]
    for _ in 1:nruns
        CUDA.synchronize()
        t = CUDA.@elapsed f!(outs, ins, scratch)
        push!(times, t)
    end
    return sort(times)[cld(nruns, 2)]
end

function max_absdiff_sample(A::CuArray{T,3}, B::CuArray{T,3}; ncheck=8) where {T}
    nb = min(size(A, 3), ncheck)
    A_h = Array(@view A[:, :, 1:nb])
    B_h = Array(@view B[:, :, 1:nb])
    return maximum(abs.(A_h .- B_h))
end

function max_output_diff(ref_outs, test_outs; ncheck=8)
    return maximum((
        max_absdiff_sample(ref_outs.R_dn, test_outs.R_dn; ncheck=ncheck),
        max_absdiff_sample(ref_outs.T_dn, test_outs.T_dn; ncheck=ncheck),
        max_absdiff_sample(ref_outs.R_up, test_outs.R_up; ncheck=ncheck),
        max_absdiff_sample(ref_outs.T_up, test_outs.T_up; ncheck=ncheck),
        max_absdiff_sample(ref_outs.J_dn, test_outs.J_dn; ncheck=ncheck),
        max_absdiff_sample(ref_outs.J_up, test_outs.J_up; ncheck=ncheck),
    ))
end

function free_namedtuple!(nt)
    for arr in values(nt)
        CUDA.unsafe_free!(arr)
    end
end

function bundled_smem_bytes(N, FT)
    return 2 * N * N * sizeof(FT) + N * sizeof(FT) + N * sizeof(Int32)
end

function estimate_case_bytes(N, batch, FT)
    matrix_elems = N * N * batch
    vector_elems = N * batch

    # Simultaneously allocated per case:
    # inputs   ->  9 matrices + 4 vectors
    # outputs  -> 12 matrices + 6 vectors   (3 output tuples)
    # scratch  -> 11 matrices + 4 vectors   (current + reduced + bundled)
    nmats = 32
    nvecs = 14
    return (nmats * matrix_elems + nvecs * vector_elems) * sizeof(FT)
end

gib(bytes) = bytes / 2.0^30

function should_skip_case(N, batch, FT)
    GC.gc()
    CUDA.reclaim()
    est = estimate_case_bytes(N, batch, FT)
    free = CUDA.free_memory()
    return est > MEMORY_FRACTION_LIMIT * free, est, free
end

# ============================================================================
# Main benchmark
# ============================================================================
function run_benchmark()
    matrix_sizes = [4, 8, 12, 16, 24, 32, 48]
    batch_sizes = [500, 2000, 5000, 10000, 50000]
    float_types = [Float32, Float64]

    starting_math_mode = CUDA.math_mode()
    if FORCE_PEDANTIC
        CUDA.math_mode!(CUDA.PEDANTIC_MATH)
    end

    println("=" ^ 164)
    println("Full interaction_11 Benchmark: current vs reduced-algebra vs bundled solve kernel")
    println("GPU: ", CUDA.name(CUDA.device()))
    println("cuBLAS math mode: ", CUDA.math_mode(), FORCE_PEDANTIC ? " (forced pedantic; TF32 disabled)" : " (default; TF32 may be used for Float32)")
    println("Memory guard: skip if estimated case footprint exceeds ", round(100 * MEMORY_FRACTION_LIMIT; digits=0), "% of current free GPU memory")
    println("=" ^ 164)

    try
        for FT in float_types
            println("\n", "-" ^ 164)
            @printf("%-8s | %4s | %6s | %12s | %12s | %12s | %8s | %8s | %10s | %10s\n",
                    "Type", "N", "Batch",
                    "Current ms", "Reduced ms", "Bundled ms",
                    "Red spd", "Bnd spd", "Max diff", "SMem KB")
            println("-" ^ 164)

            for N in matrix_sizes
                for batch in batch_sizes
                    skip, est, free = should_skip_case(N, batch, FT)
                    if skip
                        @printf("%-8s | %4d | %6d | %12s | %12s | %12s | %8s | %8s | %10s | %10.1f\n",
                                string(FT), N, batch,
                                "SKIP(mem)", "SKIP(mem)", "SKIP(mem)",
                                "-", "-", "-", bundled_smem_bytes(N, FT) / 1024)
                        println("           estimated footprint = $(round(gib(est); digits=2)) GiB, free now = $(round(gib(free); digits=2)) GiB")
                        continue
                    end

                    ins = make_inputs(FT, N, batch)

                    current_outs = make_outputs(ins)
                    reduced_outs = make_outputs(ins)
                    bundled_outs = make_outputs(ins)

                    current_scratch = make_current_scratch(ins)
                    reduced_scratch = make_reduced_scratch(ins)
                    bundled_scratch = make_bundled_scratch(ins)

                    t_current = time_interaction(cublas_current_interaction11!, current_outs, ins, current_scratch)
                    t_reduced = time_interaction(cublas_reduced_interaction11!, reduced_outs, ins, reduced_scratch)
                    t_bundled = time_interaction(ka_interaction11_bundled!, bundled_outs, ins, bundled_scratch)

                    cublas_reduced_interaction11!(reduced_outs, ins, reduced_scratch)
                    ka_interaction11_bundled!(bundled_outs, ins, bundled_scratch)
                    err = max_output_diff(reduced_outs, bundled_outs)

                    red_spd = t_current / t_reduced
                    bnd_spd = t_current / t_bundled
                    smem_kb = bundled_smem_bytes(N, FT) / 1024

                    @printf("%-8s | %4d | %6d | %12.4f | %12.4f | %12.4f | %7.2fx | %7.2fx | %10.2e | %10.1f\n",
                            string(FT), N, batch,
                            t_current * 1000, t_reduced * 1000, t_bundled * 1000,
                            red_spd, bnd_spd, err, smem_kb)

                    free_namedtuple!(current_scratch)
                    free_namedtuple!(reduced_scratch)
                    free_namedtuple!(bundled_scratch)

                    free_namedtuple!(current_outs)
                    free_namedtuple!(reduced_outs)
                    free_namedtuple!(bundled_outs)
                    free_namedtuple!(ins)
                end
            end
        end

        println("\n", "=" ^ 164)
        println("Reduced ms: one-inverse algebraic rewrite using the identity-based formulation.")
        println("Bundled ms: same reduced math executed inside one KA kernel per batch element.")
        println("Max diff compares bundled output against the reduced-algebra cuBLAS reference.")
        println("=" ^ 164)
    finally
        CUDA.math_mode!(starting_math_mode)
    end
end

run_benchmark()
