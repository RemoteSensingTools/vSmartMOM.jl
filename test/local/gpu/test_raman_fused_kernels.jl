#=
================================================================================
 test_raman_fused_kernels.jl — fused vs unfused (batched_mul) regression suite
================================================================================
Guards the GPU-fused RRS kernels in src/CoreRT/CoreKernel/fused_raman_kernels.jl
against EVER diverging from the original batched_mul math. Two layers:

  1. PER-KERNEL (the 8 fused launchers): each `apply_fused_*!` is run on the GPU
     and compared against an INDEPENDENT batched reference (plain `*` per (n₁,Δn)
     loop) computed here from the kernel's documented formula. Random operands,
     so this isolates each kernel's arithmetic from the rest of the solver. If a
     kernel's formula drifts, exactly that testset fails.

  2. INTEGRATION (the two real entry points): `doubling_inelastic!` and
     `interaction!` are run GPU (→ fused path, since NquadN fits) and CPU
     (→ batched_mul fallback, KA can't run the @synchronize kernels there) on
     bit-identical inputs and compared. This also covers the call-site WIRING
     (argument order) that the per-kernel tests can't.

GPU-only (the fused path is GPU-only by construction); skipped without CUDA.
f64 must agree to ~machine precision; f32 to a loose tol (different accumulation
order — see the fused_raman_kernels.jl header "ACCURACY / GOLDENS").
================================================================================
=#
using Test
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.Architectures
using LinearAlgebra
using Random

const _CUDA_OK = try
    using CUDA
    CUDA.functional()
catch
    false
end

relerr(a, b) = maximum(abs.(Array(a) .- Array(b))) / max(eps(Float64), maximum(abs.(Float64.(Array(b)))))

@testset "Raman fused kernels: fused ≡ batched_mul" begin
    if !_CUDA_OK
        @test_skip "no functional CUDA — fused kernels are GPU-only"
    else
        CUDA.allowscalar(false)
        N, nSpec, nRaman = 12, 40, 8          # N=12 fits every fused predicate (f64 & f32)
        shifts = collect(-3:4)                 # length nRaman, mixed sign; some n₀ out of range
        inr(n0) = 1 <= n0 <= nSpec
        gpu = CuArray

        # tolerances per float type
        tols = ((Float64, 1e-10), (Float32, 5e-3))

        # ---------- per-kernel reference helpers (batched, plain `*`) ----------
        # slice accessors keep the references readable: M4(X,n1,Δn)=X[:,:,n1,Δn] etc.
        M4(X, n1, Δn) = @view X[:, :, n1, Δn]
        M3(X, n)      = @view X[:, :, n]

        @testset "$nm ($FT)" for (nm, FT, tol) in
                [(nm, FT, tol) for nm in
                    ["doubling_nloop2","doubling_jsrc","inter_Rmp","inter_Tmm",
                     "inter_Tpp","inter_Rpm2","inter_Jm","inter_Jp"]
                 for (FT, tol) in tols]

            Random.seed!(20260619)
            m4() = FT.(0.05 .* rand(N, N, nSpec, nRaman))
            s4() = FT.(0.05 .* rand(N, 1, nSpec, nRaman))
            m3() = FT.(0.05 .* rand(N, N, nSpec))
            s3() = FT.(0.05 .* rand(N, 1, nSpec))

            if nm == "doubling_nloop2"
                IET0, IER0 = m4(), m4(); TTGP, GP, R, T = m3(), m3(), m3(), m3()
                refIET, refIER = copy(IET0), copy(IER0)
                for Δn in 1:nRaman, n1 in 1:nSpec
                    n0 = n1 + shifts[Δn]; inr(n0) || continue
                    iet, ier = M4(IET0,n1,Δn), M4(IER0,n1,Δn)
                    Q = ier*M3(R,n0) + M3(R,n1)*ier
                    GPt0 = M3(GP,n0)*M3(T,n0); GPr0 = M3(GP,n0)*M3(R,n0)
                    refIET[:,:,n1,Δn] .= M3(TTGP,n1)*(iet + Q*GPt0) + iet*GPt0
                    refIER[:,:,n1,Δn] .= ier + (iet*GPr0)*M3(T,n0) +
                                         M3(TTGP,n1)*(M3(R,n1)*iet + (ier + Q*GPr0)*M3(T,n0))
                end
                gIET, gIER = gpu(copy(IET0)), gpu(copy(IER0))
                CoreRT.apply_fused_doubling_nloop2!(gIET, gIER, gpu(TTGP), gpu(GP), gpu(R), gpu(T), gpu(shifts), GPU())
                CUDA.synchronize()
                @test relerr(gIET, refIET) < tol
                @test relerr(gIER, refIER) < tol

            elseif nm == "doubling_jsrc"
                ieJ0p0, ieJ0m0 = s4(), s4()
                tt, r = m3(), m3(); ier, iet = m4(), m4()
                tmp1, tmp2, J1m, J0p3 = s3(), s3(), s3(), s3(); expk = FT.(0.9 .+ 0.05 .* rand(nSpec))
                refP, refM = copy(ieJ0p0), copy(ieJ0m0)
                for Δn in 1:nRaman, n1 in 1:nSpec
                    n0 = n1 + shifts[Δn]; inr(n0) || continue
                    J0p, J0m = M4(ieJ0p0,n1,Δn), M4(ieJ0m0,n1,Δn)
                    ek = expk[n0]; J1p = J0p*ek; J1m_ = J0m*ek
                    ie = M4(ier,n1,Δn); it = M4(iet,n1,Δn)
                    Q = M3(r,n1)*ie + ie*M3(r,n0)
                    X = J0p + M3(r,n1)*J1m_ + ie*M3(J1m,n0) + Q*M3(tmp1,n0)
                    Y = J1m_ + ie*M3(J0p3,n0) + M3(r,n1)*J0p + Q*M3(tmp2,n0)
                    refP[:,:,n1,Δn] .= J1p + M3(tt,n1)*X + it*M3(tmp1,n0)
                    refM[:,:,n1,Δn] .= J0m + M3(tt,n1)*Y + it*M3(tmp2,n0)
                end
                gP, gM = gpu(copy(ieJ0p0)), gpu(copy(ieJ0m0))
                CoreRT.apply_fused_doubling_jsrc!(gP, gM, gpu(tt), gpu(r), gpu(ier), gpu(iet),
                    gpu(tmp1), gpu(tmp2), gpu(J1m), gpu(J0p3), gpu(expk), gpu(shifts), GPU())
                CUDA.synchronize()
                @test relerr(gP, refP) < tol
                @test relerr(gM, refM) < tol

            elseif nm == "inter_Rmp"
                T01,TINV,r,Rpm,Tpp = m3(),m3(),m3(),m3(),m3()
                ier,ieRpm,ieTpp,ieTmm,ieRmp = m4(),m4(),m4(),m4(),m4()
                O = zeros(FT,N,N,nSpec,nRaman)
                for Δn in 1:nRaman, n1 in 1:nSpec
                    n0=n1+shifts[Δn]; inr(n0)||continue
                    ie,ieR,ieT = M4(ier,n1,Δn),M4(ieRpm,n1,Δn),M4(ieTpp,n1,Δn)
                    W = M3(T01,n1)*(ie*M3(Rpm,n0)+M3(r,n1)*ieR) + M4(ieTmm,n1,Δn); Wt = W*M3(TINV,n0)
                    O[:,:,n1,Δn] .= M4(ieRmp,n1,Δn) + M3(T01,n1)*(ie*M3(Tpp,n0)+M3(r,n1)*ieT) + (Wt*M3(r,n0))*M3(Tpp,n0)
                end
                gO = gpu(zeros(FT,N,N,nSpec,nRaman))
                CoreRT.apply_fused_interaction_Rmp!(gO, gpu(T01),gpu(TINV),gpu(r),gpu(Rpm),gpu(Tpp),
                    gpu(ier),gpu(ieRpm),gpu(ieTpp),gpu(ieTmm),gpu(ieRmp), gpu(shifts), GPU())
                CUDA.synchronize(); @test relerr(gO, O) < tol

            elseif nm == "inter_Tmm"
                T01,TINV,r,Rpm,tmm = m3(),m3(),m3(),m3(),m3()
                ier,ieRpm,ieTmm,ietmm = m4(),m4(),m4(),m4()
                O = zeros(FT,N,N,nSpec,nRaman)
                for Δn in 1:nRaman, n1 in 1:nSpec
                    n0=n1+shifts[Δn]; inr(n0)||continue
                    ie,ieR = M4(ier,n1,Δn),M4(ieRpm,n1,Δn)
                    W = M3(T01,n1)*(ie*M3(Rpm,n0)+M3(r,n1)*ieR) + M4(ieTmm,n1,Δn); Wt = W*M3(TINV,n0)
                    O[:,:,n1,Δn] .= M3(T01,n1)*M4(ietmm,n1,Δn) + Wt*M3(tmm,n0)
                end
                gO = gpu(zeros(FT,N,N,nSpec,nRaman))
                CoreRT.apply_fused_interaction_Tmm!(gO, gpu(T01),gpu(TINV),gpu(r),gpu(Rpm),gpu(tmm),
                    gpu(ier),gpu(ieRpm),gpu(ieTmm),gpu(ietmm), gpu(shifts), GPU())
                CUDA.synchronize(); @test relerr(gO, O) < tol

            elseif nm == "inter_Tpp"
                T21,TINV,r,Rpm,Tpp = m3(),m3(),m3(),m3(),m3()
                ier,ieRpm,ieTpp,iet = m4(),m4(),m4(),m4()
                O = zeros(FT,N,N,nSpec,nRaman)
                for Δn in 1:nRaman, n1 in 1:nSpec
                    n0=n1+shifts[Δn]; inr(n0)||continue
                    ie,ieR = M4(ier,n1,Δn),M4(ieRpm,n1,Δn)
                    W = M3(T21,n1)*(ieR*M3(r,n0)+M3(Rpm,n1)*ie) + M4(iet,n1,Δn); Wt = W*M3(TINV,n0)
                    O[:,:,n1,Δn] .= M3(T21,n1)*M4(ieTpp,n1,Δn) + Wt*M3(Tpp,n0)
                end
                gO = gpu(zeros(FT,N,N,nSpec,nRaman))
                CoreRT.apply_fused_interaction_Tpp!(gO, gpu(T21),gpu(TINV),gpu(r),gpu(Rpm),gpu(Tpp),
                    gpu(ier),gpu(ieRpm),gpu(ieTpp),gpu(iet), gpu(shifts), GPU())
                CUDA.synchronize(); @test relerr(gO, O) < tol

            elseif nm == "inter_Rpm2"
                T21,TINV,r,Rpm,tmm = m3(),m3(),m3(),m3(),m3()
                ier,ieRpm,ietmm,iet,ierpm = m4(),m4(),m4(),m4(),m4()
                O = zeros(FT,N,N,nSpec,nRaman)
                for Δn in 1:nRaman, n1 in 1:nSpec
                    n0=n1+shifts[Δn]; inr(n0)||continue
                    ie,ieR = M4(ier,n1,Δn),M4(ieRpm,n1,Δn)
                    W = M3(T21,n1)*(ieR*M3(r,n0)+M3(Rpm,n1)*ie) + M4(iet,n1,Δn); Wt = W*M3(TINV,n0)
                    A2 = ieR*M3(tmm,n0) + M3(Rpm,n1)*M4(ietmm,n1,Δn)
                    O[:,:,n1,Δn] .= M4(ierpm,n1,Δn) + M3(T21,n1)*A2 + (Wt*M3(Rpm,n0))*M3(tmm,n0)
                end
                gO = gpu(zeros(FT,N,N,nSpec,nRaman))
                CoreRT.apply_fused_interaction_Rpm2!(gO, gpu(T21),gpu(TINV),gpu(r),gpu(Rpm),gpu(tmm),
                    gpu(ier),gpu(ieRpm),gpu(ietmm),gpu(iet),gpu(ierpm), gpu(shifts), GPU())
                CUDA.synchronize(); @test relerr(gO, O) < tol

            elseif nm == "inter_Jm"
                T01,TINV,r,Rpm = m3(),m3(),m3(),m3()
                ier,ieRpm,ieTmm = m4(),m4(),m4()
                ieJ0m,ieJ0p,addJ0m = s4(),s4(),s4(); J0p,addj0m = s3(),s3()
                O = zeros(FT,N,1,nSpec,nRaman)
                for Δn in 1:nRaman, n1 in 1:nSpec
                    n0=n1+shifts[Δn]; inr(n0)||continue
                    ie,ieR = M4(ier,n1,Δn),M4(ieRpm,n1,Δn)
                    W1 = M3(T01,n1)*(ie*M3(Rpm,n0)+M3(r,n1)*ieR) + M4(ieTmm,n1,Δn); W1t = W1*M3(TINV,n0)
                    v0 = M3(addj0m,n0) + M3(r,n0)*M3(J0p,n0)
                    inner = ie*M3(J0p,n0) + M3(r,n1)*M4(ieJ0p,n1,Δn) + M4(addJ0m,n1,Δn)
                    O[:,:,n1,Δn] .= M4(ieJ0m,n1,Δn) + M3(T01,n1)*inner + W1t*v0
                end
                gO = gpu(zeros(FT,N,1,nSpec,nRaman))
                CoreRT.apply_fused_interaction_Jm!(gO, gpu(T01),gpu(TINV),gpu(ier),gpu(r),gpu(Rpm),
                    gpu(ieRpm),gpu(ieTmm),gpu(ieJ0m),gpu(J0p),gpu(ieJ0p),gpu(addJ0m),gpu(addj0m), gpu(shifts), GPU())
                CUDA.synchronize(); @test relerr(gO, O) < tol

            else  # inter_Jp
                T21,TINV,r,Rpm = m3(),m3(),m3(),m3()
                ier,ieRpm,iet = m4(),m4(),m4()
                ieJ0p,addieJ0p,addieJ0m = s4(),s4(),s4(); J0p,addj0m = s3(),s3()
                O = zeros(FT,N,1,nSpec,nRaman)
                for Δn in 1:nRaman, n1 in 1:nSpec
                    n0=n1+shifts[Δn]; inr(n0)||continue
                    ie,ieR = M4(ier,n1,Δn),M4(ieRpm,n1,Δn)
                    W2 = M3(T21,n1)*(ieR*M3(r,n0)+M3(Rpm,n1)*ie) + M4(iet,n1,Δn); W2t = W2*M3(TINV,n0)
                    v0 = M3(J0p,n0) + M3(Rpm,n0)*M3(addj0m,n0)
                    inner = M4(ieJ0p,n1,Δn) + ieR*M3(addj0m,n0) + M3(Rpm,n1)*M4(addieJ0m,n1,Δn)
                    O[:,:,n1,Δn] .= M4(addieJ0p,n1,Δn) + M3(T21,n1)*inner + W2t*v0
                end
                gO = gpu(zeros(FT,N,1,nSpec,nRaman))
                CoreRT.apply_fused_interaction_Jp!(gO, gpu(T21),gpu(TINV),gpu(ier),gpu(r),gpu(Rpm),
                    gpu(ieRpm),gpu(iet),gpu(ieJ0p),gpu(J0p),gpu(addieJ0p),gpu(addieJ0m),gpu(addj0m), gpu(shifts), GPU())
                CUDA.synchronize(); @test relerr(gO, O) < tol
            end
        end

        # ---------- INTEGRATION: full kernels, fused (GPU) ≡ fallback (CPU) ----------
        # Covers call-site wiring + the whole doubling/interaction chains. f64 only
        # (cross-device machine-precision agreement; f32 cross-device is looser).
        @testset "integration: $nm fused(GPU) ≡ fallback(CPU)" for nm in ["doubling_inelastic!","interaction!"]
            FT = Float64; Nq = 15; nS = 30; nR = 16
            half = nR ÷ 2
            mkRS() = InelasticScattering.RRS(
                n2 = InelasticScattering.getRamanAtmoConstants(FT(1e7/760), FT(250))[1],
                o2 = InelasticScattering.getRamanAtmoConstants(FT(1e7/760), FT(250))[2],
                greek_raman = InelasticScattering.GreekCoefs(fill(FT(1),1),fill(FT(1),1),fill(FT(1),1),fill(FT(1),1),fill(FT(1),1),fill(FT(1),1)),
                fscattRayl = [FT(1)], ϖ_Cabannes = [FT(1)], ϖ_λ₁λ₀ = zeros(FT,nR),
                i_λ₁λ₀ = collect(-half:(nR-half-1)), Z⁻⁺_λ₁λ₀ = zeros(FT,1,1), Z⁺⁺_λ₁λ₀ = zeros(FT,1,1),
                i_ref = nS÷2, n_Raman = nR, F₀ = zeros(FT,3,nS), SIF₀ = zeros(FT,3,nS))
            pol = CoreRT.Stokes_IQU()
            RS = mkRS()
            fillrand!(L, seed) = (Random.seed!(seed); for f in fieldnames(typeof(L))
                v = getfield(L,f); v isa AbstractArray && eltype(v)<:AbstractFloat && copyto!(v, 0.05 .* rand(Float64, size(v))); end; L)

            if nm == "doubling_inelastic!"
                aC = fillrand!(CoreRT.make_added_layer(RS, FT, Array, (Nq,Nq), nS), 7)
                aG = fillrand!(CoreRT.make_added_layer(RS, FT, CuArray, (Nq,Nq), nS), 7)
                expk = fill(FT(0.9), nS)
                IsC = Diagonal(Array(Diagonal{FT}(ones(Nq))));  IsG = Diagonal(CuArray(Diagonal{FT}(ones(Nq))))
                CoreRT.doubling_inelastic!(RS, pol, true, copy(expk), 4, aC, IsC, CPU())
                CoreRT.doubling_inelastic!(RS, pol, true, CuArray(copy(expk)), 4, aG, IsG, GPU()); CUDA.synchronize()
                for f in (:r⁻⁺,:t⁺⁺,:ier⁻⁺,:iet⁺⁺,:ieJ₀⁺,:ieJ₀⁻,:j₀⁺,:j₀⁻)
                    @test relerr(getfield(aG,f), getfield(aC,f)) < 1e-11
                end
            else
                SI = CoreRT.ScatteringInterface_11()
                cC = fillrand!(CoreRT.make_composite_layer(RS, FT, Array, (Nq,Nq), nS), 3)
                aC = fillrand!(CoreRT.make_added_layer(RS, FT, Array, (Nq,Nq), nS), 5)
                cG = fillrand!(CoreRT.make_composite_layer(RS, FT, CuArray, (Nq,Nq), nS), 3)
                aG = fillrand!(CoreRT.make_added_layer(RS, FT, CuArray, (Nq,Nq), nS), 5)
                IsC = Diagonal(Array(Diagonal{FT}(ones(Nq))));  IsG = Diagonal(CuArray(Diagonal{FT}(ones(Nq))))
                CoreRT.interaction!(RS, SI, true, cC, aC, IsC; workspace=nothing)
                CoreRT.interaction!(RS, SI, true, cG, aG, IsG; workspace=nothing); CUDA.synchronize()
                for f in (:R⁻⁺,:T⁺⁺,:R⁺⁻,:T⁻⁻,:ieR⁻⁺,:ieT⁺⁺,:ieR⁺⁻,:ieT⁻⁻,:ieJ₀⁺,:ieJ₀⁻,:J₀⁺,:J₀⁻)
                    @test relerr(getfield(cG,f), getfield(cC,f)) < 1e-11
                end
            end
        end
    end
end
