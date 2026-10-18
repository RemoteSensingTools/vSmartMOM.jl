#=
Polarized Cox–Munk geometry, conventions and Fourier storage (release audit
finding B3). Every reference here is independent of the Cox–Munk code:

1. The Stokes basis of `JonesReference` is first shown to be vSmartMOM's own,
   by reproducing the atmospheric Rayleigh phase matrix that the solver uses.
2. The facet Mueller matrix is compared with that Jones calculation (random
   geometries, complex indices, both precisions) and with stored VLIDORT
   2.8.3 kernel values (test/reference/coxmunk_vlidort_kernels.jl).
3. The surface Fourier moments are required to follow the storage
   convention of the atmospheric moments `compute_Z_moments`, because the
   adding method combines both in the same operator slots.
4. The reflection operator maps incident column j to reflected row i.
=#
using Test, LinearAlgebra, Random, StaticArrays
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering

const CM = vSmartMOM.CoreRT
include("coxmunk_jones_reference.jl")
include(joinpath("reference", "coxmunk_vlidort_kernels.jl"))
using .JonesReference

per_axis_variance(U) = (0.003 + 0.00512U) / 2
surface(U::FT, n; shadow = false) where FT =
    CoxMunkSurface{FT}(wind_speed = U, n_water = Complex{FT}(n),
                       include_whitecaps = false, shadowing = shadow)

"Real-space matrix from vSmartMOM Fourier moments stored in the Z⁻⁺ convention."
function synthesize(moment, m_max, Δϕ)
    M = zeros(4, 4)
    for m in 0:m_max, i in 1:4, j in 1:4
        M[i, j] += moment(m)[i, j] * CM._azimuthal_kernel(i, j, m, Δϕ)
    end
    return M
end

@testset "Stokes convention anchor: Rayleigh Z⁻⁺ = dipole Jones matrix" begin
    greek = Scattering.get_greek_rayleigh(0.0)
    for (μ₀, μv, Δϕ) in ((cosd(30), cosd(45), deg2rad(20)), (cosd(60), cosd(10), deg2rad(135)),
                         (cosd(15), cosd(70), deg2rad(250)))
        Z(m) = Scattering.compute_Z_moments(Stokes_IQUV(), [μv, μ₀], greek, m)[2][1:4, 5:8]
        @test synthesize(Z, 2, Δϕ) ≈ JonesReference.rayleigh(μ₀, μv, Δϕ) atol = 1e-12
    end
end

@testset "Stokes convention anchor without the Fourier kernel" begin
    # Unpolarized incidence: the single-scattering column of the code's
    # real-space Rayleigh matrix fixes Q and U of the Jones basis directly.
    greek = Scattering.get_greek_rayleigh(0.0)
    for (μ₀, μv, Δϕ) in ((cosd(30), cosd(45), deg2rad(20)), (cosd(70), cosd(25), deg2rad(300)))
        col = collect(Scattering.phase_matrix_first_column(greek, μ₀, μv, Δϕ, Val(4)))
        @test col ≈ JonesReference.rayleigh(μ₀, μv, Δϕ)[:, 1] atol = 1e-12
    end
end

@testset "Full-matrix reciprocity" begin
    # Reference-free physics: M(μᵢ, μᵣ, Δϕ) = Δ Mᵀ(μᵣ, μᵢ, −Δϕ) Δ, Δ = diag(1, 1, −1, 1)
    # (Hovenier 1969). An argument-order or reference-plane error breaks it.
    rng = MersenneTwister(7)
    Δ = Diagonal([1.0, 1.0, -1.0, 1.0])
    for _ in 1:200
        μᵢ, μᵣ, Δϕ = 0.05 + 0.95rand(rng), 0.05 + 0.95rand(rng), 2π * rand(rng)
        n = complex(1.2 + 0.4rand(rng), 0.2rand(rng))
        surf = surface(12rand(rng), n; shadow = true)
        M = CM.coxmunk_brdf_mueller(surf, 4, μᵢ, μᵣ, Δϕ; n_water = n)
        Mt = CM.coxmunk_brdf_mueller(surf, 4, μᵣ, μᵢ, -Δϕ; n_water = n)
        @test M ≈ Δ * transpose(Mt) * Δ rtol = 1e-12 atol = 1e-14 * abs(M[1, 1])
    end
end

@testset "Float32 at grazing, nearly parallel rays" begin
    # cos θ_local must not suffer cancellation (review finding on the B3 change).
    for (μᵢ, μᵣ, Δϕ) in ((1e-4, 2e-4, 0.0), (5e-3, 3e-3, 0.0), (0.02, 0.03, 0.005))
        M64 = CM.coxmunk_brdf_mueller(surface(5.0, complex(1.33, 0.0)), 3, μᵢ, μᵣ, Δϕ;
                                      n_water = complex(1.33, 0.0))
        M32 = CM.coxmunk_brdf_mueller(surface(5f0, ComplexF32(1.33)), 3, Float32(μᵢ),
                                      Float32(μᵣ), Float32(Δϕ); n_water = ComplexF32(1.33))
        @test maximum(abs.(M32 .- M64)) / M64[1, 1] < 1e-4
    end
end

@testset "Facet Mueller matrix vs Jones reference ($FT)" for FT in (Float64, Float32)
    rng = MersenneTwister(20261017)
    tol = FT == Float64 ? 1e-11 : 2e-4
    for _ in 1:400
        μᵢ, μᵣ, Δϕ = 0.05 + 0.95rand(rng), 0.05 + 0.95rand(rng), 2π * rand(rng)
        n = complex(1.2 + 0.4rand(rng), rand(rng) < 0.5 ? 0.0 : 0.4rand(rng))
        U = 12rand(rng)
        ref = JonesReference.coxmunk(μᵢ, μᵣ, Δϕ, n, per_axis_variance(U))
        ref[1, 1] > 1e-3 || continue          # skip far tails, relative error is ill-defined
        M = CM.coxmunk_brdf_mueller(surface(FT(U), Complex{FT}(n)), 4,
                                    FT(μᵢ), FT(μᵣ), FT(Δϕ); n_water = Complex{FT}(n))
        @test maximum(abs.(M .- ref)) / ref[1, 1] < tol
        for ns in 1:3                           # truncation of the full product
            Mn = CM.coxmunk_brdf_mueller(surface(FT(U), Complex{FT}(n)), ns,
                                         FT(μᵢ), FT(μᵣ), FT(Δϕ); n_water = Complex{FT}(n))
            @test maximum(abs.(Mn .- M[1:ns, 1:ns])) <= 4eps(FT) * abs(M[1, 1])
        end
    end
end

@testset "Facet Mueller matrix vs VLIDORT 2.8.3 kernels" begin
    for c in COXMUNK_VLIDORT_KERNELS
        M = CM.coxmunk_brdf_mueller(surface(c.U, c.n; shadow = c.shadow), 4,
                                    cosd(c.sza), cosd(c.vza), deg2rad(c.raz); n_water = c.n)
        if iszero(c.M[1, 1])
            # VLIDORT sets the slope PDF to zero once its exponent exceeds 88.
            @test abs(M[1, 1]) < exp(-88)
        else
            # VLIDORT clamps μ = 1 to 1 − 1e-12, and its shadowing uses an erfc
            # approximation (DERFC_E, fractional error < 1.2e-7).
            @test M ≈ c.M rtol = (c.shadow ? 1e-7 : 1e-9) atol = 1e-12 * abs(c.M[1, 1])
        end
    end
end

@testset "Fourier storage convention matches atmospheric Z⁻⁺" begin
    # Decompose the SAME real-space matrix (vSmartMOM's Rayleigh matrix) with
    # the surface's azimuthal kernel and compare with compute_Z_moments.
    greek = Scattering.get_greek_rayleigh(0.0)
    μ₀, μv = cosd(30), cosd(45)
    ϕ, w = CM.CanopyOptics.gauleg(100, 0.0, Float64(π))
    for m in 0:2
        Zm = Scattering.compute_Z_moments(Stokes_IQUV(), [μv, μ₀], greek, m)[2][1:4, 5:8]
        R = zeros(4, 4)
        for (ϕk, wk) in zip(ϕ, w), i in 1:4, j in 1:4
            R[i, j] += wk * JonesReference.rayleigh(μ₀, μv, ϕk)[i, j] *
                       CM._azimuthal_kernel(i, j, m, ϕk)
        end
        @test (m == 0 ? 1 : 2) / π .* R ≈ Zm atol = 1e-12
    end
end

@testset "Reflection operator: incident column, reflected row, full synthesis" begin
    n = complex(1.33, 0.0)
    surf = surface(5.0, n; shadow = true)
    μ = [0.35, 0.6, 0.85]
    m_max = 60
    moments = [CM.reflectance(surf, Stokes_IQUV(), μ, m; n_water = n) for m in 0:m_max]
    for (i, j, Δϕ) in ((1, 3, 0.4), (3, 1, 2.0), (2, 3, 4.5))
        block(m) = moments[m + 1][(i - 1) * 4 .+ (1:4), (j - 1) * 4 .+ (1:4)] ./ π
        expected = CM.coxmunk_brdf_mueller(surf, 4, μ[j], μ[i], Δϕ; n_water = n)
        @test synthesize(block, m_max, Δϕ) ≈ expected rtol = 1e-6 atol = 1e-8
    end
end
