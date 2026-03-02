#=
Unit tests for the Cox-Munk ocean surface BRDF implementation.

Tests cover:
  - Fresnel coefficients and Mueller matrix (analytical checks)
  - Stokes rotation matrix properties
  - Water refractive index lookup table
  - Cox-Munk helper functions (slope variance, whitecap fraction, shadow masking)
  - Core BRDF Mueller matrix (reciprocity, energy conservation, limiting cases)
  - Polarized Fourier decomposition
  - Stokes convention verification
  - Jacobian finite-difference consistency
=#

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering
using Test
using LinearAlgebra
using StaticArrays

# Access internals via the CoreRT module
const CM = vSmartMOM.CoreRT

# ═══════════════════════════════════════════════════════════════════════
# 1. Fresnel coefficients
# ═══════════════════════════════════════════════════════════════════════

@testset "Fresnel coefficients" begin
    n_water = complex(1.33, 0.0)

    @testset "Normal incidence" begin
        r_s, r_p = CM.fresnel_coefficients(n_water, 1.0)
        R_normal = abs2(r_s)
        @test R_normal ≈ ((1.33 - 1) / (1.33 + 1))^2 atol = 1e-10
        # At normal incidence |r_s| = |r_p| (r_s = -r_p by sign convention)
        @test abs2(r_s) ≈ abs2(r_p) atol = 1e-12
    end

    @testset "Brewster angle" begin
        θ_B = atan(1.33)
        r_s_B, r_p_B = CM.fresnel_coefficients(n_water, cos(θ_B))
        @test abs(r_p_B) < 1e-10
    end

    @testset "Grazing incidence" begin
        r_s_g, r_p_g = CM.fresnel_coefficients(n_water, 0.001)
        @test abs2(r_s_g) > 0.99
        @test abs2(r_p_g) > 0.99
    end

    @testset "Complex refractive index" begin
        n_complex = complex(1.33, 0.01)
        r_s, r_p = CM.fresnel_coefficients(n_complex, 0.5)
        @test isfinite(abs2(r_s))
        @test isfinite(abs2(r_p))
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 2. Fresnel Mueller matrix
# ═══════════════════════════════════════════════════════════════════════

@testset "Fresnel Mueller matrix" begin
    n_water = complex(1.33, 0.0)
    r_s, r_p = CM.fresnel_coefficients(n_water, 0.7)
    M = CM.fresnel_mueller(r_s, r_p, 4)

    @testset "Block structure" begin
        # Off-diagonal blocks should be zero
        @test M[1, 3] ≈ 0.0 atol = 1e-15
        @test M[1, 4] ≈ 0.0 atol = 1e-15
        @test M[2, 3] ≈ 0.0 atol = 1e-15
        @test M[2, 4] ≈ 0.0 atol = 1e-15
        @test M[3, 1] ≈ 0.0 atol = 1e-15
        @test M[3, 2] ≈ 0.0 atol = 1e-15
        @test M[4, 1] ≈ 0.0 atol = 1e-15
        @test M[4, 2] ≈ 0.0 atol = 1e-15
    end

    @testset "Symmetry" begin
        @test M[1, 2] ≈ M[2, 1]
        @test M[3, 3] ≈ M[4, 4]
        @test M[3, 4] ≈ -M[4, 3]
    end

    @testset "Energy conservation" begin
        @test M[1, 1] >= abs(M[1, 2])
        @test M[1, 1] <= 1.0
    end

    @testset "Consistent sizes" begin
        M1 = CM.fresnel_mueller(r_s, r_p, 1)
        @test size(M1) == (1, 1)
        M3 = CM.fresnel_mueller(r_s, r_p, 3)
        @test size(M3) == (3, 3)
        M4 = CM.fresnel_mueller(r_s, r_p, 4)
        @test size(M4) == (4, 4)
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 3. Stokes rotation matrix
# ═══════════════════════════════════════════════════════════════════════

@testset "Stokes rotation matrix" begin
    @testset "Identity at zero angle" begin
        L0 = CM.stokes_rotation_matrix(0.0, 4)
        @test L0 ≈ SMatrix{4,4}(I) atol = 1e-15
    end

    @testset "Composition L(a)·L(b) = L(a+b)" begin
        a, b = 0.3, 0.7
        La = CM.stokes_rotation_matrix(a, 4)
        Lb = CM.stokes_rotation_matrix(b, 4)
        Lab = CM.stokes_rotation_matrix(a + b, 4)
        @test La * Lb ≈ Lab atol = 1e-12
    end

    @testset "Inverse: L(-a) = L(a)⁻¹" begin
        a = 1.2
        La = CM.stokes_rotation_matrix(a, 4)
        La_inv = CM.stokes_rotation_matrix(-a, 4)
        @test La * La_inv ≈ SMatrix{4,4}(I) atol = 1e-12
    end

    @testset "Orthogonality" begin
        L = CM.stokes_rotation_matrix(0.8, 4)
        # L^T L should be identity (orthogonal for the QU subblock)
        @test L' * L ≈ SMatrix{4,4}(I) atol = 1e-12
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 4. Water refractive index lookup
# ═══════════════════════════════════════════════════════════════════════

@testset "Water refractive index" begin
    @testset "Visible range" begin
        n550 = CM.water_refractive_index(550.0)
        @test real(n550) ≈ 1.333 atol = 0.005
        @test imag(n550) < 1e-6  # nearly transparent in visible
    end

    @testset "Near-IR absorption" begin
        n1500 = CM.water_refractive_index(1500.0)
        @test imag(n1500) > imag(CM.water_refractive_index(550.0))
    end

    @testset "UV edge" begin
        n200 = CM.water_refractive_index(200.0)
        @test real(n200) > 1.35  # higher dispersion in UV
    end

    @testset "Boundary clamping" begin
        # Should not error at extremes
        n_lo = CM.water_refractive_index(100.0)
        n_hi = CM.water_refractive_index(10000.0)
        @test isfinite(real(n_lo))
        @test isfinite(real(n_hi))
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 5. Cox-Munk helper functions
# ═══════════════════════════════════════════════════════════════════════

@testset "Cox-Munk helpers" begin
    @testset "Slope variance" begin
        @test CM.wind_to_sigma2(0.0) ≈ 0.003 atol = 1e-10
        @test CM.wind_to_sigma2(10.0) ≈ 0.003 + 0.00512 * 10 atol = 1e-10
    end

    @testset "Whitecap fraction" begin
        @test CM.whitecap_fraction(0.0) ≈ 0.0 atol = 1e-15
        @test CM.whitecap_fraction(10.0) ≈ 2.95e-6 * 10.0^3.52 rtol = 1e-6
        # At high wind: f_wc should be a few percent
        @test 0.01 < CM.whitecap_fraction(15.0) < 0.10
    end

    @testset "Slope PDF normalization" begin
        σ² = CM.wind_to_sigma2(5.0)
        # Numerical integration of the PDF over slopes should ≈ 1
        N = 200
        z_max = 4 * sqrt(σ²)
        dz = 2 * z_max / N
        integral = 0.0
        for i in 1:N
            zx = -z_max + (i - 0.5) * dz
            for j in 1:N
                zy = -z_max + (j - 0.5) * dz
                integral += CM.cox_munk_pdf(zx, zy, σ²) * dz^2
            end
        end
        @test integral ≈ 1.0 atol = 0.01
    end

    @testset "Shadow masking" begin
        σ² = CM.wind_to_sigma2(5.0)
        # Near-zenith: shadowing ≈ 1
        S_zenith = CM.shadow_factor(0.95, 0.95, σ²)
        @test S_zenith > 0.9
        # Near-grazing: shadowing < 1
        S_grazing = CM.shadow_factor(0.1, 0.1, σ²)
        @test S_grazing < S_zenith
        @test S_grazing > 0.0
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 6. Core BRDF Mueller matrix
# ═══════════════════════════════════════════════════════════════════════

@testset "BRDF Mueller matrix" begin
    surf = CoxMunkSurface(wind_speed=5.0, n_water=complex(1.33, 0.0),
                           include_whitecaps=false)

    @testset "Finiteness" begin
        M = CM.coxmunk_brdf_mueller(surf, 4, 0.6, 0.8, 1.2;
                                      n_water=complex(1.33, 0.0))
        @test all(isfinite.(M))
    end

    @testset "Reciprocity: M[1,1](μᵢ,μᵣ) = M[1,1](μᵣ,μᵢ)" begin
        μ_i, μ_r, dϕ = 0.6, 0.8, 1.2
        M_fwd = CM.coxmunk_brdf_mueller(surf, 4, μ_i, μ_r, dϕ;
                                          n_water=complex(1.33, 0.0))
        M_rev = CM.coxmunk_brdf_mueller(surf, 4, μ_r, μ_i, dϕ;
                                          n_water=complex(1.33, 0.0))
        @test M_fwd[1, 1] ≈ M_rev[1, 1] rtol = 1e-6
    end

    @testset "Specular peak: low wind > high wind at specular geometry" begin
        surf_low  = CoxMunkSurface(wind_speed=1.0, n_water=complex(1.33, 0.0),
                                    include_whitecaps=false)
        surf_high = CoxMunkSurface(wind_speed=20.0, n_water=complex(1.33, 0.0),
                                    include_whitecaps=false)
        # At specular geometry (μᵢ = μᵣ, dϕ = π for backscatter convention)
        # Actually for specular: facet is horizontal, so μᵢ = μᵣ, dϕ = 0
        M_low  = CM.coxmunk_brdf_mueller(surf_low,  1, 0.7, 0.7, 0.0;
                                           n_water=complex(1.33, 0.0))
        M_high = CM.coxmunk_brdf_mueller(surf_high, 1, 0.7, 0.7, 0.0;
                                           n_water=complex(1.33, 0.0))
        @test M_low[1, 1] > 5 * M_high[1, 1]
    end

    @testset "Energy conservation" begin
        # Hemispherical integral of M[1,1] should be ≤ 1
        μ_in = 0.5
        N_μ = 30
        N_ϕ = 60
        integral = 0.0
        dμ = 1.0 / N_μ
        dϕ = 2π / N_ϕ
        for i in 1:N_μ
            μ_out = (i - 0.5) * dμ
            for j in 1:N_ϕ
                ϕ = (j - 0.5) * dϕ
                M = CM.coxmunk_brdf_mueller(surf, 1, μ_out, μ_in, ϕ;
                                              n_water=complex(1.33, 0.0))
                integral += M[1, 1] * μ_out * dμ * dϕ
            end
        end
        @test integral <= 1.0 + 0.02  # small tolerance for quadrature
    end

    @testset "Whitecap contribution" begin
        surf_wc = CoxMunkSurface(wind_speed=12.0, n_water=complex(1.33, 0.0),
                                  include_whitecaps=true, whitecap_albedo=0.22)
        surf_no_wc = CoxMunkSurface(wind_speed=12.0, n_water=complex(1.33, 0.0),
                                     include_whitecaps=false)
        M_wc   = CM.coxmunk_brdf_mueller(surf_wc,    1, 0.3, 0.8, 2.0;
                                           n_water=complex(1.33, 0.0))
        M_no   = CM.coxmunk_brdf_mueller(surf_no_wc, 1, 0.3, 0.8, 2.0;
                                           n_water=complex(1.33, 0.0))
        # With whitecaps, the off-specular regions get a Lambertian floor
        @test M_wc[1, 1] != M_no[1, 1]
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 7. Fourier decomposition
# ═══════════════════════════════════════════════════════════════════════

@testset "Fourier decomposition" begin
    surf = CoxMunkSurface(wind_speed=5.0, n_water=complex(1.33, 0.0))
    μ = Float64[0.3, 0.5, 0.7, 0.9]

    @testset "Scalar (Stokes_I)" begin
        pol_I = Stokes_I{Float64}()
        R0 = CM.reflectance(surf, pol_I, μ, 0; n_water=complex(1.33, 0.0))
        @test size(R0) == (4, 4)
        @test all(isfinite.(R0))
        # m=0 scalar should be all positive
        @test all(R0 .>= -1e-10)
    end

    @testset "Higher moments decrease" begin
        pol_I = Stokes_I{Float64}()
        R0 = CM.reflectance(surf, pol_I, μ, 0; n_water=complex(1.33, 0.0))
        # Use m=10 (high enough that the 2× Fourier factor for m>0 doesn't dominate)
        R10 = CM.reflectance(surf, pol_I, μ, 10; n_water=complex(1.33, 0.0))
        # Normalize out the 2× Fourier convention factor for m>0
        @test maximum(abs.(R10)) / 2 < maximum(abs.(R0))
    end

    @testset "Polarized I-Q coupling (Stokes_IQUV)" begin
        pol_IQUV = Stokes_IQUV{Float64}()
        R0 = CM.reflectance(surf, pol_IQUV, μ, 0; n_water=complex(1.33, 0.0))
        @test size(R0) == (16, 16)
        @test all(isfinite.(R0))

        # I-Q coupling should be nonzero (the whole point of polarized BRDF!)
        has_IQ = false
        for iq in 1:4   # quadrature points
            # Check (I,Q) element at this quadrature point
            idx_I = (iq - 1) * 4 + 1
            idx_Q = (iq - 1) * 4 + 2
            for jq in 1:4
                jdx_I = (jq - 1) * 4 + 1
                if abs(R0[idx_Q, jdx_I]) > 1e-12
                    has_IQ = true
                end
            end
        end
        @test has_IQ
    end

    @testset "Stokes_IQU" begin
        pol_IQU = Stokes_IQU{Float64}()
        R0 = CM.reflectance(surf, pol_IQU, μ, 0; n_water=complex(1.33, 0.0))
        @test size(R0) == (12, 12)
        @test all(isfinite.(R0))
    end
end

# ═══════════════════════════════════════════════════════════════════════
# 8. Stokes convention verification
# ═══════════════════════════════════════════════════════════════════════

@testset "Stokes convention" begin
    # At near-specular geometry, Fresnel reflection produces Q < 0
    # (partially polarized perpendicular to the plane of incidence)
    surf = CoxMunkSurface(wind_speed=3.0, n_water=complex(1.33, 0.0),
                           include_whitecaps=false)

    # Near-specular geometry
    μ = 0.7
    M = CM.coxmunk_brdf_mueller(surf, 4, μ, μ, 0.0; n_water=complex(1.33, 0.0))

    # M[1,1] = (rs² + rp²)/2, M[2,1] = (rs² - rp²)/2
    # Since rs > rp at non-normal incidence: M[2,1] > 0
    # But the sign of Q/I in the reflected Stokes vector = M[2,1]/M[1,1]
    # For water: rs² > rp² → M[2,1] > 0, meaning the reflected Q is positive
    # in the Fresnel sense, which corresponds to polarization perpendicular to
    # the plane of incidence in the convention where Q > 0 = parallel.
    # Actually: for our Mueller matrix, M[1,2] = M[2,1] = (rs²-rp²)/2
    # and rs > rp for incidence angle < Brewster → M[1,2] > 0.
    @test M[1, 1] > 0
    @test M[1, 2] ≈ M[2, 1] atol = 1e-10  # symmetric
    # The magnitude of polarization should be significant
    @test abs(M[1, 2]) / M[1, 1] > 0.01
end

# ═══════════════════════════════════════════════════════════════════════
# 9. Jacobian consistency (finite difference)
# ═══════════════════════════════════════════════════════════════════════

@testset "Jacobian finite-difference" begin
    U = 5.0
    ε = 1e-4
    n_w = complex(1.33, 0.0)
    μ = Float64[0.3, 0.5, 0.7]
    pol = Stokes_I{Float64}()

    # Reflectance at U, U+ε, U-ε
    surf   = CoxMunkSurface(wind_speed=U,     n_water=n_w)
    surf_p = CoxMunkSurface(wind_speed=U + ε, n_water=n_w)
    surf_m = CoxMunkSurface(wind_speed=U - ε, n_water=n_w)

    R   = CM.reflectance(surf,   pol, μ, 0; n_water=n_w)
    R_p = CM.reflectance(surf_p, pol, μ, 0; n_water=n_w)
    R_m = CM.reflectance(surf_m, pol, μ, 0; n_water=n_w)

    dR_fd = (R_p - R_m) / (2ε)

    # The derivative should be nonzero (BRDF depends on wind speed)
    @test maximum(abs.(dR_fd)) > 0.0

    # Check that the sign makes sense: increasing wind spreads the glint,
    # so diagonal elements of R at m=0 should generally change
    @test any(abs.(dR_fd) .> 1e-6)
end

# ═══════════════════════════════════════════════════════════════════════
# 10. Analytical vs finite-difference derivative comparison
# ═══════════════════════════════════════════════════════════════════════

@testset "Analytical vs FD derivatives" begin
    n_w = complex(1.33, 0.0)

    @testset "Pointwise BRDF derivative (coxmunk_brdf_mueller_and_deriv)" begin
        surf = CoxMunkSurface(wind_speed=5.0, n_water=n_w, include_whitecaps=true)
        ε = 1e-5

        for (μᵢ, μᵣ, dϕ) in [(0.6, 0.8, 1.2), (0.3, 0.9, 0.5), (0.7, 0.7, 0.01)]
            M, dM_anal = CM.coxmunk_brdf_mueller_and_deriv(surf, 4, μᵢ, μᵣ, dϕ; n_water=n_w)

            # FD: perturb wind speed
            surf_p = CoxMunkSurface(wind_speed=surf.wind_speed + ε, n_water=n_w,
                                     include_whitecaps=true)
            surf_m = CoxMunkSurface(wind_speed=surf.wind_speed - ε, n_water=n_w,
                                     include_whitecaps=true)
            M_p = CM.coxmunk_brdf_mueller(surf_p, 4, μᵢ, μᵣ, dϕ; n_water=n_w)
            M_m = CM.coxmunk_brdf_mueller(surf_m, 4, μᵢ, μᵣ, dϕ; n_water=n_w)
            dM_fd = (M_p - M_m) / (2ε)

            # Forward value should match
            @test M ≈ CM.coxmunk_brdf_mueller(surf, 4, μᵢ, μᵣ, dϕ; n_water=n_w) atol=1e-14

            # Analytical ≈ FD for each element
            for si in 1:4, sj in 1:4
                if abs(dM_fd[si, sj]) > 1e-15
                    @test dM_anal[si, sj] ≈ dM_fd[si, sj] rtol=1e-3
                else
                    @test abs(dM_anal[si, sj]) < 1e-10
                end
            end
        end
    end

    @testset "Pointwise BRDF derivative — no whitecaps" begin
        surf = CoxMunkSurface(wind_speed=8.0, n_water=n_w, include_whitecaps=false)
        ε = 1e-5
        μᵢ, μᵣ, dϕ = 0.5, 0.6, 0.8

        M, dM_anal = CM.coxmunk_brdf_mueller_and_deriv(surf, 4, μᵢ, μᵣ, dϕ; n_water=n_w)
        surf_p = CoxMunkSurface(wind_speed=surf.wind_speed + ε, n_water=n_w, include_whitecaps=false)
        surf_m = CoxMunkSurface(wind_speed=surf.wind_speed - ε, n_water=n_w, include_whitecaps=false)
        M_p = CM.coxmunk_brdf_mueller(surf_p, 4, μᵢ, μᵣ, dϕ; n_water=n_w)
        M_m = CM.coxmunk_brdf_mueller(surf_m, 4, μᵢ, μᵣ, dϕ; n_water=n_w)
        dM_fd = (M_p - M_m) / (2ε)

        for si in 1:4, sj in 1:4
            if abs(dM_fd[si, sj]) > 1e-15
                @test dM_anal[si, sj] ≈ dM_fd[si, sj] rtol=1e-3
            else
                @test abs(dM_anal[si, sj]) < 1e-10
            end
        end
    end

    @testset "Fourier reflectance derivative (reflectance_and_deriv)" begin
        μ = Float64[0.3, 0.5, 0.7]

        for (U, wc) in [(5.0, true), (10.0, false), (3.0, true)]
            surf = CoxMunkSurface(wind_speed=U, n_water=n_w, include_whitecaps=wc)

            for m in [0, 1, 3]
                for pol in [Stokes_I{Float64}(), Stokes_IQUV{Float64}()]
                    R_anal, dR_anal = CM.reflectance_and_deriv(surf, pol, μ, m; n_water=n_w)
                    _, dR_fd = CM.reflectance_fd_deriv(surf, pol, μ, m; n_water=n_w)

                    # Check agreement element-by-element
                    max_fd = maximum(abs.(dR_fd))
                    if max_fd > 1e-12
                        # Relative tolerance on the maximum element
                        @test maximum(abs.(dR_anal - dR_fd)) / max_fd < 0.01
                    else
                        @test maximum(abs.(dR_anal)) < 1e-10
                    end
                end
            end
        end
    end

    @testset "Derivative helper functions" begin
        # cox_munk_pdf_dσ² via FD
        σ² = 0.03
        zx, zy = 0.1, -0.05
        ε = 1e-7
        dP_anal = CM.cox_munk_pdf_dσ²(zx, zy, σ²)
        dP_fd = (CM.cox_munk_pdf(zx, zy, σ² + ε) - CM.cox_munk_pdf(zx, zy, σ² - ε)) / (2ε)
        @test dP_anal ≈ dP_fd rtol=1e-5

        # shadow_factor_dσ² via FD
        μᵢ, μᵣ = 0.6, 0.7
        dS_anal = CM.shadow_factor_dσ²(μᵢ, μᵣ, σ²)
        dS_fd = (CM.shadow_factor(μᵢ, μᵣ, σ² + ε) - CM.shadow_factor(μᵢ, μᵣ, σ² - ε)) / (2ε)
        @test dS_anal ≈ dS_fd rtol=1e-4

        # whitecap_fraction_deriv via FD
        U = 7.0
        ε_U = 1e-6
        df_anal = CM.whitecap_fraction_deriv(U)
        df_fd = (CM.whitecap_fraction(U + ε_U) - CM.whitecap_fraction(U - ε_U)) / (2ε_U)
        @test df_anal ≈ df_fd rtol=1e-5
    end
end
