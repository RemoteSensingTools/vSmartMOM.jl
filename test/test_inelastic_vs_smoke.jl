# =============================================================================
# test_inelastic_vs_smoke.jl
#
# Smoke tests for the VS_0to1 / VS_1to0 Raman paths that were previously
# broken by UndefVarError:
#
#   Bug 1 (raman_atmo_prop.jl): compute_optical_Rayl! did not exist;
#          fixed to compute_optical_Rayl (no bang, returns value).
#
#   Bug 2 (inelastic_helper.jl): compute_optical_RS! for VS_1to0 called
#          get_greek_raman (discarding return) and compute_ϖ_Cabannes!
#          (undefined bang variant); both are now commented out, matching
#          the VS_0to1 pattern.
#
#   Bug 3 (apply_lineshape.jl / inelastic_helper.jl): cMassMolIE was
#          1.66053873e-30 (1000x too small — Doppler HWHM 31.6x too large);
#          corrected to 1.66053873e-27 kg.  Comment on cMassMol in
#          inelastic_helper.jl incorrectly said "grams"; fixed to "kilograms".
#
# These tests exercise the call paths that raised UndefVarError, verifying
# they now execute without error. They are executability smoke tests only,
# not physics-validation tests.
# =============================================================================

using Test
using vSmartMOM.InelasticScattering

@testset "VS inelastic smoke tests" begin

    # -------------------------------------------------------------------------
    # Shared setup: molecular constants and a minimal wavenumber grid.
    # Use the same central frequency as the O2-A band tests (~13170 cm⁻¹).
    # -------------------------------------------------------------------------
    ν̃ = 13170.0          # representative wavenumber [cm⁻¹]
    effT = 250.0          # representative atmospheric temperature [K]
    λ₀ = 1e7 / ν̃         # incident wavelength [nm]

    # N₂ and O₂ molecular constants (this also exercises getRamanAtmoConstants)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

    # A small wavenumber grid centred on the band (needed by compute_optical_RS!)
    # For VS_0to1 the shift is ~2330 cm⁻¹; put the grid well away from the
    # incident wavelength so the VRS lines actually fall inside.
    Δν_vs = 2330.0
    grid_vs = collect((ν̃ + Δν_vs - 100.0):0.5:(ν̃ + Δν_vs + 100.0))

    # -------------------------------------------------------------------------
    # Bug 1 fix: compute_optical_Rayl (no bang) must be callable and return a
    # positive scalar cross-section.
    # -------------------------------------------------------------------------
    @testset "Bug 1: compute_optical_Rayl is callable and returns a scalar" begin
        σ_rayl = InelasticScattering.compute_optical_Rayl(λ₀, n2, o2)
        @test isfinite(σ_rayl)
        @test σ_rayl > 0
    end

    # -------------------------------------------------------------------------
    # Bug 2 fix: compute_optical_RS! for VS_0to1 and VS_1to0 must not raise
    # UndefVarError.  Construct minimal placeholder instances (dispatch only —
    # the cross-section arrays are computed internally and returned, so the
    # struct field values don't affect the computation inside compute_optical_RS!).
    # -------------------------------------------------------------------------
    FT = Float64
    greek_placeholder = InelasticScattering.GreekCoefs(
        [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)])
    dummy_F₀   = zeros(FT, 1, 1)
    dummy_SIF₀ = zeros(FT, 1, 1)

    vs_0to1 = InelasticScattering.VS_0to1(
        n2           = n2,
        o2           = o2,
        greek_raman  = greek_placeholder,
        fscattRayl   = FT(1),
        ϖ_Cabannes   = FT(1),
        ϖ_λ₁λ₀      = zeros(FT, 1),
        i_λ₁λ₀      = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀    = zeros(FT, 1, 1),
        Z⁺⁺_λ₁λ₀    = zeros(FT, 1, 1),
        dτ₀          = FT(0),
        dτ₀_λ        = FT(0),
        k_Rayl_scatt = FT(1),
        n_Raman      = 0,
        F₀           = dummy_F₀,
        SIF₀         = dummy_SIF₀,
    )

    vs_1to0 = InelasticScattering.VS_1to0(
        n2           = n2,
        o2           = o2,
        greek_raman  = greek_placeholder,
        fscattRayl   = FT(1),
        ϖ_Cabannes   = FT(1),
        ϖ_λ₁λ₀      = zeros(FT, 1),
        i_λ₁λ₀      = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀    = zeros(FT, 1, 1),
        Z⁺⁺_λ₁λ₀    = zeros(FT, 1, 1),
        dτ₀          = FT(0),
        dτ₀_λ        = FT(0),
        k_Rayl_scatt = FT(1),
        n_Raman      = 0,
        F₀           = dummy_F₀,
        SIF₀         = dummy_SIF₀,
    )

    @testset "Bug 2: compute_optical_RS! VS_0to1 — no UndefVarError" begin
        # Returns (index_VRS, σ_VRS, index_RVRS, σ_RVRS)
        i_vrs, σ_vrs, i_rvrs, σ_rvrs =
            InelasticScattering.compute_optical_RS!(vs_0to1, grid_vs, λ₀, n2, o2)
        # All arrays must be finite
        @test all(isfinite, σ_vrs)
        @test all(isfinite, σ_rvrs)
    end

    @testset "Bug 2: compute_optical_RS! VS_1to0 — no UndefVarError" begin
        # Anti-Stokes: shift is negative — grid centred below incident ν̃
        grid_as = collect((ν̃ - Δν_vs - 100.0):0.5:(ν̃ - Δν_vs + 100.0))
        i_vrs, σ_vrs, i_rvrs, σ_rvrs =
            InelasticScattering.compute_optical_RS!(vs_1to0, grid_as, λ₀, n2, o2)
        @test all(isfinite, σ_vrs)
        @test all(isfinite, σ_rvrs)
    end

    # -------------------------------------------------------------------------
    # Bug 3 fix: cMassMolIE corrected from 1e-30 to 1e-27 kg.
    # We can't directly inspect the module-level constant after load, but we can
    # verify that apply_lineshape_! (which uses cMassMolIE) can be called
    # without producing non-finite output when given a single narrow line.
    # -------------------------------------------------------------------------
    @testset "Bug 3: apply_lineshape_! produces finite output (cMassMolIE fix)" begin
        # Single synthetic Raman transition at +100 cm⁻¹ from band centre
        Δνᵢ_in = [100.0]
        σᵢ_in  = [1.0e-30]   # plausible cross-section prefactor [cm²]
        out_grid = collect(-200.0:0.5:200.0)
        σ_out    = zeros(length(out_grid))
        pressure    = 500.0   # hPa
        temperature = 250.0   # K
        molMass     = 28.0    # N₂ molecular mass [dimensionless, in units of u]

        InelasticScattering.apply_lineshape_!(
            Δνᵢ_in, σᵢ_in,
            λ₀,
            out_grid, σ_out,
            pressure, temperature, molMass)

        @test all(isfinite, σ_out)
        @test any(>(0), σ_out)   # some non-zero contribution expected
    end

end
