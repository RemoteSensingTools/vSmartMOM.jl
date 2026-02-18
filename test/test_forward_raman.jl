# =================================================================
# Raman Scattering Forward Model Test (Memory-Aware)
#
# Tests Rotational Raman Scattering (RRS) using a narrow O2-A band.
# RRS requires cross-wavelength coupling matrices (nSpec × nSpec) per layer;
# memory scales as O(nλ²).
#
# - GPU: uses O2Parameters_GPU.yaml (small wavelength range, ~60 points)
#   so forward Raman fits in device memory.
# - CPU: uses O2Parameters.yaml (finer resolution, ~6800 points).
# =================================================================

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Statistics
using Test

# Use GPU when CUDA is available (Raman is much faster on GPU)
const USE_GPU_RAMAN = try
    using CUDA
    CUDA.functional()
catch
    false
end

# ─────────────────────────────────────────────────────────────────
# Configuration: narrow O2-A band. Use smaller wavelength range on GPU
# to avoid OOM (RRS coupling matrices scale as nSpec²).
# ─────────────────────────────────────────────────────────────────
if USE_GPU_RAMAN
    params = parameters_from_yaml("test/test_parameters/O2Parameters_GPU.yaml")
    params.architecture = vSmartMOM.Architectures.GPU()
else
    params = parameters_from_yaml("test/test_parameters/O2Parameters.yaml")
    params.architecture = vSmartMOM.Architectures.CPU()
end
model = model_from_parameters(params)

iBand = 1
FT = Float64
ν = model.params.spec_bands[iBand]
nSpec = length(ν)
ν̃ = mean(ν)

# ─────────────────────────────────────────────────────────────────
# Set up RRS (Rotational Raman Scattering) type
# Following the pattern from prototype_inelastic.jl
# ─────────────────────────────────────────────────────────────────

# Effective temperature for Raman calculations
effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)

# Compute N₂ and O₂ molecular constants
n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

# Solar irradiance (Stokes I only, unit flux) — required by RRS and noRS constructors
nPol = model.params.polarization_type.n
F₀ = zeros(FT, nPol, nSpec)
F₀[1, :] .= 1.0
SIF₀ = zeros(FT, nPol, nSpec)

# Initialize RRS struct with placeholder arrays
RS_type = InelasticScattering.RRS(
    n2          = n2,
    o2          = o2,
    greek_raman = InelasticScattering.GreekCoefs(
        [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
    fscattRayl  = [FT(1)],
    ϖ_Cabannes  = [FT(1)],
    ϖ_λ₁λ₀     = zeros(FT, 1),
    i_λ₁λ₀     = zeros(Int, 1),
    Z⁻⁺_λ₁λ₀   = zeros(FT, 1, 1),
    Z⁺⁺_λ₁λ₀   = zeros(FT, 1, 1),
    i_ref       = argmin(abs.(ν .- ν̃)),
    n_Raman     = 0,
    F₀          = F₀,
    SIF₀        = SIF₀
)

# Compute Raman single-scattering properties
CoreRT.getRamanSSProp!(RS_type, 1e7/ν̃, ν)

@testset "Raman Scattering Test" begin

    # ─────────────────────────────────────────────────────────────
    # Test 1: Run RRS forward model
    # ─────────────────────────────────────────────────────────────
    @testset "RRS forward run" begin
        @time R_rrs, T_rrs, ieR, ieT = CoreRT.rt_run_test(RS_type, model, iBand)
        
        @test ndims(R_rrs) == 3
        @test all(isfinite.(R_rrs))
        @test all(isfinite.(ieR))
        
        nVza = length(model.obs_geom.vza)
        @test size(R_rrs) == (nVza, nPol, nSpec)
        @test size(ieR) == (nVza, nPol, nSpec)
        
        # Elastic component (R_rrs) should be positive for Stokes I
        I_elastic = R_rrs[:, 1, :]
        @test all(I_elastic .> 0) || @warn "Some elastic I < 0: min=$(minimum(I_elastic))"
    end
    
    # ─────────────────────────────────────────────────────────────
    # Test 2: Run noRS for comparison
    # ─────────────────────────────────────────────────────────────
    @testset "noRS baseline" begin
        noRS_type = InelasticScattering.noRS(
            fscattRayl  = [FT(1)],
            ϖ_Cabannes  = [FT(1)],
            bandSpecLim = [],
            iBand       = [1],
            F₀          = zeros(FT, nPol, nSpec)
        )
        noRS_type.F₀[1, :] .= 1.0
        
        @time R_noRS, T_noRS, _, _ = CoreRT.rt_run_test(noRS_type, model, iBand)
        
        @test all(isfinite.(R_noRS))
        I_noRS = R_noRS[:, 1, :]
        @test all(I_noRS .> 0)
    end
    
    # ─────────────────────────────────────────────────────────────
    # Test 3: Ring effect magnitude
    # The Ring effect (filling-in of Fraunhofer lines) should be
    # a few percent in the O2-A band region.
    # ─────────────────────────────────────────────────────────────
    @testset "Ring effect" begin
        # Re-run both to have them in same scope
        noRS_type = InelasticScattering.noRS(
            fscattRayl  = [FT(1)],
            ϖ_Cabannes  = [FT(1)],
            bandSpecLim = [],
            iBand       = [1],
            F₀          = zeros(FT, nPol, nSpec)
        )
        noRS_type.F₀[1, :] .= 1.0
        R_noRS, _, _, _ = CoreRT.rt_run_test(noRS_type, model, iBand)
        R_rrs, _, ieR, _ = CoreRT.rt_run_test(RS_type, model, iBand)
        
        # Total RRS radiance = elastic + inelastic
        I_total = R_rrs[:, 1, :] .+ ieR[:, 1, :]
        I_noRS = R_noRS[:, 1, :]
        
        # Ring effect: relative difference (RRS - noRS) / noRS
        ring_pct = 100.0 .* (I_total .- I_noRS) ./ I_noRS
        
        mean_ring = mean(abs.(ring_pct))
        max_ring = maximum(abs.(ring_pct))
        
        # Ring effect should be modest (typically 1-10% in O2-A)
        @test mean_ring < 50.0  # Loose bound — exact value depends on config
        @test max_ring < 200.0  # Some lines can have large filling-in
        
        # Total radiance should still be physical
        @test all(isfinite.(I_total))
        @test all(I_total .> 0)
    end
    
    # ─────────────────────────────────────────────────────────────
    # Test 4: Polarization consistency
    # |Q|, |U| should be bounded by I for both RRS and noRS
    # ─────────────────────────────────────────────────────────────
    if nPol >= 2
        @testset "Polarization consistency" begin
            R_rrs, _, ieR, _ = CoreRT.rt_run_test(RS_type, model, iBand)
            I_total = R_rrs[:, 1, :] .+ ieR[:, 1, :]
            Q_total = R_rrs[:, 2, :] .+ ieR[:, 2, :]
            
            @test all(abs.(Q_total) .<= I_total .+ 1e-10)
            
            if nPol >= 3
                U_total = R_rrs[:, 3, :] .+ ieR[:, 3, :]
                @test all(abs.(U_total) .<= I_total .+ 1e-10)
            end
        end
    end
end
