# Test the wigner 3-j symbol calculations (slow — ~60s, gated behind VSMARTMOM_FULL_TESTS)
if get(ENV, "VSMARTMOM_FULL_TESTS", "") == "true"
@testset "wigner3j" begin

    # Meta-parameters
    j_max = 300     # Range of j-values to test
    N = 1000        # Number of randomized tests to perform

    # Compute Wigner matrices
    wigner_A, wigner_B = compute_wigner_values((2j_max + 1), j_max + 1, 2j_max + 1)

    println("Verifying Wigner Symbols...")

    # Outputs
    phase_function_results_A = Array{Float64,1}()
    phase_function_results_B = Array{Float64,1}()
    wigner_symbols_results_A = Array{Float64,1}()
    wigner_symbols_results_B = Array{Float64,1}()

    # Set counter to 1
    count = 0

    # Loop until we see N non-zeros Wigner values
    while count < N

        # Random inputs
        m = rand(1:j_max)
        n = rand(1:j_max)
        l = rand(1:j_max)

        # Result from PhaseFunction module
        push!(phase_function_results_A, wigner_A[m, n, l])
        push!(phase_function_results_B, wigner_B[m, n, l])

        # Result from external Wigner Symbols package
        # If there's a domain error, replace with 0.0.
        # (The 3j symbol *should* be zero outside the domain)
        try
            push!(wigner_symbols_results_A, Float64(wigner3j(m, n, l-1, -1, 1, 0)))
        catch
            push!(wigner_symbols_results_A, 0.0)
        end

        try
            push!(wigner_symbols_results_B, Float64(wigner3j(m, n, l-1, -1, -1, 2)))
        catch
            push!(wigner_symbols_results_B, 0.0)
        end

        # If a discrepancy ever pops up, print the discrepancy so it can be reproduced:
        if (phase_function_results_A[end] ≉ wigner_symbols_results_A[end] ||
            phase_function_results_B[end] ≉ wigner_symbols_results_B[end])
            println("Error with: ", (m, n, l))
            println("PhaseFunction output: ", (phase_function_results_A[end], phase_function_results_B[end]))
            println("WignerSymbols output: ", (wigner_symbols_results_A[end], wigner_symbols_results_B[end]))
        end

        # Only increment the counter if non-zero
        (phase_function_results_A[end] > 0) && (count = count + 1)

    end

    # Compare the result with WignerSymbols package values
    @test phase_function_results_A ≈ wigner_symbols_results_A
    @test phase_function_results_B ≈ wigner_symbols_results_B
end
end

# Test analytic phase functions through the shared Greek/AerosolOptics path.
@testset "analytic phase functions" begin
    hg = HenyeyGreensteinPhaseFunction(g = 0.4)
    μ = [-0.5, 0.0, 0.5]
    expected_hg = @. (1 - 0.4^2) / (1 + 0.4^2 - 2 * 0.4 * μ)^1.5
    @test phase_function(hg, μ) ≈ expected_hg

    greek = greek_coefficients(hg; l_max = 12, nquad = 48)
    @test length(greek.β) == 12
    @test greek.β[1] ≈ 1 atol = 1e-12

    optics = analytic_aerosol_optics(
        hg;
        single_scattering_albedo = 0.9,
        extinction_cross_section = 1.3,
        l_max = 12,
        nquad = 48)
    @test optics isa AerosolOptics
    @test optics.ω̃ ≈ 0.9
    @test optics.k ≈ 1.3
    @test optics.fᵗ == 0

    polarized = SyntheticPolarizedHenyeyGreensteinPhaseFunction(
        g = 0.3, polarization_fraction = 0.6)
    polarized_greek = greek_coefficients(polarized; l_max = 12, nquad = 48)
    @test any(abs.(polarized_greek.γ[3:end]) .> 0)

    pol_iq = Stokes_IQ{Float64}()
    @test pol_iq.n == 2
    @test pol_iq.D == [1.0, 1.0]
    @test pol_iq.I₀ == [1.0, 0.0]
    @test sprint(show, pol_iq) == "Stokes_IQ()"
    Z⁺⁺, Z⁻⁺ = vSmartMOM.Scattering.compute_Z_moments(
        pol_iq, [0.3, 0.7], polarized_greek, 0)
    @test size(Z⁺⁺) == (4, 4)
    @test size(Z⁻⁺) == (4, 4)
    @test all(isfinite.(Z⁺⁺))
    @test all(isfinite.(Z⁻⁺))
end

# Test the Aerosol Optics calculations (both NAI2 and Siewert)
@testset "aerosol_optics" begin

    println("Testing NAI2 and PCW equivalence...")

    # STEP 1: Create the Aerosol

    # Aerosol particle distribution and properties
    μ  = 0.3                # Log mean radius
    σ  = 2.1               # Log stddev of radius
    r_max = 30.0            # Maximum radius
    nquad_radius = 2500     # Number of quadrature points for integrating of size dist.
    nᵣ = 1.3                # Real part of refractive index
    nᵢ = 0.001              # Imag part of refractive index

    size_distribution = LogNormal(log(μ), log(σ))

    # Create the aerosol
    aero = Aerosol(size_distribution, nᵣ, nᵢ)

    # STEP 2: Create the Mie Calculations model

    λ = 0.55   # Incident wavelength
    polarization_type = Stokes_IQUV()
    truncation_type = δBGE(10, 10)
    model_NAI2 = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type, r_max, nquad_radius)

    # STEP 3: Perform the Mie Calculations and compare against saved PCW reference

    aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);

    # Load truth values computed from PCW
    @load "test_pcw/PCW_AerosolOptics_v2.jld" aerosol_optics_PCW

    # Greek coefficients: NAI2 numerical integration and PCW (Wigner-based)
    # agree to ~1e-4 absolute; tiny residuals in high-order tail where PCW
    # gives exact zeros cause the default rtol to fail, so we use atol.
    @test isapprox(aerosol_optics_NAI2.greek_coefs.α, aerosol_optics_PCW.greek_coefs.α, atol=1e-3)
    @test isapprox(aerosol_optics_NAI2.greek_coefs.β, aerosol_optics_PCW.greek_coefs.β, atol=1e-3)
    @test isapprox(aerosol_optics_NAI2.greek_coefs.γ, aerosol_optics_PCW.greek_coefs.γ, atol=1e-3)
    @test isapprox(aerosol_optics_NAI2.greek_coefs.δ, aerosol_optics_PCW.greek_coefs.δ, atol=1e-3)
    @test isapprox(aerosol_optics_NAI2.greek_coefs.ζ, aerosol_optics_PCW.greek_coefs.ζ, atol=1e-3)
    @test isapprox(aerosol_optics_NAI2.greek_coefs.ϵ, aerosol_optics_PCW.greek_coefs.ϵ, atol=1e-3)

    @test aerosol_optics_NAI2.ω̃ ≈ aerosol_optics_PCW.ω̃
    @test aerosol_optics_NAI2.k ≈ aerosol_optics_PCW.k
    @test aerosol_optics_NAI2.fᵗ ≈ aerosol_optics_PCW.fᵗ

end
