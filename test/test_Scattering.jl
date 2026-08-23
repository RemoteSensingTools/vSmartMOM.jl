# Standalone guard: when this file is run directly (not included from runtests.jl),
# the imports from runtests.jl are missing. Add them here so standalone execution works.
if !@isdefined(phase_function)
    using Test
    using vSmartMOM
    using vSmartMOM.Scattering
    using Distributions
    using JLD2
end

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

    zero_lin = zeros(4, length(polarized_greek.β))
    lin_greek = linGreekCoefs(zero_lin, copy(zero_lin), copy(zero_lin),
                              copy(zero_lin), copy(zero_lin), copy(zero_lin))
    Z⁺⁺_lin, Z⁻⁺_lin, dZ⁺⁺, dZ⁻⁺ =
        vSmartMOM.Scattering.compute_Z_moments(
            pol_iq, [0.3, 0.7], polarized_greek, lin_greek, 0)
    @test Z⁺⁺_lin ≈ Z⁺⁺
    @test Z⁻⁺_lin ≈ Z⁻⁺
    @test size(dZ⁺⁺) == (4, 4, 4)
    @test size(dZ⁻⁺) == (4, 4, 4)
    @test all(iszero, dZ⁺⁺)
    @test all(iszero, dZ⁻⁺)
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

# Verify that truncate_phase(δBGE) forward-cone exclusion is applied.
# Use the saved PCW reference (a large, forward-peaked aerosol with ~760 Greek
# modes) and l_max=10.  The normal equations are then 750+ observations × 10
# unknowns — strongly overdetermined and numerically stable on every platform.
# Using l_max ≈ length(β) (as the original test did with a small aerosol)
# makes the system nearly rank-deficient after the forward-cone exclusion
# removes a handful of GL points, causing c₀ to blow up outside [0,1].
@testset "δBGE forward-cone exclusion" begin
    @load "test_pcw/PCW_AerosolOptics_v2.jld" aerosol_optics_PCW
    aerosol_optics_PCW = AerosolOptics(
        greek_coefs=aerosol_optics_PCW.greek_coefs,
        ω̃=aerosol_optics_PCW.ω̃, k=aerosol_optics_PCW.k,
        fᵗ=aerosol_optics_PCW.fᵗ)
    # aerosol_optics_PCW carries raw Greek coefficients (fᵗ=1 sentinel).
    l_max = 10   # << length(β) ≈ 760 → strongly overdetermined → stable

    aop_tr0  = Scattering.truncate_phase(δBGE(l_max, 0.0),  aerosol_optics_PCW)
    aop_tr10 = Scattering.truncate_phase(δBGE(l_max, 10.0), aerosol_optics_PCW)

    # SSA and extinction unchanged by truncation:
    @test aop_tr0.ω̃ ≈ aerosol_optics_PCW.ω̃
    @test aop_tr0.k  ≈ aerosol_optics_PCW.k
    # fᵗ = 1 - c₀ must be non-negative; with l_max=10 the 10-mode fit of the
    # large aerosol's 760-mode series gives fᵗ ≈ 0.33 even at Δ_angle=0.
    @test aop_tr0.fᵗ  ≥ -1e-10
    @test 0 ≤ aop_tr10.fᵗ ≤ 1
    # fᵗ changes when the exclusion cone is applied (proves iμ subset is used):
    @test aop_tr10.fᵗ > 0.1   # Δ_angle=10° should give non-trivial truncation for large aerosol
    println("  δBGE fᵗ: Δ_angle=0° → $(round(aop_tr0.fᵗ, sigdigits=6)), Δ_angle=10° → $(round(aop_tr10.fᵗ, sigdigits=6))")
end

@testset "δBGE uniformly normalizes polarized matrix elements" begin
    @load "test_pcw/PCW_AerosolOptics_v2.jld" aerosol_optics_PCW
    aerosol_optics_PCW = AerosolOptics(
        greek_coefs=aerosol_optics_PCW.greek_coefs,
        ω̃=aerosol_optics_PCW.ω̃, k=aerosol_optics_PCW.k,
        fᵗ=aerosol_optics_PCW.fᵗ)
    mod = δBGE(10, 10.0)
    μ, wμ = Scattering.gausslegendre(
        length(aerosol_optics_PCW.greek_coefs.β))
    outside = findall(x -> x < cosd(mod.Δ_angle), μ)
    original = Scattering.reconstruct_phase(
        aerosol_optics_PCW.greek_coefs, μ)

    for (truncator, lowconf) in ((Scattering.truncate_phase, false),
                                 (Scattering.truncate_phase_lowconf, true))
        truncated = truncator(mod, aerosol_optics_PCW)
        reconstructed = Scattering.reconstruct_phase(
            truncated.greek_coefs, μ)
        c₀ = 1 - truncated.fᵗ

        # The fit targets the original, unnormalised f12/f34. After removal
        # of the delta spike, the retained smooth matrix must therefore obey
        # c₀ Fᵗ ≈ F for every fitted element outside the forward cone.
        for field in (:f₁₂, :f₃₄)
            truth = getproperty(original, field)[outside]
            scaled = c₀ .* getproperty(reconstructed, field)[outside]
            unscaled = getproperty(reconstructed, field)[outside]
            floor = fill(1e-8, length(truth))
            denominator = max.(abs.(truth), floor)
            # Main solver objective: Σw (residual/y)^2. Low-confidence
            # solver applies W=diag(w/|y|) before least squares, hence W².
            weights = lowconf ? (wμ[outside] ./ denominator).^2 :
                                wμ[outside] ./ denominator.^2
            objective(x) = sum(weights .* (x .- truth).^2)
            @test objective(scaled) ≤ objective(unscaled)
        end
        @test truncated.greek_coefs.β[1] ≈ 1 atol=1e-10
    end
end

@testset "δBGE polarized normalization tangent matches finite difference" begin
    @load "test_pcw/PCW_AerosolOptics_v2.jld" aerosol_optics_PCW
    raw = AerosolOptics(greek_coefs=aerosol_optics_PCW.greek_coefs,
        ω̃=aerosol_optics_PCW.ω̃, k=aerosol_optics_PCW.k,
        fᵗ=aerosol_optics_PCW.fᵗ)
    g = raw.greek_coefs
    fields = (:α, :β, :γ, :δ, :ϵ, :ζ)
    directions = Dict(f => 0.01 .* getproperty(g, f) for f in fields)
    tangent(a) = vcat(reshape(a, 1, :), zeros(3, length(a)))
    lg = linGreekCoefs((tangent(directions[f]) for f in fields)...)
    raw_lin = linAerosolOptics(lin_greek_coefs=lg,
        ω̃̇=zeros(4), k̇=zeros(4), ḟᵗ=zeros(4))

    # Δ=0 makes the forward and legacy linearized fit domains identical.
    mod = δBGE(10, 0.0)
    _, lin_out = Scattering.truncate_phase(mod, raw, raw_lin)
    # Central-difference step. The δ-BGE normal equations are ill-conditioned,
    # so the FD noise floor scales as O(eps·cond/h) and even flips with
    # compiler flags (Pkg.test's --check-bounds=yes changes SIMD summation
    # order); with the tangent direction 0.01·g, h = 1e-2 gives a 1e-4
    # relative step — large enough to dominate that noise, small enough that
    # the O(h²) truncation stays far below the tolerance below.
    h = 1e-2
    function perturbed(sign)
        gp = GreekCoefs((getproperty(g, f) .+
                         sign * h .* directions[f] for f in fields)...)
        Scattering.truncate_phase(mod,
            AerosolOptics(greek_coefs=gp, ω̃=raw.ω̃, k=raw.k, fᵗ=raw.fᵗ))
    end
    plus, minus = perturbed(1), perturbed(-1)
    for (field, tangent_field) in ((:γ, :γ̇), (:ϵ, :ϵ̇))
        fd = (getproperty(plus.greek_coefs, field) .-
              getproperty(minus.greek_coefs, field)) ./ (2h)
        analytic = getproperty(lin_out.lin_greek_coefs, tangent_field)[1, :]
        @test analytic ≈ fd rtol=2e-4 atol=2e-7
    end
end
