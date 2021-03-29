# Test the wigner 3-j symbol calculations
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

# Test the Aerosol Optics calculations (both NAI2 and Siewert)
@testset "aerosol_optics" begin

    println("Testing NAI2 and PCW equivalence...")

    # STEP 1: Create the Aerosol

    τ

    # Aerosol particle distribution and properties 
    μ  = 0.3                # Log mean radius
    σ  = 6.82               # Log stddev of radius
    r_max = 30.0            # Maximum radius
    nquad_radius = 2500     # Number of quadrature points for integrating of size dist.
    nᵣ = 1.3                # Real part of refractive index
    nᵢ = 0.001              # Imag part of refractive index
    size_distribution = LogNormal(log(μ), log(σ))

    # Create the aerosol
    aero = make_univariate_aerosol(size_distribution, r_max, nquad_radius, nᵣ, nᵢ)

    # STEP 2: Create the Mie Calculations model

    λ = 0.55   # Incident wavelength
    polarization_type = Stokes_IQUV()
    truncation_type = δBGE(10, 10)
    model_NAI2 = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type)

    ### 
    ### NOTE: Temporarily removing PCW tests because they are too heavy to run on Travis
    ### 

    # Get saved wigner matrices
    # ftp = FTP("ftp://fluo.gps.caltech.edu/XYZT_hitran/")
    # println("Downloading full Wigner values...")
    # download(ftp, "wigner_values.jld", "/tmp/wigner_values.jld");

    # println("Loading full Wigner values...")
    # wigner_A, wigner_B = load_wigner_values("/home/rjeyaram/RadiativeTransfer/src/Scattering/Mie/wigner_values.jld")
    # model_PCW = make_mie_model(PCW(), aero, λ, polarization_type, truncation_type, wigner_A, wigner_B)


    # STEP 3: Perform the Mie Calculations and compare the results

    aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);
    # aerosol_optics_PCW = compute_aerosol_optical_properties(model_PCW);

    # Load truth values computed from PCW
    @load "helper/PCW_AerosolOptics.jld" aerosol_optics_PCW

    @test aerosol_optics_NAI2 ≈ aerosol_optics_PCW

    println("Testing aerosol_optical autodiff...")

    # Test whether autodiff works
    compute_aerosol_optical_properties(model_NAI2, autodiff=true);

end