# Test the wigner 3-j symbol calculations
@testset "wigner3j" begin

    # Meta-parameters
    j_max = 300     # Range of j-values to test
    N = 1000        # Number of randomized tests to perform

    # Compute Wigner matrices
    wigner_A, wigner_B = PhaseFunction.compute_wigner_values((2j_max + 1), j_max + 1, 2j_max + 1)

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