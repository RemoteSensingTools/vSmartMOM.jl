# Test the wigner 3-j symbol calculations
@testset "wigner3j" begin

    # Meta-parameters
    j_max = 1000
    N = 1000

    # Testing all three cases
    for m_set in [[-1,-1,0]]

        # Set counter to 1
        n = 0

        # Run 1000 non-zero tests that compare our wigner 3j symbols to the WignerSymbols pkg
        while true
            # Random inputs
            m = rand(1:j_max)
            n = rand(1:j_max)
            l = rand(1:j_max)

            # Result from PhaseFunction module 
            res = wigner!(m, n, l, m_set[1], m_set[2], mset[3])

            # Compare the result with WignerSymbols package value
            @test res ≈ Float64(wigner3j(m, n, l, -1, 1, 0))

            # Only increment if non-zero
            if res > 0 
                n = n + 1
            end

            # Break if we've reached 1000 successful tests
            n == N ? break : nothing
        end
    end

    # Testing -1, 1, 0 case

end