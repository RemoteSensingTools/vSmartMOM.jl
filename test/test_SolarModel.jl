@testset "planck_spectra" begin

    ν_grid = collect((1e7/778):0.015:(1e7/755))
    T = 3777

    println("Verifying planck function correctness...")

    Ts = [290, 1000, 3777]
    peak_ν = [568.693, 1961.01, 7406.74]
    peak_λ = [9.99225, 2.89775, 0.76721]
    peak_L_ν = [138.636, 5684.38, 306284]
    peak_L_λ = [8.40098, 4095.81, 3.14829e6]

    for i in 1:length(Ts)

        T_curr = Ts[i]

        peak_ν_true = peak_ν[i]
        peak_L_true = peak_L_ν[i]

        peak_calc = findmax(planck_spectrum_wn(T_curr)[:,2])

        @test isapprox(peak_calc[1], peak_L_true, rtol=0.0001, atol=1)
        @test isapprox(peak_calc[2], peak_ν_true, rtol=0.0001, atol=1)
        @test isapprox(peak_L_λ[i], planck_spectrum_wl(T_curr, [peak_λ[i]])[1], rtol=0.0001, atol=1)

    end

end