@testset "planck_spectra" begin

    ν_grid = collect((1e7/778):0.015:(1e7/755))
    T = 3777

    println("Verifying wavelength/wavenumber correctness...")
    
    planck_ν = planck_spectrum(T, ν_grid)
    planck_λ = planck_spectrum(T, reverse(1e7 ./ ν_grid), wavelength_flag=true)
    
    @test planck_ν ≈ reverse(planck_λ)

    println("Verifying planck law correctness...")

    Ts = [290, 1000, 3777]
    peak_ν = [568.693, 1961.01, 7406.74]
    peak_L = [138.636, 5684.38, 306284]

    for i in 1:length(Ts)

        T_curr = Ts[i]

        peak_ν_true = peak_ν[i]
        peak_L_true = peak_L[i]

        peak_calc = findmax(planck_spectrum(T_curr)[:,2])

        @test isapprox(peak_calc[1], peak_L_true, rtol=0.0001, atol=1)
        @test isapprox(peak_calc[2], peak_ν_true, rtol=0.0001, atol=1)
    end

end