@testset "planck_spectra" begin

    ν_grid = collect((1e7/778):0.015:(1e7/755))
    T = 3777

    println("Verifying wavelength/wavenumber correctness...")
    
    planck_ν = planck_spectrum(T, ν_grid)[:,2]
    planck_λ = planck_spectrum(T, reverse(1e7 ./ ν_grid), wavelength_flag=true)[:,2]
    
    @test planck_ν ≈ reverse(planck_λ)

end