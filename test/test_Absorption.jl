# Test that a hitran fixed-width file is correctly parsed
@testset "read_hitran" begin

    println("Testing read_hitran...")

    ####
    #### These tests check that the correct values & types are extracted from the file,
    #### especially that molecule/isotope match and ν between ν_min and ν_max)
    ####

    CO2_test_file = "test_profiles/testCO2.data"

    test_ht = Absorption.read_hitran(CO2_test_file, mol=2, iso=1, ν_min=6000, ν_max=6400)
    @test test_ht.mol == [2, 2, 2, 2] && test_ht.mol[1] isa Int64
    @test test_ht.iso == [1, 1, 1, 1] && test_ht.iso[1] isa Int64
    @test test_ht.νᵢ == [6000.542970, 6286.403343, 6317.417493, 6380.824116] && test_ht.νᵢ[1] isa Float64
    @test test_ht.Sᵢ == [ 1.098E-28, 9.843E-30, 5.613E-27, 1.809E-30] && test_ht.Sᵢ[1] isa Float64
    @test test_ht.Aᵢ == [ 9.993e-08, 1.179e-08, 1.324e-05, 1.601e-02] && test_ht.Aᵢ[1] isa Float64
    @test test_ht.γ_air == [.0880, .0687, .0682, .0671] && test_ht.γ_air[1] isa Float64
    @test test_ht.γ_self == [0.118, 0.087, 0.081, 0.073] && test_ht.γ_self[1] isa Float64
    @test test_ht.E″ == [7.8043, 464.1717, 639.6004, 3798.2095] && test_ht.E″[1] isa Float64
    @test test_ht.n_air == [0.77, 0.76, 0.76, 0.73] && test_ht.n_air[1] isa Float64
    @test test_ht.δ_air == [-.004342, -.007362, -.007443, -.007669] && test_ht.δ_air[1] isa Float64
    @test test_ht.global_upper_quanta == ["       4 1 1 03", "       2 2 2 12", "       2 2 2 12", "       4 2 2 12"] && test_ht.global_upper_quanta[1] isa String
    @test test_ht.global_lower_quanta == ["       0 0 0 01", "       0 0 0 01", "       0 0 0 01", "       1 2 2 01"] && test_ht.global_lower_quanta[1] isa String
    @test test_ht.local_upper_quanta == ["               ", "               ", "               ", "               "] && test_ht.local_upper_quanta[1] isa String
    @test test_ht.local_lower_quanta == ["     Q  4e     ", "     Q 34e     ", "     R 40e     ", "     R 51f     "] && test_ht.local_lower_quanta[1] isa String
    @test test_ht.ierr == ["367774", "367764", "367764", "367774"] && test_ht.ierr[1] isa String
    @test test_ht.iref == ["2029 5 4 5 7", "2029 5 4 5 7", "2029 5 4 5 7", "2029 5 4 5 7"] && test_ht.iref[1] isa String
    @test test_ht.line_mixing_flag == [" ", " ", " ", " "] && test_ht.line_mixing_flag[1] isa String
    @test test_ht.g′ == [9.0, 69.0, 83.0, 105.0] && test_ht.g′[1] isa Float64
    @test test_ht.g″ == [9.0, 69.0, 81.0, 103.0] && test_ht.g″[1] isa Float64

    ####
    #### These tests check that the optional params work correctly
    ####

    # Not specifying the molecule #
    test_ht = Absorption.read_hitran(CO2_test_file, iso=1, ν_min=6000, ν_max=6400)
    @test test_ht.mol == [1, 2, 2, 2, 2]
    @test test_ht.iso == [1, 1, 1, 1, 1]
    @test test_ht.νᵢ == [6286.403343, 6000.542970, 6286.403343, 6317.417493, 6380.824116]
    @test test_ht.g″ == [69.0, 9.0, 69.0, 81.0, 103.0]

    # Not specifying the isotope #
    test_ht = Absorption.read_hitran(CO2_test_file, mol=2, ν_min=6000, ν_max=6400)
    @test test_ht.mol == [2, 2, 2, 2, 2]
    @test test_ht.iso == [2, 1, 1, 1, 1]
    @test test_ht.νᵢ == [6000.542970, 6000.542970, 6286.403343, 6317.417493, 6380.824116]
    @test test_ht.g″ == [9.0, 9.0, 69.0, 81.0, 103.0]

    # Not specifying the molecule # OR isotope #
    test_ht = Absorption.read_hitran(CO2_test_file, ν_min=6000, ν_max=6400)
    @test test_ht.mol == [1, 2, 2, 2, 2, 2]
    @test test_ht.iso == [1, 2, 1, 1, 1, 1]
    @test test_ht.νᵢ == [6286.403343, 6000.542970, 6000.542970, 6286.403343, 6317.417493, 6380.824116]
    @test test_ht.g″ == [69.0, 9.0, 9.0, 69.0, 81.0, 103.0]

    # Not specifying ν_min
    test_ht = Absorption.read_hitran(CO2_test_file, mol=2, iso=1, ν_max=6400)
    @test length(test_ht.mol) == 9

    # Not specifying ν_max
    test_ht = Absorption.read_hitran(CO2_test_file, mol=2, iso=1, ν_min=6000)
    @test length(test_ht.mol) == 7

    # Not specifying ν_min OR ν_max
    test_ht = Absorption.read_hitran(CO2_test_file, mol=2, iso=1)
    @test length(test_ht.mol) == 12

end

# Test that absorption cross sections are calculated correctly 
# (Not using pre-saved interpolator)

@testset "absorption_cross_section_hitran" begin

    println("Testing absorption_cross_section_hitran...")

    ####
    #### These tests check that for one molecule (CO2), over a temperature/
    #### pressure grid and in a specific band (6000<ν<6400), our Voigt 
    #### broadening implementation closely matches the HAPI output under
    #### the same conditions
    ####

    # Temperature and pressure grids. 
    # Note that these are pre-defined for testing -- the values are defined in generateHapiTests.py
    # To test against a different range, you must change the grids in that 
    # file and rerun it. Then, you can change it here. 

    pressures = 250:250:1250
    temperatures = 100:75:400 

    # Get the test data
    CO2_file = artifact("CO2")
    test_ht = Absorption.read_hitran(CO2_file, mol=2, iso=1, ν_min=6000, ν_max=6400)

    grid = 6000:0.01:6400;

    # Threshold -- our value must be within ϵ of the HAPI value
    ϵ = 3.6e-27

    # Create a HitranModel 
    model = make_hitran_model(test_ht, Voigt(), CEF=HumlicekWeidemann32SDErrorFunction())

    # Loop over every temperature/pressure combo and test that the results match HAPI
    @showprogress 1 "Testing HAPI equivalence (On CO2 Band)..." for temp in temperatures
        for pres in pressures
            jl_cs = absorption_cross_section(model, grid, pres, temp)
            py_cs = array_type(default_architecture)(readdlm("test_profiles/Voigt_CO2_T" * string(temp) * "_P" * string(pres) * ".csv"))
            Δcs = abs.(jl_cs - py_cs)
            @test maximum(Δcs) < ϵ
        end
    end

    #### 
    #### Now test HAPI equivalence with other molecules
    ####

    names = ["H2O", "CO2", "O3", "N2O", "CO"]
    molecules = [1, 2, 3, 4, 5]
    isotopes = [1, 1, 1, 1, 1]

    # Loop over every temperature/pressure combo and test that the results match HAPI
    # (Doing this for other molecules)
    @showprogress 1 "Testing HAPI equivalence (On Other Molecules)..." for name in names

        pres = 1000
        temp = 250

        # Get the test data
        test_ht = Absorption.read_hitran(artifact(name), iso=1, ν_min=6000, ν_max=6400)
        # Create a HitranModel 
        model = make_hitran_model(test_ht, Voigt(), CEF=HumlicekWeidemann32SDErrorFunction())

        jl_cs = absorption_cross_section(model, grid, pres, temp)
        py_cs = array_type(default_architecture)(readdlm("test_profiles/Voigt_" * name * "_T250_P1000.csv"))
        Δcs = abs.(jl_cs - py_cs)
        @test maximum(Δcs) < ϵ
    end
end

# Test that absorption cross sections if calculated using wavenumbers 

@testset "absorption_cross_section_wavelengths" begin

    println("Testing absorption_cross_section_wavelengths...")

    hitran_data = read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=6000, ν_max=6400)
    model = make_hitran_model(hitran_data, Voigt(), wing_cutoff = 40, CEF=HumlicekWeidemann32SDErrorFunction(), architecture=CPU())
    grid = 6000:0.01:6400
    @test absorption_cross_section(model, grid, 1000.1, 296.1) ≈ reverse(absorption_cross_section(model, reverse(1e7 ./ collect(grid)), 1000.1, 296.1, wavelength_flag=true));
end

# Test that absorption cross sections are calculated correctly 
# using a new interpolator

@testset "absorption_cross_section_interpolator" begin

    println("Testing absorption_cross_section_interpolator...")

    # Get the test data
    test_ht = Absorption.read_hitran(artifact("CO2"), iso=1, ν_min=6000, ν_max=6400)

    # Pressure and temperature grids
    pressures = 250:250:1250
    temperatures = 100:75:400 

    # Wavelength grid
    ν_grid = 6000:0.01:6400

    # Create the Interpolation Model
    interp_model = make_interpolation_model(test_ht, Voigt(), ν_grid, pressures, 
                                            temperatures, CEF=HumlicekWeidemann32SDErrorFunction()) 

    # Threshold -- our value must be within ϵ of the HAPI value
    ϵ = 3.6e-27

    # Loop over every temperature/pressure combo and test that the results match HAPI
    for temp in temperatures
        for pres in pressures
            jl_cs = absorption_cross_section(interp_model, ν_grid, pres, temp)
            py_cs = readdlm("test_profiles/Voigt_CO2_T" * string(temp) * "_P" * string(pres) * ".csv")
            Δcs = abs.(jl_cs - py_cs)
            @test maximum(Δcs) < ϵ
        end
    end

end 

# Test that checks whether the cross-section auto-differentiation works

@testset "absorption_cross_section_autodiff" begin

    println("Testing absorption_cross_section_autodiff...")
    
    # Load HITRAN data
    hitran_data = read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=6000, ν_max=6400)

    # Create the model with parameters
    model = make_hitran_model(hitran_data, Voigt(), wing_cutoff = 40, CEF=HumlicekWeidemann32SDErrorFunction(), architecture=CPU())

    # Compute the cross-section with autodifferentiation
    value, derivs = absorption_cross_section(model, 6000:0.01:6400, 1000.1, 296.1, autodiff=true);

end