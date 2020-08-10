using RadiativeTransfer.CrossSection
using Test
using DelimitedFiles
using Statistics
using Plots

# Test the Cross Section module
@testset "RadiativeTransfer.CrossSection" begin

    # Test that a hitran fixed-width file is correctly parsed
    @testset "read_hitran" begin

        ####
        #### These tests check that the correct values & types are extracted from the file,
        #### especially that molecule/isotope match and ν between ν_min and ν_max)
        ####

        test_ht = CrossSection.read_hitran("helper/testCO2.data", mol=2, iso=1, ν_min=6000, ν_max=6400)
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
        test_ht = CrossSection.read_hitran("helper/testCO2.data", iso=1, ν_min=6000, ν_max=6400)
        @test test_ht.mol == [1, 2, 2, 2, 2]
        @test test_ht.iso == [1, 1, 1, 1, 1]
        @test test_ht.νᵢ == [6286.403343, 6000.542970, 6286.403343, 6317.417493, 6380.824116]
        @test test_ht.g″ == [69.0, 9.0, 69.0, 81.0, 103.0]

        # Not specifying the isotope #
        test_ht = CrossSection.read_hitran("helper/testCO2.data", mol=2, ν_min=6000, ν_max=6400)
        @test test_ht.mol == [2, 2, 2, 2, 2]
        @test test_ht.iso == [2, 1, 1, 1, 1]
        @test test_ht.νᵢ == [6000.542970, 6000.542970, 6286.403343, 6317.417493, 6380.824116]
        @test test_ht.g″ == [9.0, 9.0, 69.0, 81.0, 103.0]

        # Not specifying the molecule # OR isotope #
        test_ht = CrossSection.read_hitran("helper/testCO2.data", ν_min=6000, ν_max=6400)
        @test test_ht.mol == [1, 2, 2, 2, 2, 2]
        @test test_ht.iso == [1, 2, 1, 1, 1, 1]
        @test test_ht.νᵢ == [6286.403343, 6000.542970, 6000.542970, 6286.403343, 6317.417493, 6380.824116]
        @test test_ht.g″ == [69.0, 9.0, 9.0, 69.0, 81.0, 103.0]

        # Not specifying ν_min
        test_ht = CrossSection.read_hitran("helper/testCO2.data", mol=2, iso=1, ν_max=6400)
        @test length(test_ht.mol) == 9

        # Not specifying ν_max
        test_ht = CrossSection.read_hitran("helper/testCO2.data", mol=2, iso=1, ν_min=6000)
        @test length(test_ht.mol) == 7

        # Not specifying ν_min OR ν_max
        test_ht = CrossSection.read_hitran("helper/testCO2.data", mol=2, iso=1)
        @test length(test_ht.mol) == 12

    end

    # Test that absorption cross sections are calculated correctly 
    # (Not using pre-saved interpolator)

    @testset "absorption_cross_section" begin

        ####
        #### These tests check that for one molecule (CO2), over a temperature/
        #### pressure grid and in a specific band (6000<ν<6400), our Voigt 
        #### broadening implementation closely matches the HAPI output under
        #### the same conditions
        ####

        # Temperature and pressure grids. 
        # Note that these are static -- the values are defined in generateHapiTests.py
        # To test against a different range, you must change the grids in that 
        # file and rerun it. Then, you can change it here. 

        test_ht = CrossSection.read_hitran("helper/CO2.data", ν_min=6000, ν_max=6400)

        temperatures = [100, 175, 250, 325, 400]
        pressures = [250, 500, 750, 1000, 1250]

        grid = collect(6000:0.01:6400);
        CEF = ErfcHumliErrorFunctionVoigt()

        # Threshold -- our value must be within ϵ of the HAPI value
        ϵ = 3.5e-27

        for temp in temperatures
            println(temp)
            for pres in pressures
                jl_cs = absorption_cross_section(Voigt(CEF),test_ht,grid,false,pres,temp,40,vmr=0)
                py_cs = readdlm("helper/Voigt_CO2_T" * string(temp) * "_P" * string(pres) * ".csv")
                Δcs = abs.(jl_cs - py_cs)
                @test maximum(Δcs) < ϵ
            end
        end

    end


end
