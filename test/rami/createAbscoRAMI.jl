using vSmartMOM.Absorption
using Unitful

# Pressure and temperature grids
pressures = 1e-10:25.0:1050
temperatures = 180.0:5:340 
Δλ = 0.1

bands = ([1e7/456.0, 1e7/533.0],
        [1e7/538.0, 1e7/583.0],
        [1e7/646.0, 1e7/684.0],
        [1e7/837.0, 1e7/881.0],
        [1e7/1539.0, 1e7/1682.0],
        [1e7/2078.0, 1e7/2320.0])

molecules = ("CO2", "H2O", "O3", "N2O", "CO", "CH4", "O2")

for iMolec in eachindex(molecules)
    specFile = artifact(molecules[iMolec])
    for iBand in eachindex(bands)
        start = minimum(bands[iBand])-1
        stop  = maximum(bands[iBand])+1
        @show specFile, start, stop
        outFile = "rami_spectroscopy_" * molecules[iMolec] * "_" * string(iBand) * ".jld2"
        @show outFile
        try
            test_ht = Absorption.read_hitran(specFile, ν_min=start-20, ν_max=stop+20)
            
            # Wavelength grid
            ν_grid = start:Δλ:stop
            # Create the Interpolation Model
            interp_model = make_interpolation_model(test_ht, Voigt(), ν_grid, pressures, 
                                        temperatures, CEF=HumlicekWeidemann32SDErrorFunction()) 
            save_interpolation_model(interp_model, )
        catch 
            println("No hitran?")
        end
    end
end

# Get the test data
#CO2_file = artifact("CO2")
#test_ht = Absorption.read_hitran(CO2_file, mol=2, iso=1, ν_min=6000, ν_max=6400)

#grid = 6000:0.01:6400;