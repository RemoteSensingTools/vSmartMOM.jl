using vSmartMOM
using vSmartMOM.Absorption


folder = "/net/fluo/data2/data/Rami_CS_database/"
# Pressure and temperature grids
pressures = 1e-5:25.0:1050
temperatures = 180.0:10:340 
Δν = 0.001

bands = ([1e7/456.0,  1e7/533.0],
         [1e7/538.0,  1e7/583.0],
         [1e7/646.0,  1e7/684.0],
         [1e7/837.0,  1e7/881.0],
         [1e7/1539.0, 1e7/1682.0],
         [1e7/2078.0, 1e7/2320.0])

molecules = ("CO2", "H2O", "O3", "N2O", "CO", "CH4", "O2")


for iMolec in eachindex(molecules)
    specFile = artifact(molecules[iMolec])
    for iBand in eachindex(bands)
        start = minimum(bands[iBand])-1
        stop  = maximum(bands[iBand])+1
        @show specFile, start, stop
        outFile = folder * "rami_spectroscopy_" * molecules[iMolec] * "_" * string(iBand) * ".jld2"
        @show outFile
        try
            test_ht = Absorption.read_hitran(specFile, ν_min=start-20, ν_max=stop+20)
            
            # Wavelength grid
            ν_grid = start:Δν:stop
            # Create the Interpolation Model
            interp_model = make_interpolation_model(test_ht, Voigt(), ν_grid, pressures, 
                                        temperatures, CEF=HumlicekWeidemann32SDErrorFunction()) 
            save_interpolation_model(interp_model, outFile)
        catch e
            println(e)
        end
    end
end

# Ozone separately: #####################
using DelimitedFiles
o3 = readdlm("//net/fluo/data2/data/Rami_CS_database/serdyuchenkogorshelev5digits.dat", skipstart=10000)
ν_grid    = o3[1,1]:0.01:o3[end,1]
# Remove <0 values:
o3[o3.<0] .= 0
cs_matrix = o3[:,end:-1:2]
d1,d2 = size(cs_matrix)
#cs_matrix = reshape(cs_matrix,(d1,1,d2))
cs_matrix_2 = zeros(d1,2,d2);
cs_matrix_2[:,1,:] = cs_matrix
cs_matrix_2[:,2,:] = cs_matrix

# Create Interpolation model from Interpolations.jl:
itp = interpolate(cs_matrix_2,  (BSpline(Linear()),BSpline(Constant()),BSpline(Linear())))
# Use flat extrapolation:
etpf = extrapolate(itp, Flat())
# Get the molecule and isotope numbers from the HitranTable
mol = -1
iso = -1
p_grid2 = 0.0:1100:1100.0
t_grid2 = 193.0:10:293.0
# Return the interpolation model with all the proper parameters
itp_model =  InterpolationModel(etpf, mol, iso,  ν_grid, p_grid2, t_grid2)
outFile = folder * "rami_spectroscopy_o3_all.jld2"
save_interpolation_model(itp_model, outFile)
