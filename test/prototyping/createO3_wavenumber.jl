# Ozone separately: #####################
using DelimitedFiles
folder = "/net/fluo/data2/data/Rami_CS_database/"
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
#itp_model =  InterpolationModel(etpf, mol, iso,  ν_grid, p_grid2, t_grid2)
sitp = Interpolations.scale(itp, ν_grid, p_grid2, t_grid2)

ν  = 1e7/500.0:0.2:1e7/350
cs_matrix_nu  = sitp(1e7./ν, p_grid2, t_grid2)
itp_nu = interpolate(cs_matrix_nu,  (BSpline(Linear()),BSpline(Constant()),BSpline(Linear())))
etpf_nu = extrapolate(itp_nu, Flat())
itp_model2 =  vSmartMOM.Absorption.InterpolationModel(etpf_nu, mol, iso,  ν, p_grid2, t_grid2)

outFile = folder * "rami_spectroscopy_o3_all_nu.jld2"
vSmartMOM.Absorption.save_interpolation_model(itp_model2, outFile)