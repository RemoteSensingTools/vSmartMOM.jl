using Revise
using Plots
using Pkg
# Pkg.activate(".")
using vSmartMOM
using vSmartMOM.Architectures
using vSmartMOM.Absorption


pressures = 0.0001:25:1150.0
temperatures = 160:10:360.0

# O2 ABSCO:
file = "/net/fluo/data2/data/ABSCO_CS_Database/v5.2_final/o2_v52.hdf"
file_out = "/net/fluo/data2/data/ABSCO_CS_Database/v5.2_final/o2_v52.jld2"
a = CoreRT.loadAbsco(file; scale=1.0);

ν_grid = a.ν[1]:0.01:a.ν[end]

model_interp_O2 = make_interpolation_model(a, ν_grid, pressures, temperatures)
save_interpolation_model(model_interp_O2, file_out)

# H2O ABSCO:
file = "/net/fluo/data2/data/ABSCO_CS_Database/v5.2_final/h2o_v52.hdf"
file_out = "/net/fluo/data2/data/ABSCO_CS_Database/v5.2_final/wh2o_v52.jld2"
a = CoreRT.loadAbsco(file; scale=1.0);
grid_weak = a.ν[a.ν.>6000 .&& a.ν.<12000]

ν_grid = grid_weak[1]:0.01:grid_weak[end]
model_interp_H2O = make_interpolation_model(a, ν_grid, pressures, temperatures)
save_interpolation_model(model_interp_H2O, file_out)

# CO2 ABSCO:
file = "/net/fluo/data2/data/ABSCO_CS_Database/v5.2_final/co2_v52.hdf"
file_out = "/net/fluo/data2/data/ABSCO_CS_Database/v5.2_final/wco2_v52.jld2"
a = CoreRT.loadAbsco(file; scale=1.0);
grid_weak = a.ν[a.ν.>6000 .&& a.ν.<12000]

ν_grid = grid_weak[1]:0.01:grid_weak[end]

model_interp_CO2 = make_interpolation_model(a, ν_grid, pressures, temperatures)
save_interpolation_model(model_interp_CO2, file_out)

# H2O ABSCO:
file = "/net/fluo/data2/data/ABSCO_CS_Database/v5.2_final/h2o_v52.hdf"
file_out = "/net/fluo/data2/data/ABSCO_CS_Database/v5.2_final/sh2o_v52.jld2"
a = CoreRT.loadAbsco(file; scale=1.0);
grid_strong = a.ν[a.ν.<6000]

ν_grid = grid_strong[1]:0.01:grid_strong[end]
model_interp_H2O = make_interpolation_model(a, ν_grid, pressures, temperatures)
save_interpolation_model(model_interp_H2O, file_out)

# CO2 ABSCO:
file = "/net/fluo/data2/data/ABSCO_CS_Database/v5.2_final/co2_v52.hdf"
file_out = "/net/fluo/data2/data/ABSCO_CS_Database/v5.2_final/sco2_v52.jld2"
a = CoreRT.loadAbsco(file; scale=1.0);
strong = a.ν[a.ν.<6000]

ν_grid = grid_strong[1]:0.01:grid_strong[end]

model_interp_CO2 = make_interpolation_model(a, ν_grid, pressures, temperatures)
save_interpolation_model(model_interp_CO2, file_out)