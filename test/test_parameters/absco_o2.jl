using Test

o2 = Dataset("/net/fluo/data2/data/ABSCO_CS_Database/v5.1_final/o2_v51.hdf")
cs_table    = o2["Gas_07_Absorption"]
wavenumber  = o2["Wavenumber"][:]
Temperature = o2["Temperature"][:]
Pressure    = o2["Pressure"][:]
Broadener   = o2["Broadener_01_VMR"][:]

i = 10
j = 64
PP = Pressure[j]/100
TT = Temperature[i,j]

hitran_data      = read_hitran(artifact("O2"), iso=1);
absorption_model = make_hitran_model(hitran_data, Voigt(), CEF=HumlicekWeidemann32SDErrorFunction(), vmr=0.0, wing_cutoff=150)
hitran_cs = absorption_cross_section(absorption_model, wavenumber, PP, TT)

a = vSmartMOM.loadAbsco("/net/fluo/data2/data/ABSCO_CS_Database/v5.1_final/o2_v51.hdf")
# Define interpolation table grid:
ν_grid = a.ν[1]:0.01:a.ν[end]
pressures = PP-10:10:PP+10.1
temperatures = TT-10:10:TT+10.1

model_interp = make_interpolation_model(a,  ν_grid, pressures, temperatures)
abscoInter_cs = absorption_cross_section(model_interp, ν_grid, PP, TT)

@test 1e20*abscoInter_cs ≈ 1e20*cs_table[:,1,i,j]
