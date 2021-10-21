o2 = Dataset("/net/fluo/data2/data/ABSCO_CS_Database/v5.1_final/o2_v51.hdf")
cs_table    = o2["Gas_07_Absorption"]
wavenumber  = o2["Wavenumber"][:]
Temperature = o2["Temperature"][:]
Pressure    = o2["Pressure"][:]
Broadener   = o2["Broadener_01_VMR"][:]

i = 10
j = 64
hitran_data      = read_hitran(artifact("O2"), iso=1);
absorption_model = make_hitran_model(hitran_data, Voigt(), CEF=HumlicekWeidemann32SDErrorFunction(), vmr=0.0)
hitran_cs = absorption_cross_section(absorption_model, wavenumber, Pressure[j]/100, Temperature[i,j])

a = vSmartMOM.loadAbsco("/net/fluo/data2/data/ABSCO_CS_Database/v5.1_final/o2_v51.hdf")
ν_grid = a.ν[1]:0.01:a.ν[end]
pressures = 1:20:1101.0
temperatures = 160:10:360.0
model_interp = make_interpolation_model(a, Voigt(), ν_grid, pressures, temperatures)
abscoInter_cs = absorption_cross_section(model_interp, wavenumber, Pressure[j]/100, Temperature[i,j])


