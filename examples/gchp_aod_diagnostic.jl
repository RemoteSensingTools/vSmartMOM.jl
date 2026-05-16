using vSmartMOM

gchp_path = get(ENV, "VSMARTMOM_GCHP_FILE",
                "/home/cfranken/data/GEOSChem.Custom.20190702_0000z.nc4")
out_path = get(ENV, "VSMARTMOM_GCHP_AOD_OUT",
               joinpath(pwd(), "gchp_aod_20190702_0000z.nc4"))

isfile(gchp_path) || error("GCHP file not found: $gchp_path")

write_gchp_aod_diagnostic(out_path, gchp_path;
                          wavelengths_um = [0.760, 1.600, 2.200],
                          aerosol_scheme = :auto,
                          ri_round_digits = tryparse(Int, get(ENV, "VSMARTMOM_AOD_RI_ROUND_DIGITS", "")))

println(out_path)
