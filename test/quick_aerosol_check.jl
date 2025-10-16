"""
Quick TOMAS Data Inspector

Simple script to quickly check what's in the GEOSChem-TOMAS file.
Run interactively or as a script.

Usage:
    julia --project=. test/quick_aerosol_check.jl
"""

using NCDatasets
using Printf

# File and location
ncfile = "GEOSChem.Custom.20190702_0000z.nc4"
idx, idy, idf = 12, 12, 1  # Middle of grid, face 1

println("\n" * "="^60)
println("Quick TOMAS Aerosol Data Check")
println("="^60 * "\n")

ds = NCDataset(ncfile)

# Basic info
println("📍 Location: idx=$idx, idy=$idy, face=$idf")
println("   Lat: $(round(ds["lats"][idf,idy,idx], digits=2))°")
println("   Lon: $(round(ds["lons"][idf,idy,idx], digits=2))°\n")

# Aerosol species
species = ["DUST", "SS", "SF", "ECIL", "ECOB", "OCIL", "OCOB", "NK", "AW"]

println("🔬 Aerosol Species Summary:")
println("   Species | Bins | Max VMR    | Units")
println("   " * "-"^45)

for sp in species
    # Check first bin
    var = "SpeciesConcVV_$(sp)01"
    if haskey(ds, var)
        # Read all 15 bins at surface level (lev=1, which is BOA)
        max_vmr = 0.0
        for bin in 1:15
            vname = "SpeciesConcVV_$(sp)$(lpad(bin,2,'0'))"
            if haskey(ds, vname)
                data = ds[vname][idf, idy, idx, :, 1]
                max_vmr = max(max_vmr, maximum(data))
            end
        end
        
        units = ds[var].attrib["units"]
        @printf("   %-7s | %4d | %10.2e | %s\n", sp, 15, max_vmr, units)
    else
        println("   $sp: Not found")
    end
end

# Atmospheric structure
println("\n🌍 Atmospheric Structure:")
n_lev = ds.dim["lev"]
sp_val = ds["Met_PS2WET"][idf, idy, idx, 1]
T_surf = ds["Met_T"][idf, idy, idx, 1, 1]
T_top = ds["Met_T"][idf, idy, idx, end, 1]

println("   Vertical levels: $n_lev")
println("   Surface pressure: $(round(sp_val, digits=1)) hPa")
println("   Surface temp: $(round(T_surf, digits=1)) K")
println("   TOA temp: $(round(T_top, digits=1)) K")

# TOMAS bin info
println("\n📏 TOMAS-15 Bin Structure:")
diam_nm = [10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 
           20480, 40960, 81920, 163840, 327680]
radius_um = diam_nm ./ 2000
println("   Bin 1:  $(radius_um[1]) - $(radius_um[2]) μm radius")
println("   Bin 8:  $(radius_um[8]) - $(radius_um[9]) μm radius (middle)")
println("   Bin 15: $(radius_um[15]) - $(radius_um[16]) μm radius")

close(ds)

println("\n" * "="^60)
println("✅ Data check complete!")
println("\nTo explore visually, run:")
println("   julia --project=. test/test_aerosol_exploration.jl")
println("="^60 * "\n")
