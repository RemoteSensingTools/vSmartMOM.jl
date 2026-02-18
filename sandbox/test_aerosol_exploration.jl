"""
TOMAS Aerosol Data Exploration

This script explores the GEOSChem-TOMAS aerosol data structure to understand:
1. Size distributions per species
2. Vertical profiles
3. Data characteristics before implementing full integration

Run with:
    julia --project=. test/test_aerosol_exploration.jl
"""

using NCDatasets
using Plots
using Statistics
using Printf

# File path
ncfile = "GEOSChem.Custom.20190702_0000z.nc4"
@assert isfile(ncfile) "NetCDF file not found: $ncfile"

println("="^70)
println("TOMAS Aerosol Data Exploration")
println("="^70)
println()

# Open dataset
ds = NCDataset(ncfile)

# Test location (middle of grid, face 1)
idx, idy, idf = 12, 12, 1
println("📍 Location: idx=$idx, idy=$idy, face=$idf")

# Get coordinates
lat = ds["lats"][idf, idy, idx]
lon = ds["lons"][idf, idy, idx]
println("   Latitude:  $(round(lat, digits=2))°")
println("   Longitude: $(round(lon, digits=2))°")
println()

# TOMAS bin definitions (standard TOMAS-15)
# Dry diameter edges in nm
diam_edges_nm = [
    10.0, 20.0, 40.0, 80.0, 160.0, 320.0, 640.0, 1280.0,
    2560.0, 5120.0, 10240.0, 20480.0, 40960.0, 81920.0, 163840.0, 327680.0
]

# Convert to radius in μm
radius_edges_um = diam_edges_nm ./ 2000.0  # nm → μm, diameter → radius
radius_centers_um = sqrt.(radius_edges_um[1:end-1] .* radius_edges_um[2:end])

println("🔬 TOMAS-15 Size Bins:")
println("   Bin  | Radius Range (μm)     | Center (μm)")
println("   " * "-"^50)
for i in 1:15
    @printf("   %2d   | %8.4f - %8.4f | %8.4f\n", 
            i, radius_edges_um[i], radius_edges_um[i+1], radius_centers_um[i])
end
println()

# Aerosol species to explore
aerosol_types = [
    ("DUST", "Mineral Dust"),
    ("SS", "Sea Salt"),
    ("SF", "Sulfate"),
    ("ECIL", "BC Hydrophilic"),
    ("ECOB", "BC Hydrophobic"),
    ("OCIL", "OC Hydrophilic"),
    ("OCOB", "OC Hydrophobic"),
    ("NK", "Nitrate/Potassium"),
    ("AW", "Aerosol Water")
]

# Get pressure and temperature profiles
n_lev = ds.dim["lev"]
dp = ds["Met_DELP"][idf, idy, idx, :, 1]
sp = ds["Met_PS2WET"][idf, idy, idx, 1]

# Pressure at half-levels (boundaries): [Surface, ..., TOA]
# GCHP stores bottom-to-top
pressure_half = [sp; sp .- cumsum(dp)]
# Flip to TOA→BOA for vSmartMOM convention
pressure_half_toa2boa = reverse(pressure_half)

# Layer center pressures
pressure_centers = (pressure_half_toa2boa[1:end-1] .+ pressure_half_toa2boa[2:end]) ./ 2

# Temperature profile (flip to TOA→BOA)
temperature = reverse(ds["Met_T"][idf, idy, idx, :, 1])

println("🌍 Atmospheric Profile:")
println("   Levels: $n_lev")
println("   Surface Pressure: $(round(sp, digits=1)) hPa")
println("   Top Pressure: $(round(minimum(pressure_half), digits=3)) hPa")
println("   Temperature range: $(round(minimum(temperature), digits=1)) - $(round(maximum(temperature), digits=1)) K")
println()

# Function to read a full species profile (all 15 bins)
function read_species_profile(ds, species_prefix, idx, idy, idf, n_lev)
    # Preallocate: 15 bins × n_lev layers
    concentration = zeros(15, n_lev)
    
    for bin_idx in 1:15
        var_name = "SpeciesConcVV_$(species_prefix)$(lpad(bin_idx, 2, '0'))"
        if haskey(ds, var_name)
            # VMR [mol/mol dry], stored BOA→TOA, flip to TOA→BOA
            concentration[bin_idx, :] = reverse(ds[var_name][idf, idy, idx, :, 1])
        end
    end
    
    return concentration  # [15 bins × n_lev layers]
end

println("📊 Reading Aerosol Species Profiles...")
println()

# Store all species data
species_data = Dict{String, Matrix{Float64}}()

for (prefix, name) in aerosol_types
    conc = read_species_profile(ds, prefix, idx, idy, idf, n_lev)
    species_data[prefix] = conc
    
    # Statistics
    total_conc = sum(conc)
    max_conc = maximum(conc)
    
    if max_conc > 1e-30
        println("   ✓ $name ($prefix):")
        println("      Total concentration: $((@sprintf "%.2e" total_conc)) mol/mol")
        println("      Max concentration:   $((@sprintf "%.2e" max_conc)) mol/mol")
        
        # Find layer with max concentration
        max_idx = argmax(vec(conc))
        bin_max, lev_max = divrem(max_idx - 1, n_lev) .+ (1, 1)
        println("      Peak at: Bin $(bin_max) (r=$(round(radius_centers_um[bin_max], digits=3)) μm), " *
                "Level $(lev_max) (p=$(round(pressure_centers[lev_max], digits=1)) hPa)")
    else
        println("   ○ $name ($prefix): negligible")
    end
end

close(ds)

println()
println("="^70)
println("Creating Visualizations...")
println("="^70)
println()

# Create output directory
mkpath("test/aerosol_exploration_output")

# Plot 1: Vertical profiles of total concentration per species
println("📈 Plot 1: Vertical Profiles (Total Concentration)")
p1 = plot(
    xlabel = "Total Concentration [mol/mol]",
    ylabel = "Pressure [hPa]",
    yflip = true,
    yscale = :log10,
    xscale = :log10,
    title = "Aerosol Vertical Profiles",
    legend = :bottomright,
    size = (800, 600)
)

for (prefix, name) in aerosol_types
    conc = species_data[prefix]
    total_per_layer = sum(conc, dims=1)[1, :]  # Sum over all bins
    
    if maximum(total_per_layer) > 1e-30
        plot!(p1, total_per_layer, pressure_centers, 
              label = name, linewidth = 2)
    end
end

savefig(p1, "test/aerosol_exploration_output/01_vertical_profiles.png")
println("   Saved: 01_vertical_profiles.png")

# Plot 2: Size distributions at different altitudes
println("📈 Plot 2: Size Distributions at Different Altitudes")

# Select interesting pressure levels (boundary layer, free troposphere, upper troposphere)
interesting_levels = [
    findfirst(p -> p < 900, pressure_centers),   # ~900 hPa (boundary layer)
    findfirst(p -> p < 500, pressure_centers),   # ~500 hPa (mid troposphere)
    findfirst(p -> p < 200, pressure_centers),   # ~200 hPa (upper troposphere)
]
interesting_levels = filter(!isnothing, interesting_levels)

# Plot major species
major_species = ["DUST", "SS", "SF", "ECIL"]

for species_prefix in major_species
    conc = species_data[species_prefix]
    
    if maximum(conc) < 1e-30
        continue
    end
    
    p2 = plot(
        xlabel = "Radius [μm]",
        ylabel = "Concentration [mol/mol]",
        xscale = :log10,
        yscale = :log10,
        title = "Size Distribution: $(species_prefix)",
        legend = :topright,
        size = (800, 600)
    )
    
    for lev_idx in interesting_levels
        conc_at_level = conc[:, lev_idx]
        p_label = round(pressure_centers[lev_idx], digits=0)
        
        plot!(p2, radius_centers_um, conc_at_level,
              label = "$(p_label) hPa",
              marker = :circle,
              linewidth = 2)
    end
    
    savefig(p2, "test/aerosol_exploration_output/02_size_dist_$(species_prefix).png")
    println("   Saved: 02_size_dist_$(species_prefix).png")
end

# Plot 3: 2D heatmap (altitude vs size) for each major species
println("📈 Plot 3: Altitude-Size Heatmaps")

for species_prefix in major_species
    conc = species_data[species_prefix]
    
    if maximum(conc) < 1e-30
        continue
    end
    
    # Take log10 for better visualization, avoid log(0)
    conc_log = log10.(conc .+ 1e-35)
    
    p3 = heatmap(
        1:15,  # Bin number
        pressure_centers,
        conc_log',
        xlabel = "Bin Number",
        ylabel = "Pressure [hPa]",
        title = "$(species_prefix): log₁₀(Concentration)",
        yflip = true,
        yscale = :log10,
        colorbar_title = "log₁₀[mol/mol]",
        size = (800, 600),
        color = :viridis
    )
    
    savefig(p3, "test/aerosol_exploration_output/03_heatmap_$(species_prefix).png")
    println("   Saved: 03_heatmap_$(species_prefix).png")
end

# Plot 4: Comparison of all species at a single altitude
println("📈 Plot 4: Species Comparison at Boundary Layer")

# Boundary layer level (~900 hPa)
bl_level = findfirst(p -> p < 900, pressure_centers)
if !isnothing(bl_level)
    p4 = plot(
        xlabel = "Radius [μm]",
        ylabel = "Concentration [mol/mol]",
        xscale = :log10,
        yscale = :log10,
        title = "Size Distributions at $(round(pressure_centers[bl_level], digits=0)) hPa",
        legend = :topright,
        size = (800, 600)
    )
    
    for (prefix, name) in aerosol_types
        conc = species_data[prefix]
        conc_at_bl = conc[:, bl_level]
        
        if maximum(conc_at_bl) > 1e-30
            plot!(p4, radius_centers_um, conc_at_bl,
                  label = name,
                  marker = :circle,
                  linewidth = 2)
        end
    end
    
    savefig(p4, "test/aerosol_exploration_output/04_all_species_comparison.png")
    println("   Saved: 04_all_species_comparison.png")
end

# Plot 5: Total aerosol mass/number vs altitude
println("📈 Plot 5: Column Integrated Aerosol")

p5 = plot(
    xlabel = "Column Integrated [mol/mol × layers]",
    ylabel = "Pressure [hPa]",
    yflip = true,
    yscale = :log10,
    xscale = :log10,
    title = "Cumulative Aerosol Distribution",
    legend = :bottomright,
    size = (800, 600)
)

for (prefix, name) in aerosol_types
    conc = species_data[prefix]
    
    # Cumulative sum from TOA to surface
    cumulative = cumsum(sum(conc, dims=1)[1, :])
    
    if maximum(cumulative) > 1e-30
        plot!(p5, cumulative, pressure_centers,
              label = name,
              linewidth = 2)
    end
end

savefig(p5, "test/aerosol_exploration_output/05_cumulative_profile.png")
println("   Saved: 05_cumulative_profile.png")

println()
println("="^70)
println("✅ Exploration Complete!")
println("="^70)
println()
println("Outputs saved to: test/aerosol_exploration_output/")
println()
println("Key Findings:")
println("  • Size range: $(round(radius_edges_um[1], digits=5)) - $(round(radius_edges_um[end], digits=2)) μm")
println("  • Number of bins: 15")
println("  • Vertical levels: $n_lev")
println("  • Active species: $(count(prefix -> maximum(species_data[prefix]) > 1e-30, first.(aerosol_types)))")
println()
println("Next Steps:")
println("  1. Review the generated plots")
println("  2. Identify which species/sizes dominate")
println("  3. Decide on fitting strategy (if needed)")
println("  4. Implement optical property calculations")
println()
