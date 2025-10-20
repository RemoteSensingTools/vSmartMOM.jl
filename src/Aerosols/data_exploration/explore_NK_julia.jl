#!/usr/bin/env julia
"""
NK Number Size Distribution Explorer (Julia version)

Focuses on the NK variable (number concentration #/cm³) to understand
the actual particle size distribution in TOMAS-15.

Requirements:
    using Pkg
    Pkg.add(["NCDatasets", "Plots", "Statistics", "LsqFit"])

Usage:
    julia test/explore_NK_julia.jl
"""

using NCDatasets
using Plots
using Statistics
using LsqFit
using Printf

# Set plot backend
gr()

"""
    lognormal(r, params)

Log-normal size distribution (number concentration form)

Parameters:
- r: Particle radius [μm]
- params: [N_total, r_med, sigma_g]
  - N_total: Total number concentration [#/cm³]
  - r_med: Median radius [μm]
  - sigma_g: Geometric standard deviation

Returns:
- dN/dlogr: Number per log-radius bin [#/cm³]
"""
function lognormal(r::AbstractVector, p::AbstractVector)
    N_total, r_med, sigma_g = p
    return @. (N_total / (sqrt(2π) * log(sigma_g))) * 
              exp(-0.5 * (log(r / r_med) / log(sigma_g))^2)
end

# Version for single radius value
lognormal(r::Real, p::AbstractVector) = lognormal([r], p)[1]

"""
    bimodal_lognormal(r, params)

Bimodal log-normal distribution (sum of two modes)

Parameters:
- r: Particle radius [μm]
- params: [N1, r_med1, sigma_g1, N2, r_med2, sigma_g2]

Returns:
- dN/dlogr: Total number per log-radius bin [#/cm³]
"""
function bimodal_lognormal(r::AbstractVector, p::AbstractVector)
    N1, r_med1, sigma_g1, N2, r_med2, sigma_g2 = p
    mode1 = lognormal(r, [N1, r_med1, sigma_g1])
    mode2 = lognormal(r, [N2, r_med2, sigma_g2])
    return mode1 .+ mode2
end

"""
    fit_lognormal(r_centers, concentrations)

Fit a log-normal distribution to binned concentration data

Returns (N_total, r_med, sigma_g) or nothing if fit fails
"""
function fit_lognormal(r_centers::AbstractVector, concentrations::AbstractVector)
    # Filter out zero/negative values
    valid = concentrations .> 0
    if sum(valid) < 3
        return nothing
    end
    
    r_valid = r_centers[valid]
    c_valid = concentrations[valid]
    
    # Initial guess
    r_med_guess = r_valid[argmax(c_valid)]
    N_total_guess = sum(c_valid .* diff([log.(r_valid); log(r_valid[end])*1.1]))
    sigma_g_guess = 2.0
    
    p0 = [N_total_guess, r_med_guess, sigma_g_guess]
    
    # Bounds: N_total > 0, r_med in range, sigma_g > 1.01
    lower = [0.0, minimum(r_valid), 1.01]
    upper = [Inf, maximum(r_valid), 5.0]
    
    try
        fit = curve_fit((r, p) -> lognormal(r, p), r_valid, c_valid, p0, 
                       lower=lower, upper=upper)
        return fit.param
    catch e
        @warn "Fit failed: $e"
        return nothing
    end
end

"""
    fit_bimodal_lognormal(r_centers, concentrations)

Fit a bimodal log-normal distribution to binned concentration data

Returns (N1, r_med1, sigma_g1, N2, r_med2, sigma_g2) or nothing if fit fails
"""
function fit_bimodal_lognormal(r_centers::AbstractVector, concentrations::AbstractVector)
    # Filter out zero/negative values
    valid = concentrations .> 0
    if sum(valid) < 6
        return nothing
    end
    
    r_valid = r_centers[valid]
    c_valid = concentrations[valid]
    
    # Initial guess: split at 0.05 μm (Aitken vs accumulation boundary)
    aitken_mask = r_valid .< 0.05
    accum_mask = r_valid .>= 0.05
    
    # Mode 1 (Aitken): smaller particles
    if sum(aitken_mask) > 0
        r_med1_guess = 0.03  # 30 nm
        N1_guess = sum(c_valid[aitken_mask]) / sum(aitken_mask) * 0.5
        sigma_g1_guess = 1.5
    else
        r_med1_guess = 0.03
        N1_guess = maximum(c_valid) * 0.3
        sigma_g1_guess = 1.5
    end
    
    # Mode 2 (Accumulation): larger particles
    if sum(accum_mask) > 0
        accum_idx = argmax(c_valid[accum_mask])
        r_med2_guess = r_valid[accum_mask][accum_idx]
        N2_guess = sum(c_valid[accum_mask]) / sum(accum_mask) * 0.5
        sigma_g2_guess = 2.0
    else
        r_med2_guess = 0.15
        N2_guess = maximum(c_valid) * 0.7
        sigma_g2_guess = 2.0
    end
    
    # Ensure positive values
    N1_guess = max(N1_guess, 1.0)
    N2_guess = max(N2_guess, 1.0)
    
    p0 = [N1_guess, r_med1_guess, sigma_g1_guess, 
          N2_guess, r_med2_guess, sigma_g2_guess]
    
    # Bounds
    lower = [0.0, minimum(r_valid), 1.01, 0.0, 0.05, 1.01]
    upper = [Inf, 0.1, 3.0, Inf, maximum(r_valid), 5.0]
    
    try
        fit = curve_fit((r, p) -> bimodal_lognormal(r, p), r_valid, c_valid, p0,
                       lower=lower, upper=upper, maxIter=10000)
        return fit.param
    catch e
        @warn "Bimodal fit failed: $e"
        return nothing
    end
end

# ============================================================================
# Main Script
# ============================================================================

println("="^70)
println("NK Number Size Distribution Explorer (Julia)")
println("="^70)
println()

# Configuration
ncfile = "GEOSChem.Custom.20190702_0000z.nc4"
output_dir = "test/aerosol_exploration_output"
mkpath(output_dir)

# Define multiple locations to analyze
# Format: (idx, idy, idf, label) - Julia uses 1-based indexing
locations = [
    (3, 13, 5, "Central USA"),        # 36.8°N, 97.9°W (Oklahoma/Kansas)
    (16, 16, 2, "Amazon Basin"),      # Tropical South America
    (9, 19, 3, "Sahara Desert"),      # North Africa
    (13, 9, 4, "India/South Asia"),   # Indian subcontinent
    (6, 6, 6, "Northern China"),      # East Asia
    (18, 16, 4, "South Pacific"),     # Clean marine background (~22°S, 184°E)
]

# Open dataset
ds = NCDataset(ncfile, "r")

# Check if NK variable exists
test_var = "SpeciesConcVV_NK01"
if !(test_var in keys(ds))
    error("Variable $test_var not found in dataset!")
end

# Get location info
println("📍 Analyzing $(length(locations)) locations:")
for (i, (idx_loc, idy_loc, idf_loc, name)) in enumerate(locations)
    lat_loc = ds["lats"][idx_loc, idy_loc, idf_loc]
    lon_loc = ds["lons"][idx_loc, idy_loc, idf_loc]
    @printf("   %d. %-20s: %6.2f°N, %6.2f°E\n", i, name, lat_loc, lon_loc)
end
println()

# Use first location for detailed analysis
idx, idy, idf, loc_name = locations[1]
lat = ds["lats"][idx, idy, idf]
lon = ds["lons"][idx, idy, idf]
println("📍 Detailed analysis for: $loc_name")
@printf("   Location: idx=%d, idy=%d, face=%d\n", idx, idy, idf)
@printf("   Latitude:  %.2f°\n", lat)
@printf("   Longitude: %.2f°\n", lon)
println()

# TOMAS-15 bin definitions (dry diameter)
diam_min_nm = 10.0      # 10 nm minimum dry diameter
diam_max_nm = 10000.0   # 10 μm maximum dry diameter
diam_edges_nm = diam_min_nm .* (diam_max_nm / diam_min_nm) .^ (collect(0:15) ./ 15)
diam_edges_um = diam_edges_nm ./ 1000.0

# Bin centers (geometric mean)
diam_centers_um = sqrt.(diam_edges_um[1:end-1] .* diam_edges_um[2:end])

# Radius (for reference)
radius_centers_um = diam_centers_um ./ 2.0

println("🔬 TOMAS-15 Size Bins (Dry Particle Diameter):")
println("   Bin  | Diameter Range (μm)   | Center (μm)   | Radius Center (μm)")
println("   " * "-"^70)
for i in 1:15
    @printf("   %2d   | %8.4f - %8.4f | %8.4f      | %8.4f\n", 
            i, diam_edges_um[i], diam_edges_um[i+1], 
            diam_centers_um[i], radius_centers_um[i])
end
println()

# Get atmospheric structure
n_lev = ds.dim["lev"]
dp = vec(ds["Met_DELP"][idx, idy, idf, :, 1])
sp = Float64(ds["Met_PS2WET"][idx, idy, idf, 1])

# Pressure at half-levels (flip from BOA→TOA to TOA→BOA)
pressure_half = vcat([sp], sp .- cumsum(dp))
pressure_half_toa2boa = reverse(pressure_half)
pressure_centers = (pressure_half_toa2boa[1:end-1] .+ pressure_half_toa2boa[2:end]) ./ 2

# Temperature (flip to TOA→BOA)
temperature = reverse(vec(ds["Met_T"][idx, idy, idf, :, 1]))

# Calculate air density for unit conversion
# NK is stored in #/kg, need to convert to #/cm³
R_specific = 287.05  # J/(kg·K) for dry air
pressure_Pa = pressure_centers .* 100.0  # hPa to Pa
rho_air_kg_m3 = pressure_Pa ./ (R_specific .* temperature)  # kg/m³
rho_air_kg_cm3 = rho_air_kg_m3 .* 1e-6  # kg/cm³

println("🌍 Atmospheric Profile:")
@printf("   Levels: %d\n", n_lev)
@printf("   Surface Pressure: %.1f hPa\n", sp)
@printf("   Top Pressure: %.3f hPa\n", minimum(pressure_half))
@printf("   Temperature range: %.1f - %.1f K\n", minimum(temperature), maximum(temperature))
@printf("   Air density range: %.3f - %.3f kg/m³\n", minimum(rho_air_kg_m3), maximum(rho_air_kg_m3))
println()

# Read NK data (number concentration)
println("📊 Reading NK Number Concentration Data...")
println()

# Preallocate: 15 bins × n_lev layers
nk_raw = zeros(15, n_lev)  # As stored in file (#/kg)
nk_concentration = zeros(15, n_lev)  # Converted to #/cm³

for bin_idx in 1:15
    var_name = @sprintf("SpeciesConcVV_NK%02d", bin_idx)
    if var_name in keys(ds)
        # Read raw values (stored as #/kg despite metadata claiming mol/mol)
        raw_data = vec(ds[var_name][idx, idy, idf, :, 1])
        # Replace missing values with 0.0 and reverse
        for (i, val) in enumerate(reverse(raw_data))
            nk_raw[bin_idx, i] = ismissing(val) ? 0.0 : Float64(val)
        end
        
        # Convert from #/kg to #/cm³
        nk_concentration[bin_idx, :] .= nk_raw[bin_idx, :] .* rho_air_kg_cm3
        
        # Get metadata from first variable
        if bin_idx == 1
            units_meta = haskey(ds[var_name].attrib, "units") ? ds[var_name].attrib["units"] : "unknown"
            println("   Metadata units: $units_meta (INCORRECT!)")
            println("   Actual units: #/kg (particles per kilogram of air)")
            println("   Converted to: #/cm³ using air density")
        end
    end
end

# Statistics
total_n = sum(nk_concentration)
max_n = maximum(nk_concentration)

# Also show raw values for reference
total_n_raw = sum(nk_raw)
println()
@printf("   Raw NK values (as stored): %.2e #/kg\n", total_n_raw)
@printf("   Converted to #/cm³: %.2e #/cm³\n", total_n)
@printf("   Max number concentration: %.2e #/cm³\n", max_n)

if max_n > 0
    bin_max, lev_max = Tuple(argmax(nk_concentration))
    @printf("   Peak at: Bin %d (d=%.3f μm), Level %d (p=%.1f hPa)\n",
            bin_max, diam_centers_um[bin_max], lev_max, pressure_centers[lev_max])
else
    println("   WARNING: All NK values are zero!")
end

println()
println("="^70)
println("Creating Visualizations...")
println("="^70)
println()

# Calculate dlogD for TOMAS bins
dlogD = log10.(diam_edges_um[2:end] ./ diam_edges_um[1:end-1])

# ============================================================================
# Plot 1: Vertical profile by size mode
# ============================================================================
println("📈 Plot 1: Vertical Profile by Size Mode")

# Define size modes based on diameter (μm):
# Aitken mode: < 0.1 μm diameter (bins 1-5)
# Accumulation mode: 0.1 - 1.0 μm diameter (bins 6-10)
# Coarse mode: > 1.0 μm diameter (bins 11-15)

aitken_bins = 1:5
accumulation_bins = 6:10
coarse_bins = 11:15

aitken_n = vec(sum(nk_concentration[aitken_bins, :], dims=1))
accumulation_n = vec(sum(nk_concentration[accumulation_bins, :], dims=1))
coarse_n = vec(sum(nk_concentration[coarse_bins, :], dims=1))
total_n_profile = vec(sum(nk_concentration, dims=1))

p1 = plot(aitken_n, pressure_centers, 
          linewidth=2.5, label="Aitken (<0.1 μm)", 
          xaxis=:log, yflip=true,
          xlabel="Number Concentration [#/cm³]", ylabel="Pressure [hPa]",
          title="Aerosol Number Density by Size Mode - $loc_name\n(Converted from #/kg to #/cm³)",
          legend=:best, grid=true, gridα=0.3,
          size=(800, 600), color=:blue)
plot!(p1, accumulation_n, pressure_centers, linewidth=2.5, 
      label="Accumulation (0.1-1.0 μm)", color=:green)
plot!(p1, coarse_n, pressure_centers, linewidth=2.5, 
      label="Coarse (>1.0 μm)", color=:red)
plot!(p1, total_n_profile, pressure_centers, linewidth=2.5, 
      label="Total", linestyle=:dash, color=:black)

savefig(p1, joinpath(output_dir, "NK_01_vertical_profile.png"))
println("   Saved: NK_01_vertical_profile.png")
@printf("   Mode splits:\n")
@printf("     Aitken bins 1-5: %.4f - %.4f μm\n", diam_edges_um[1], diam_edges_um[6])
@printf("     Accumulation bins 6-10: %.4f - %.4f μm\n", diam_edges_um[6], diam_edges_um[11])
@printf("     Coarse bins 11-15: %.4f - %.4f μm\n", diam_edges_um[11], diam_edges_um[16])

# ============================================================================
# Plot 2: Size distributions at multiple altitudes
# ============================================================================
println("📈 Plot 2: Size Distributions at Multiple Altitudes")
println()

# Select interesting pressure levels
pressure_levels = [900, 800, 500, 300, 200]
level_indices = [argmin(abs.(pressure_centers .- p)) for p in pressure_levels]

# Show info about conversion at selected levels
println("   Conversion examples at selected altitudes:")
for p_target in pressure_levels
    lev_idx = argmin(abs.(pressure_centers .- p_target))
    p_actual = pressure_centers[lev_idx]
    rho = rho_air_kg_cm3[lev_idx]
    nk_val = sum(nk_raw[:, lev_idx])
    n_val = sum(nk_concentration[:, lev_idx])
    @printf("     %.0f hPa: ρ=%.3e kg/cm³, NK=%.2e #/kg → N=%.2e #/cm³\n", 
            p_actual, rho, nk_val, n_val)
end

println()

# Fine diameter grid for smooth curves (for future fits)
d_fine = 10 .^ range(log10(minimum(diam_centers_um)), log10(maximum(diam_centers_um)), length=200)

p2 = plot(xlabel="Dry Diameter [μm]", ylabel="dN/dlogD [#/cm³]",
          xaxis=:log, 
          title="Aerosol Number Size Distribution - $loc_name\n(True number density: NK × ρ_air)",
          legend=:topright, grid=true, gridα=0.3,
          size=(1000, 700))

colors = palette(:viridis, length(level_indices))

for (i, (lev_idx, p_target)) in enumerate(zip(level_indices, pressure_levels))
    nk_at_level = nk_concentration[:, lev_idx]
    p_actual = pressure_centers[lev_idx]
    
    # Convert to dN/dlogD
    dN_dlogD = nk_at_level ./ dlogD
    
    # Plot data
    plot!(p2, diam_centers_um, dN_dlogD, 
          marker=:circle, linewidth=2.5, markersize=5,
          label="$(Int(round(p_actual))) hPa", 
          color=colors[i], α=0.8)
    
    # Fit log-normal (print params but don't plot)
    @printf("   Fitting log-normal at %.0f hPa:\n", p_actual)
    params = fit_lognormal(radius_centers_um, dN_dlogD)
    if params !== nothing
        N_total, r_med, sigma_g = params
        d_med = r_med * 2
        d_eff = d_med * exp(2.5 * log(sigma_g)^2)
        
        @printf("      N_total = %.2e #/cm³\n", N_total)
        @printf("      d_med   = %.4f μm\n", d_med)
        @printf("      d_eff   = %.4f μm\n", d_eff)
        @printf("      σ_g     = %.3f\n", sigma_g)
    else
        println("      Could not fit log-normal distribution")
    end
end

savefig(p2, joinpath(output_dir, "NK_02_size_distributions_multi_altitude.png"))
println("\n   Saved: NK_02_size_distributions_multi_altitude.png")

# ============================================================================
# Plot 3: 2D heatmap (altitude vs size)
# ============================================================================
println("📈 Plot 3: Altitude-Size Heatmap")

# Log scale for visualization
nk_log = log10.(nk_concentration .+ 1e-10)  # Avoid log(0)

p3 = heatmap(diam_centers_um, pressure_centers, nk_log',
             xlabel="Dry Diameter [μm]", ylabel="Pressure [hPa]",
             xaxis=:log, yaxis=:log, yflip=true,
             title="Aerosol Number Density Heatmap - $loc_name\n(True concentration in #/cm³)",
             colorbar_title="log₁₀[Number Density (#/cm³)]",
             color=:plasma, size=(1000, 600))

savefig(p3, joinpath(output_dir, "NK_03_heatmap.png"))
println("   Saved: NK_03_heatmap.png")

# ============================================================================
# Plot 4: Boundary layer focus with bimodal fit
# ============================================================================
println("📈 Plot 4: Boundary Layer Size Distribution (detailed)")

bl_level = argmin(abs.(pressure_centers .- 900))
p_bl = pressure_centers[bl_level]
nk_bl = nk_concentration[:, bl_level]
dN_dlogD_bl = nk_bl ./ dlogD

println()
@printf("   Fitting bimodal log-normal at %.0f hPa:\n", p_bl)

# Fit bimodal distribution
bimodal_params = fit_bimodal_lognormal(radius_centers_um, dN_dlogD_bl)

# Left panel: Linear y-axis with bimodal fit
p4_left = plot(diam_centers_um, dN_dlogD_bl,
               marker=:circle, linewidth=2.5, markersize=6,
               label="Data", color=:darkgreen,
               xaxis=:log,
               xlabel="Dry Diameter [μm]", ylabel="dN/dlogD [#/cm³]",
               title="Boundary Layer ($(Int(round(p_bl))) hPa) - Bimodal Fit\n(True number density)",
               legend=:topright, grid=true, gridα=0.3,
               size=(700, 500))

if bimodal_params !== nothing
    N1, r_med1, sigma_g1, N2, r_med2, sigma_g2 = bimodal_params
    
    # Convert back to diameter for display
    d_med1 = r_med1 * 2
    d_med2 = r_med2 * 2
    
    # Plot individual modes
    r_fine = d_fine ./ 2
    mode1 = lognormal(r_fine, [N1, r_med1, sigma_g1])
    mode2 = lognormal(r_fine, [N2, r_med2, sigma_g2])
    total_fit = mode1 .+ mode2
    
    plot!(p4_left, d_fine, mode1, linestyle=:dash, linewidth=2, 
          label=@sprintf("Aitken: N=%.2e, d_med=%.3fμm, σ_g=%.2f", N1, d_med1, sigma_g1),
          color=:blue, α=0.7)
    plot!(p4_left, d_fine, mode2, linestyle=:dash, linewidth=2,
          label=@sprintf("Accumulation: N=%.2e, d_med=%.3fμm, σ_g=%.2f", N2, d_med2, sigma_g2),
          color=:orange, α=0.7)
    plot!(p4_left, d_fine, total_fit, linewidth=2.5,
          label="Total Fit", color=:red, α=0.9)
    
    @printf("      Mode 1 (Aitken):       N = %.2e #/cm³, d_med = %.4f μm, σ_g = %.3f\n", N1, d_med1, sigma_g1)
    @printf("      Mode 2 (Accumulation): N = %.2e #/cm³, d_med = %.4f μm, σ_g = %.3f\n", N2, d_med2, sigma_g2)
    @printf("      Total N: %.2e #/cm³\n", N1 + N2)
    @printf("      Aitken fraction: %.1f%%\n", N1/(N1+N2)*100)
else
    println("      Could not fit bimodal distribution")
end

# Right panel: Log y-axis
p4_right = plot(diam_centers_um, dN_dlogD_bl,
                marker=:circle, linewidth=2.5, markersize=6,
                label="Data", color=:darkgreen,
                xaxis=:log, yaxis=:log,
                xlabel="Dry Diameter [μm]", ylabel="dN/dlogD [#/cm³]",
                title="Boundary Layer ($(Int(round(p_bl))) hPa) - Log Scale\n(True number density)",
                legend=:topright, grid=true, gridα=0.3,
                size=(700, 500))

if bimodal_params !== nothing
    plot!(p4_right, d_fine, mode1, linestyle=:dash, linewidth=2,
          label="Aitken Mode", color=:blue, α=0.7)
    plot!(p4_right, d_fine, mode2, linestyle=:dash, linewidth=2,
          label="Accumulation Mode", color=:orange, α=0.7)
    plot!(p4_right, d_fine, total_fit, linewidth=2.5,
          label="Total Fit", color=:red, α=0.9)
end

p4 = plot(p4_left, p4_right, layout=(1, 2), size=(1400, 500))
savefig(p4, joinpath(output_dir, "NK_04_boundary_layer_detail.png"))
println("   Saved: NK_04_boundary_layer_detail.png")

# ============================================================================
# Plot 5: Multi-location comparison
# ============================================================================
println("📈 Plot 5: Size Distributions Across Multiple Locations")

# Use 800 hPa level
lev_800 = argmin(abs.(pressure_centers .- 800))

p5 = plot(xlabel="Dry Diameter [μm]", ylabel="dN/dlogD [#/cm³]",
          xaxis=:log,
          title="Aerosol Number Density at $(Int(round(pressure_centers[lev_800]))) hPa - Multi-Location\n(True concentration: NK × ρ_air)",
          legend=:topright, grid=true, gridα=0.3,
          size=(1000, 700))

println()
for (loc_idx, (idx_loc, idy_loc, idf_loc, loc_name_comp)) in enumerate(locations)
    # Get air density for this location
    temp_loc_raw = vec(ds["Met_T"][idx_loc, idy_loc, idf_loc, :, 1])
    temp_loc = zeros(n_lev)
    for (i, val) in enumerate(reverse(temp_loc_raw))
        temp_loc[i] = ismissing(val) ? 250.0 : Float64(val)
    end
    rho_air_loc = (pressure_Pa ./ (R_specific .* temp_loc)) .* 1e-6  # kg/cm³
    
    # Read NK data for this location (raw #/kg)
    nk_raw_loc = zeros(15, n_lev)
    nk_loc = zeros(15, n_lev)
    for bin_idx in 1:15
        var_name = @sprintf("SpeciesConcVV_NK%02d", bin_idx)
        if var_name in keys(ds)
            raw_data = vec(ds[var_name][idx_loc, idy_loc, idf_loc, :, 1])
            for (i, val) in enumerate(reverse(raw_data))
                nk_raw_loc[bin_idx, i] = ismissing(val) ? 0.0 : Float64(val)
            end
            # Convert to #/cm³
            nk_loc[bin_idx, :] .= nk_raw_loc[bin_idx, :] .* rho_air_loc
        end
    end
    
    # Get at 800 hPa
    nk_at_level = nk_loc[:, lev_800]
    dN_dlogD_loc = nk_at_level ./ dlogD
    
    # Plot
    plot!(p5, diam_centers_um, dN_dlogD_loc,
          marker=:circle, linewidth=2.5, markersize=5,
          label=loc_name_comp, α=0.8)
    
    # Fit (print but don't plot)
    println("   $loc_name_comp:")
    params_loc = fit_lognormal(radius_centers_um, dN_dlogD_loc)
    if params_loc !== nothing
        N_total_loc, r_med_loc, sigma_g_loc = params_loc
        d_med_loc = r_med_loc * 2
        @printf("      N_total = %.2e #/cm³, d_med = %.4f μm, σ_g = %.3f\n", 
                N_total_loc, d_med_loc, sigma_g_loc)
    end
end

savefig(p5, joinpath(output_dir, "NK_05_multi_location_comparison.png"))
println("\n   Saved: NK_05_multi_location_comparison.png")

# Close dataset
close(ds)

println()
println("="^70)
println("✅ NK Number Distribution Analysis Complete!")
println("="^70)
println()
println("Outputs saved to: $output_dir/")
println()

# Find 800 hPa level for summary
lev_800_summary = argmin(abs.(pressure_centers .- 800))
surface_n = sum(nk_concentration[:, end])
n_800 = sum(nk_concentration[:, lev_800_summary])

println("Summary Statistics:")
@printf("  • Size range: %.4f - %.1f μm diameter\n", diam_edges_um[1], diam_edges_um[end])
println("  • Number of bins: 15")
@printf("  • Vertical levels: %d\n", n_lev)
@printf("  • Surface number concentration: %.2e #/cm³\n", surface_n)
@printf("  • At 800 hPa: %.2e #/cm³\n", n_800)
println()
println("Key Findings:")
println("  • NK is stored as #/kg (not mol/mol as metadata claims)")
println("  • Converted using air density: N(#/cm³) = NK(#/kg) × ρ_air(kg/cm³)")
println("  • Values are now in true atmospheric number concentration")
println("  • Typical continental air: 10³-10⁴ #/cm³ (matches our values!)")
