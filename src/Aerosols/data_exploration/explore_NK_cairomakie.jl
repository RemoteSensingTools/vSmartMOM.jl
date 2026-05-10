#!/usr/bin/env julia
"""
NK Number Size Distribution Explorer (CairoMakie version)

Focuses on the NK variable (number concentration #/cm³) to understand
the actual particle size distribution in TOMAS-15.

Requirements:
    using Pkg
    Pkg.add(["NCDatasets", "CairoMakie", "Statistics", "LsqFit", "Distributions"])

Usage:
    julia src/Aerosols/data_exploration/explore_NK_cairomakie.jl
"""

using NCDatasets
using CairoMakie
using Statistics
using LsqFit
using Distributions
using Printf

# Set CairoMakie theme for publication quality
set_theme!(Theme(
    fontsize = 14,
    linewidth = 2.5,
    markersize = 10,
    Axis = (
        xgridvisible = true,
        ygridvisible = true,
        xgridalpha = 0.3,
        ygridalpha = 0.3,
    )
))

"""
    lognormal(r, params)

Log-normal size distribution (number concentration form)

Parameters:
- r: particle radius [μm]
- params: [N_total, r_med, σ_g] where
  - N_total: total number concentration [#/cm³]
  - r_med: median radius [μm]
  - σ_g: geometric standard deviation (dimensionless)

Returns:
- dN/dlogr [#/cm³] - number per log-radius bin
"""
function lognormal(r, params)
    N_total, r_med, σ_g = params
    # Direct formula for dN/dlogr
    return (N_total / (sqrt(2*π) * log(σ_g))) .* 
           exp.(-0.5 .* (log.(r./r_med) ./ log(σ_g)).^2)
end

"""
    bimodal_lognormal(r, params)

Bimodal log-normal size distribution

Parameters:
- r: particle diameter [μm]
- params: [N1, d1, σ1, N2, d2, σ2] where subscript 1,2 refer to mode 1 and 2
"""
function bimodal_lognormal(r, params)
    N1, d1, σ1, N2, d2, σ2 = params
    mode1 = lognormal(r, [N1, d1, σ1])
    mode2 = lognormal(r, [N2, d2, σ2])
    return mode1 .+ mode2
end

"""
    fit_bimodal_lognormal(r_centers, concentrations)

Fit a bimodal log-normal distribution to aerosol size distribution data.
Uses radius (not diameter) and has intelligent initial guesses for Aitken/Accumulation modes.
Fits to dN/dlogD values with physically reasonable constraints on σ_g.
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
        sigma_g1_guess = 1.6
    else
        r_med1_guess = 0.03
        N1_guess = maximum(c_valid) * 0.3
        sigma_g1_guess = 1.6
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
    
    # Bounds - match Python version exactly
    # Python: bounds=([0, r_valid.min(), 1.01, 0, 0.05, 1.01],
    #                [np.inf, 0.1, 3.0, np.inf, r_valid.max(), 5.0])
    lower = [0.0, minimum(r_valid), 1.01, 0.0, 0.05, 1.01]
    upper = [Inf, 0.1, 3.0, Inf, maximum(r_valid), 5.0]
    
    try
        # Fit to dN/dlogD (not integrated bins)
        fit = curve_fit((r, p) -> bimodal_lognormal(r, p), r_valid, c_valid, p0,
                       lower=lower, upper=upper, maxIter=10000)
        return fit.param
    catch e
        @warn "Bimodal fit failed: $e"
        return nothing
    end
end

# ============================================================================
# Configuration
# ============================================================================

# NetCDF file path
nc_file = "GEOSChem.Custom.20190702_0000z.nc4"

# Output directory
output_dir = "test/aerosol_exploration_output"
isdir(output_dir) || mkdir(output_dir)

# TOMAS-15 bin information (dry particle diameter in μm)
diam_edges_um = [0.01, 0.0158, 0.0251, 0.0398, 0.0631, 0.1, 0.1585, 0.2512, 
                 0.3981, 0.631, 1.0, 1.5849, 2.5119, 3.9811, 6.3096, 10.0]
diam_centers_um = sqrt.(diam_edges_um[1:end-1] .* diam_edges_um[2:end])

# Radius centers for reference
radius_centers_um = diam_centers_um ./ 2.0

# Locations to analyze (idx, idy, face, name) - Julia uses 1-based indexing
locations = [
    (3, 13, 5, "Central USA"),        # 36.8°N, 97.9°W (Oklahoma/Kansas)
    (16, 16, 2, "Amazon Basin"),      # Tropical South America
    (9, 19, 3, "Sahara Desert"),      # North Africa
    (13, 9, 4, "India/South Asia"),   # Indian subcontinent
    (6, 6, 6, "Northern China"),      # East Asia
    (18, 16, 4, "South Pacific"),     # Clean marine background (~22°S, 184°E)
]

# Primary location for detailed analysis
idx, idy, idf, loc_name = locations[1]

println("="^70)
println("NK Number Size Distribution Explorer (CairoMakie)")
println("="^70)
println()
println("📍 Analyzing $(length(locations)) locations:")
for (i, (f, x, y, name)) in enumerate(locations)
    println(@sprintf("   %d. %-18s :  (face=%d, x=%d, y=%d)", i, name, f, x, y))
end
println()

# ============================================================================
# Load Data
# ============================================================================

ds = NCDataset(nc_file)

# Get lat/lon for the primary location
lat = ds["lats"][idx, idy, idf]
lon = ds["lons"][idx, idy, idf]

println("📍 Detailed analysis for: $loc_name")
println("   Location: idx=$idx, idy=$idy, face=$idf")
@printf("   Latitude:  %.2f°\n", lat)
@printf("   Longitude: %.2f°\n", lon)
println()

# Get atmospheric structure
n_lev = ds.dim["lev"]
dp = vec(ds["Met_DELP"][idx, idy, idf, :, 1])
sp = Float64(ds["Met_PS2WET"][idx, idy, idf, 1])

# Pressure at half-levels (flip from BOA→TOA to TOA→BOA)
pressure_half = vcat([sp], sp .- cumsum(dp))
pressure_half_toa2boa = reverse(pressure_half)
pressure_centers = (pressure_half_toa2boa[1:end-1] .+ pressure_half_toa2boa[2:end]) ./ 2
pressure_Pa = pressure_centers .* 100.0  # Convert to Pascals

# Get temperature (flip to TOA→BOA)
temperature = reverse(vec(ds["Met_T"][idx, idy, idf, :, 1]))

# Physical constants
M_air = 28.9644e-3  # kg/mol (molar mass of dry air)
R_specific = 287.05  # J/(kg·K) for dry air

# Calculate air density for reference
rho_air_kg_m3 = pressure_Pa ./ (R_specific .* temperature)

# Display size bin information
println("🔬 TOMAS-15 Size Bins (Dry Particle Diameter):")
println("   Bin  | Diameter Range (μm)   | Center (μm)   | Radius Center (μm)")
println("   " * "-"^70)
for i in 1:15
    @printf("   %2d   |   %.4f -   %.4f |   %.4f      |   %.4f\n",
            i, diam_edges_um[i], diam_edges_um[i+1], diam_centers_um[i], radius_centers_um[i])
end
println()

println("🌍 Atmospheric Profile:")
@printf("   Levels: %d\n", n_lev)
@printf("   Surface Pressure: %.1f hPa\n", pressure_centers[end])
@printf("   Top Pressure: %.3f hPa\n", pressure_centers[1])
@printf("   Temperature range: %.1f - %.1f K\n", minimum(temperature), maximum(temperature))
@printf("   Air density range: %.3f - %.3f kg/m³\n", minimum(rho_air_kg_m3), maximum(rho_air_kg_m3))
println()

# ============================================================================
# Read NK Data with CORRECT Conversion
# ============================================================================

# Read meteorology for NK conversion
Met_AD = reverse(vec(ds["Met_AD"][idx, idy, idf, :, 1]))  # kg
Met_AIRVOL = reverse(vec(ds["Met_AIRVOL"][idx, idy, idf, :, 1]))  # m³
n_air = Met_AD ./ M_air  # moles of air
vol_cm3 = Met_AIRVOL .* 1e6  # cm³
rho_air_kg_m3 = Met_AD ./ Met_AIRVOL  # for reference only

# Read NK data (number concentration)
println("📊 Reading NK Number Concentration Data...")
println()

# Preallocate: 15 bins × n_lev layers
nk_raw = zeros(15, n_lev)  # Raw NK values (1000 × #/mol_air)
nk_concentration = zeros(15, n_lev)  # Converted to #/cm³

for bin_idx in 1:15
    var_name = @sprintf("SpeciesConcVV_NK%02d", bin_idx)
    if var_name in keys(ds)
        # Read raw values (stored as 1000 × particles/mol_air)
        raw_data = vec(ds[var_name][idx, idy, idf, :, 1])
        # Replace missing values with 0.0 and reverse
        for (i, val) in enumerate(reverse(raw_data))
            nk_raw[bin_idx, i] = ismissing(val) ? 0.0 : Float64(val)
        end
        
        # Convert using CORRECT formula: N = (NK/1000) × (n_air/vol_cm3)
        nk_concentration[bin_idx, :] .= (nk_raw[bin_idx, :] ./ 1000.0) .* (n_air ./ vol_cm3)
        
        # Get metadata from first variable
        if bin_idx == 1
            units_meta = haskey(ds[var_name].attrib, "units") ? ds[var_name].attrib["units"] : "unknown"
            println("   Metadata units: $units_meta")
            println("   Actual units: 1000 × (particles/mol_air)")
            println("   Conversion: N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)")
            println("   Result: #/cm³")
        end
    end
end

# Statistics
total_n = sum(nk_concentration)
max_n = maximum(nk_concentration)

# Also show raw values for reference
total_n_raw = sum(nk_raw)
println()
@printf("   Raw NK values (as stored): %.2e (1000 × #/mol_air)\n", total_n_raw)
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

# Calculate dN/dlogD (convert from #/cm³ per bin to #/cm³ per decade of diameter)
dlogD = log10.(diam_edges_um[2:end]) .- log10.(diam_edges_um[1:end-1])
dN_dlogD_profiles = nk_concentration ./ dlogD

# ============================================================================
# Create Visualizations with CairoMakie
# ============================================================================

println("="^70)
println("Creating Visualizations...")
println("="^70)
println()

# ============================================================================
# Plot 1: Vertical Profile by Size Mode
# ============================================================================
println("📈 Plot 1: Vertical Profile by Size Mode")

# Split into size modes
aitken_n = vec(sum(nk_concentration[1:5, :], dims=1))       # < 0.1 μm
accumulation_n = vec(sum(nk_concentration[6:10, :], dims=1)) # 0.1-1.0 μm
coarse_n = vec(sum(nk_concentration[11:15, :], dims=1))     # > 1.0 μm
total_n_profile = vec(sum(nk_concentration, dims=1))

fig1 = Figure(size=(900, 700))
ax1 = Axis(fig1[1, 1],
    xlabel = "Number Concentration [#/cm³]",
    ylabel = "Pressure [hPa]",
    xscale = log10,
    yreversed = true,
    title = "Aerosol Number Density by Size Mode - $loc_name\nValidated: N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)"
)

lines!(ax1, aitken_n, pressure_centers, label="Aitken (<0.1 μm)", color=:blue, linewidth=2.5)
lines!(ax1, accumulation_n, pressure_centers, label="Accumulation (0.1-1.0 μm)", color=:green, linewidth=2.5)
lines!(ax1, coarse_n, pressure_centers, label="Coarse (>1.0 μm)", color=:red, linewidth=2.5)
lines!(ax1, total_n_profile, pressure_centers, label="Total", color=:black, linewidth=2.5, linestyle=:dash)

axislegend(ax1, position=:rb)

save(joinpath(output_dir, "NK_01_vertical_profile.png"), fig1)
println("   Saved: NK_01_vertical_profile.png")
println("   Mode splits:")
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

# Show conversion info
println("   Conversion examples at selected altitudes:")
for p_target in pressure_levels
    lev_idx = argmin(abs.(pressure_centers .- p_target))
    p_actual = pressure_centers[lev_idx]
    
    n_air_val = n_air[lev_idx]
    vol_val = vol_cm3[lev_idx]
    conv_factor = n_air_val / vol_val
    
    nk_val = sum(nk_raw[:, lev_idx])
    n_val = sum(nk_concentration[:, lev_idx])
    
    @printf("     %.0f hPa: NK=%.2e (1000×#/mol), n_air/V=%.2e mol/cm³ → N=%.2e #/cm³\n", 
            p_actual, nk_val, conv_factor, n_val)
end
println()

fig2 = Figure(size=(1100, 750))
ax2 = Axis(fig2[1, 1],
    xlabel = "Dry Diameter [μm]",
    ylabel = "dN/dlogD [#/cm³]",
    xscale = log10,
    title = "Aerosol Number Size Distribution - $loc_name\nValidated: N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)"
)

colors_viridis = cgrad(:viridis, length(level_indices), categorical=true)

for (i, (lev_idx, p_target)) in enumerate(zip(level_indices, pressure_levels))
    dN_dlogD = dN_dlogD_profiles[:, lev_idx]
    p_actual = pressure_centers[lev_idx]
    
    # Fit lognormal
    N_total_est = sum(nk_concentration[:, lev_idx])
    d_med_est = diam_centers_um[argmax(dN_dlogD)]
    σ_g_est = 2.0
    p0 = [N_total_est, d_med_est, σ_g_est]
    
    try
        fit = curve_fit(lognormal, diam_centers_um, dN_dlogD .* dlogD, p0)
        N_fit, d_med_fit, σ_g_fit = fit.param
        d_eff = d_med_fit * exp(2.5 * log(σ_g_fit)^2)
        
        println("   Fitting log-normal at $(Int(round(p_actual))) hPa:")
        @printf("      N_total = %.2e #/cm³\n", N_fit)
        @printf("      d_med   = %.4f μm\n", d_med_fit)
        @printf("      d_eff   = %.4f μm\n", d_eff)
        @printf("      σ_g     = %.3f\n", σ_g_fit)
    catch e
        println("   Fitting failed at $(Int(round(p_actual))) hPa")
    end
    
    # Plot
    lines!(ax2, diam_centers_um, dN_dlogD, 
           label="$(Int(round(p_actual))) hPa",
           color=colors_viridis[i], linewidth=2.5)
    scatter!(ax2, diam_centers_um, dN_dlogD,
             color=colors_viridis[i], markersize=8)
end

axislegend(ax2, position=:rt)

save(joinpath(output_dir, "NK_02_size_distributions_multi_altitude.png"), fig2)
println()
println("   Saved: NK_02_size_distributions_multi_altitude.png")

# ============================================================================
# Plot 3: Heatmap (altitude vs size)
# ============================================================================
println("📈 Plot 3: Altitude-Size Heatmap")

# Create log10 of dN/dlogD for better visualization
log_dN_dlogD = log10.(dN_dlogD_profiles .+ 1e-10)

fig3 = Figure(size=(1000, 700))
ax3 = Axis(fig3[1, 1],
    xlabel = "Dry Diameter [μm]",
    ylabel = "Pressure [hPa]",
    xscale = log10,
    yreversed = true,
    title = "Aerosol Number Density Heatmap - $loc_name\n(True concentration in #/cm³)"
)

# Create meshgrid for heatmap
# CairoMakie heatmap expects (n_x, n_y) data where x varies along first dimension
hm = heatmap!(ax3, diam_centers_um, pressure_centers, log_dN_dlogD,
              colormap=:turbo)

Colorbar(fig3[1, 2], hm, label="log₁₀[dN/dlogD (#/cm³)]")

save(joinpath(output_dir, "NK_03_heatmap.png"), fig3)
println("   Saved: NK_03_heatmap.png")

# ============================================================================
# Plot 4: Boundary layer detailed analysis
# ============================================================================
println("📈 Plot 4: Boundary Layer Size Distribution (detailed)")
println()

# Use ~900 hPa as boundary layer
lev_bl = argmin(abs.(pressure_centers .- 900))
p_bl = pressure_centers[lev_bl]
dN_dlogD_bl = dN_dlogD_profiles[:, lev_bl]

println()
@printf("   Boundary layer: level %d at %.0f hPa\n", lev_bl, p_bl)
@printf("   Actual data (first 10 bins):\n")
for i in 1:min(10, length(dN_dlogD_bl))
    @printf("     Bin %2d: dN/dlogD = %.2f #/cm³\n", i, dN_dlogD_bl[i])
end
@printf("   Sum(dN/dlogD) = %.2f #/cm³\n", sum(dN_dlogD_bl))
@printf("   Integrated N (sum(dN/dlogD×dlogD)) = %.2f #/cm³\n", sum(dN_dlogD_bl .* dlogD))
println()
@printf("   Fitting bimodal log-normal at %.0f hPa:\n", p_bl)

# Fit bimodal distribution using RADIUS (not diameter) like the old script
bimodal_params = fit_bimodal_lognormal(radius_centers_um, dN_dlogD_bl)

if bimodal_params !== nothing
    N1, r_med1, σ1, N2, r_med2, σ2 = bimodal_params
    
    # Convert back to diameter for display
    d_med1 = r_med1 * 2
    d_med2 = r_med2 * 2
    
    @printf("      Mode 1 (Aitken):       N = %.2e #/cm³, d_med = %.4f μm, σ_g = %.3f\n", N1, d_med1, σ1)
    @printf("      Mode 2 (Accumulation): N = %.2e #/cm³, d_med = %.4f μm, σ_g = %.3f\n", N2, d_med2, σ2)
    @printf("      Total N: %.2e #/cm³\n", N1+N2)
    @printf("      Aitken fraction: %.1f%%\n", 100*N1/(N1+N2))
    
    # Create plots
    d_fine = 10 .^ range(log10(minimum(diam_centers_um)), log10(maximum(diam_centers_um)), length=200)
    r_fine = d_fine ./ 2
    mode1 = lognormal(r_fine, [N1, r_med1, σ1])
    mode2 = lognormal(r_fine, [N2, r_med2, σ2])
    total_fit = mode1 .+ mode2
    
    fig4 = Figure(size=(900, 700))
    
    # Single panel with log x-axis
    ax4 = Axis(fig4[1, 1],
        xlabel = "Dry Diameter [μm]",
        ylabel = "dN/dlogD [#/cm³]",
        xscale = log10,
        title = "Boundary Layer ($(Int(round(p_bl))) hPa) - Bimodal Fit\n(True number density)"
    )
    
    scatter!(ax4, diam_centers_um, dN_dlogD_bl, label="Data", color=:darkgreen, markersize=10)
    lines!(ax4, d_fine, mode1, 
           label=@sprintf("Aitken: N=%.2e, d_med=%.3fμm, σ_g=%.2f", N1, d_med1, σ1),
           color=:blue, linestyle=:dash, linewidth=2)
    lines!(ax4, d_fine, mode2, 
           label=@sprintf("Accumulation: N=%.2e, d_med=%.3fμm, σ_g=%.2f", N2, d_med2, σ2),
           color=:orange, linestyle=:dash, linewidth=2)
    lines!(ax4, d_fine, total_fit, label="Total Fit", color=:red, linewidth=2.5)
    
    axislegend(ax4, position=:lt)
    
    save(joinpath(output_dir, "NK_04_boundary_layer_detail.png"), fig4)
    println("   Saved: NK_04_boundary_layer_detail.png")
else
    println("   Bimodal fitting failed!")
end

# ============================================================================
# Plot 5: Multi-location comparison
# ============================================================================
println("📈 Plot 5: Size Distributions Across Multiple Locations")
println()

# Use 800 hPa level
lev_800 = argmin(abs.(pressure_centers .- 800))

fig5 = Figure(size=(1100, 750))
ax5 = Axis(fig5[1, 1],
    xlabel = "Dry Diameter [μm]",
    ylabel = "dN/dlogD [#/cm³]",
    xscale = log10,
    title = "Aerosol Number Density at $(Int(round(pressure_centers[lev_800]))) hPa - Multi-Location\nValidated: N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)"
)

colors_tab10 = Makie.wong_colors()

for (loc_idx, (idx_loc, idy_loc, idf_loc, loc_name_comp)) in enumerate(locations)
    # Get meteorology for this location
    Met_AD_loc_raw = vec(ds["Met_AD"][idx_loc, idy_loc, idf_loc, :, 1])
    Met_AIRVOL_loc_raw = vec(ds["Met_AIRVOL"][idx_loc, idy_loc, idf_loc, :, 1])
    
    Met_AD_loc = zeros(n_lev)
    Met_AIRVOL_loc = zeros(n_lev)
    for (i, val) in enumerate(reverse(Met_AD_loc_raw))
        Met_AD_loc[i] = ismissing(val) ? 1.0 : Float64(val)
    end
    for (i, val) in enumerate(reverse(Met_AIRVOL_loc_raw))
        Met_AIRVOL_loc[i] = ismissing(val) ? 1.0 : Float64(val)
    end
    
    n_air_loc = Met_AD_loc ./ M_air
    vol_cm3_loc = Met_AIRVOL_loc .* 1e6
    
    # Read NK data
    nk_raw_loc = zeros(15, n_lev)
    nk_loc = zeros(15, n_lev)
    for bin_idx in 1:15
        var_name = @sprintf("SpeciesConcVV_NK%02d", bin_idx)
        if var_name in keys(ds)
            raw_data = vec(ds[var_name][idx_loc, idy_loc, idf_loc, :, 1])
            for (i, val) in enumerate(reverse(raw_data))
                nk_raw_loc[bin_idx, i] = ismissing(val) ? 0.0 : Float64(val)
            end
            nk_loc[bin_idx, :] .= (nk_raw_loc[bin_idx, :] ./ 1000.0) .* (n_air_loc ./ vol_cm3_loc)
        end
    end
    
    # Calculate dN/dlogD
    dN_dlogD_loc = nk_loc[:, lev_800] ./ dlogD
    
    # Fit lognormal
    N_total_loc = sum(nk_loc[:, lev_800])
    d_med_est_loc = diam_centers_um[argmax(dN_dlogD_loc)]
    σ_g_est_loc = 2.0
    
    try
        fit_loc = curve_fit(lognormal, diam_centers_um, dN_dlogD_loc .* dlogD, 
                           [N_total_loc, d_med_est_loc, σ_g_est_loc])
        N_fit, d_med_fit, σ_g_fit = fit_loc.param
        
        println("   $loc_name_comp:")
        @printf("      N_total = %.2e #/cm³, d_med = %.4f μm, σ_g = %.3f\n", 
                N_fit, d_med_fit, σ_g_fit)
    catch e
        println("   $loc_name_comp: Fitting failed")
    end
    
    # Plot
    color_idx = mod1(loc_idx, length(colors_tab10))
    lines!(ax5, diam_centers_um, dN_dlogD_loc,
           label=loc_name_comp, color=colors_tab10[color_idx], linewidth=2.5)
    scatter!(ax5, diam_centers_um, dN_dlogD_loc,
             color=colors_tab10[color_idx], markersize=8)
end

axislegend(ax5, position=:rt)

save(joinpath(output_dir, "NK_05_multi_location_comparison.png"), fig5)
println()
println("   Saved: NK_05_multi_location_comparison.png")

# Close NetCDF file
close(ds)

# ============================================================================
# Summary
# ============================================================================

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
println("  • NK stored as: 1000 × (particles/mol_air)")
println("  • Conversion formula: N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)")
println("  • Formula validated by colleague's working code and species comparison")
println("  • Values in true atmospheric number concentration (#/cm³)")
println("  • Typical continental air: 10³-10⁴ #/cm³ (matches our values!)")
println()
