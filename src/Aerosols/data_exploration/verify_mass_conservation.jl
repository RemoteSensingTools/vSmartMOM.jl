"""
Verify that the bimodal fit conserves mass per bin.

This script checks that for each TOMAS bin, the integral of the fitted
lognormal distribution over that bin equals the actual bin value.
"""

using NCDatasets, Statistics, Distributions, Printf

# Read the fitted parameters from the previous run
# These come from the boundary layer (897 hPa) fit
N1, d_med1, σ1 = 2.32e+03, 0.2000, 1.499     # Mode 1 (Aitken)
N2, d_med2, σ2 = 8.68e+02, 0.9670, 1.012     # Mode 2 (Accumulation)
r_med1, r_med2 = d_med1/2, d_med2/2  # Convert diameter to radius (in μm)

println("="^70)
println("Verifying Mass Conservation for Bimodal Fit")
println("="^70)
println()
println("Fit Parameters:")
println("  Mode 1 (Aitken):       N = $(N1) #/cm³, d_med = $(d_med1) μm, σ_g = $(σ1)")
println("  Mode 2 (Accumulation): N = $(N2) #/cm³, d_med = $(d_med2) μm, σ_g = $(σ2)")
println()

# Function to integrate lognormal over a bin
function integrate_lognormal_bin(r_left, r_right, N_total, r_med, σ_g)
    dist = Normal(0, 1)
    arg_right = log(r_right/r_med) / (sqrt(2) * log(σ_g))
    arg_left = log(r_left/r_med) / (sqrt(2) * log(σ_g))
    return N_total * (cdf(dist, arg_right) - cdf(dist, arg_left))
end

# Read actual data - should get EXACTLY the same data as explore script!
ds = NCDataset("GEOSChem.Custom.20190702_0000z.nc4")

# Central USA location - EXACTLY as in explore_NK_cairomakie.jl
# locations[1] = (3, 13, 5, "Central USA")  where order is (idx, idy, face, name)
idx, idy, idf = 3, 13, 5

println("Location: Central USA (idx=$idx, idy=$idy, face=$idf)")
println()

# Read meteorology EXACTLY as in explore_NK_cairomakie.jl
# Lines 253-257 of that script
n_lev = ds.dim["lev"]
dp = vec(ds["Met_DELP"][idx, idy, idf, :, 1])
sp = Float64(ds["Met_PS2WET"][idx, idy, idf, 1])

# Pressure at half-levels (flip from BOA→TOA to TOA→BOA)
pressure_half = vcat([sp], sp .- cumsum(dp))
pressure_half_toa2boa = reverse(pressure_half)
pressure_centers = (pressure_half_toa2boa[1:end-1] .+ pressure_half_toa2boa[2:end]) ./ 2

# Read meteorology for NK conversion - EXACTLY as in explore script lines 289-291
Met_AD = reverse(vec(ds["Met_AD"][idx, idy, idf, :, 1]))
Met_AIRVOL = reverse(vec(ds["Met_AIRVOL"][idx, idy, idf, :, 1]))

# Find boundary layer level (~900 hPa) - EXACTLY as in explore_NK_cairomakie.jl
lev_bl = argmin(abs.(pressure_centers .- 900))
p_bl = pressure_centers[lev_bl]

println("Checking boundary layer level: $(Int(round(p_bl))) hPa (level index $lev_bl)")
println()

# TOMAS15 bin structure (15 bins, factor of 2 spacing)
# Bins cover 0.01 to 10 μm diameter
r_min_um = 0.01 / 2  # Smallest radius (μm)
r_max_um = 10.0 / 2  # Largest radius (μm)
n_bins = 15

# Calculate bin edges (radius in μm)
# Each bin spans factor of 2^(1/3) in radius
log_r_min = log10(r_min_um)
log_r_max = log10(r_max_um)
log_r_centers = range(log_r_min, log_r_max, length=n_bins)
r_centers = 10 .^ log_r_centers

# Bin edges: each bin is 2^(1/3) wide in radius
bin_factor = 2.0^(1.0/3.0)
r_edges = r_centers ./ sqrt(bin_factor)
push!(r_edges, r_centers[end] * sqrt(bin_factor))

# Read NK data and convert to N (#/cm³) for each bin - EXACTLY as in explore script
M_air = 28.9644e-3  # kg/mol
N_bins_actual = zeros(n_bins)

for i in 1:n_bins
    nk_var = @sprintf("SpeciesConcVV_NK%02d", i)
    # Must read and reverse to match explore script
    NK_profile = reverse(vec(ds[nk_var][idx, idy, idf, :, 1]))
    NK = NK_profile[lev_bl]
    N_bins_actual[i] = (NK/1000) * (Met_AD[lev_bl]/M_air) / (Met_AIRVOL[lev_bl]*1e6)
end

close(ds)


# Calculate integrated values from fit
N_bins_fit = zeros(n_bins)
N_bins_mode1 = zeros(n_bins)
N_bins_mode2 = zeros(n_bins)

for i in 1:n_bins
    N_bins_mode1[i] = integrate_lognormal_bin(r_edges[i], r_edges[i+1], N1, r_med1, σ1)
    N_bins_mode2[i] = integrate_lognormal_bin(r_edges[i], r_edges[i+1], N2, r_med2, σ2)
    N_bins_fit[i] = N_bins_mode1[i] + N_bins_mode2[i]
end

# Print comparison table
println("-"^70)
println(@sprintf("%4s %8s %10s %10s %10s %8s", 
                 "Bin", "r_ctr", "N_actual", "N_fit", "Diff", "Error%"))
println(@sprintf("%4s %8s %10s %10s %10s %8s",
                 "", "(μm)", "(#/cm³)", "(#/cm³)", "(#/cm³)", ""))
println("-"^70)

global total_actual = 0.0
global total_fit = 0.0
global max_error = 0.0

for i in 1:n_bins
    diff = N_bins_fit[i] - N_bins_actual[i]
    if N_bins_actual[i] > 0
        error_pct = 100 * abs(diff) / N_bins_actual[i]
        global max_error = max(max_error, error_pct)
    else
        error_pct = 0.0
    end
    
    global total_actual += N_bins_actual[i]
    global total_fit += N_bins_fit[i]
    
    println(@sprintf("%4d %8.4f %10.2f %10.2f %10.2f %7.1f%%",
                     i, r_centers[i], N_bins_actual[i], N_bins_fit[i], 
                     diff, error_pct))
end

println("-"^70)
println(@sprintf("%-4s %8s %10.2f %10.2f %10.2f %7.1f%%",
                 "TOT", "", total_actual, total_fit, 
                 total_fit - total_actual,
                 100 * abs(total_fit - total_actual) / total_actual))
println("-"^70)

println()
println("Summary:")
println("  Total N (actual):  $(round(total_actual, digits=1)) #/cm³")
println("  Total N (fit):     $(round(total_fit, digits=1)) #/cm³")
println("  Overall error:     $(round(100*abs(total_fit-total_actual)/total_actual, digits=2))%")
println("  Max bin error:     $(round(max_error, digits=1))%")
println()

# Check if mass is reasonably conserved (within 5% per bin, 1% total)
if max_error < 5.0 && abs(total_fit - total_actual) / total_actual < 0.01
    println("✓ Mass conservation: EXCELLENT")
elseif max_error < 10.0 && abs(total_fit - total_actual) / total_actual < 0.05
    println("✓ Mass conservation: GOOD")
else
    println("⚠ Mass conservation: NEEDS IMPROVEMENT")
end
println()
