#!/usr/bin/env julia
"""
Test the colleague's formula with our GEOSChem file

Their formula:
    Nk[idk] = num * air / vol / 1e3
where:
    num = NK value from file [claimed to be 1e3 * num/mol dry air]
    air = Met_AD / 28.9644e-3 [mol dry air]
    vol = Met_AIRVOL * 1e6 [cm³]
"""

using NCDatasets
using Printf
using Statistics

println("="^70)
println("Testing Colleague's Formula vs Our Empirical Formula")
println("="^70)
println()

# Open the file
ncfile = "GEOSChem.Custom.20190702_0000z.nc4"
ds = NCDataset(ncfile, "r")

# Select a location (Central USA from our analysis)
idx, idy, idf = 3, 13, 5  # Julia 1-indexed
idl = 72  # Surface layer

println("Location: Central USA")
println("  Grid indices: idx=$idx, idy=$idy, face=$idf")
println("  Vertical level: $idl (surface)")
println()

# Define mass bin edges (from colleague's code)
Mo = 1e-21 * 4^-3
Xk = zeros(16)
for k in 1:16
    if k < 15
        Xk[k] = Mo * 4^(k-1)
    else
        Xk[k] = Xk[k-1] * 32
    end
end

# Convert to size bin edges (assuming ρ = 1400 kg/m³)
ρ = 1400
Dp = 1e6 .* (6 .* Xk ./ ρ ./ π) .^ (1/3)  # μm
Dpmid = [(Dp[i] + Dp[i+1]) / 2 for i in 1:15]

println("TOMAS-15 size bins:")
for i in 1:5
    @printf("  Bin %2d: %.4f - %.4f μm (center: %.4f μm)\n", 
            i, Dp[i], Dp[i+1], Dpmid[i])
end
println("  ...")
println()

# Get meteorological variables
Met_AD = ds["Met_AD"][idx, idy, idf, idl, 1]
Met_AIRVOL = ds["Met_AIRVOL"][idx, idy, idf, idl, 1]
pressure = ds["Met_PMID"][idx, idy, idf, idl, 1]
temp = ds["Met_T"][idx, idy, idf, idl, 1]

println("Meteorological variables from file:")
@printf("  Met_AD (air mass):      %.6e kg\n", Met_AD)
@printf("  Met_AIRVOL (air volume): %.6e m³\n", Met_AIRVOL)
@printf("  Pressure:                %.2f Pa (%.2f hPa)\n", pressure, pressure/100)
@printf("  Temperature:             %.2f K\n", temp)
println()

# Calculate air density
ρ_air = Met_AD / Met_AIRVOL  # kg/m³
ρ_air_cm3 = ρ_air * 1e-6     # kg/cm³

@printf("Derived air density:\n")
@printf("  ρ_air = %.6f kg/m³\n", ρ_air)
@printf("  ρ_air = %.6e kg/cm³\n", ρ_air_cm3)
println()

# Read NK data
println("Reading NK data...")
Nk_colleague = zeros(15)  # Using colleague's formula
Nk_empirical = zeros(15)  # Using our empirical formula

M_air = 28.9644e-3  # kg/mol
air = Met_AD / M_air  # mol dry air
vol = Met_AIRVOL * 1e6  # cm³

@printf("Intermediate values:\n")
@printf("  air (mol):  %.6e mol\n", air)
@printf("  vol (cm³):  %.6e cm³\n", vol)
@printf("  air/vol:    %.6e mol/cm³\n", air/vol)
println()

for idk in 1:15
    var_name = "SpeciesConcVV_NK$(lpad(string(idk), 2, "0"))"
    num = ds[var_name][idx, idy, idf, idl, 1]  # NK value from file
    
    # Colleague's formula
    Nk_colleague[idk] = num * air / vol / 1e3
    
    # Our empirical formula
    Nk_empirical[idk] = num * ρ_air_cm3
end

println("="^70)
println("RESULTS: Number concentration per bin [#/cm³]")
println("="^70)
println()

@printf("%-6s %-12s %-15s %-15s %-10s\n", 
        "Bin", "Dp (μm)", "Colleague", "Empirical", "Ratio")
println("-"^70)

for idk in 1:15
    ratio = Nk_empirical[idk] / Nk_colleague[idk]
    @printf("%-6d %-12.4f %-15.1f %-15.1f %-10.2f\n", 
            idk, Dpmid[idk], Nk_colleague[idk], Nk_empirical[idk], ratio)
end

println("-"^70)
@printf("%-18s %-15.1f %-15.1f %-10.2f\n", 
        "TOTAL:", sum(Nk_colleague), sum(Nk_empirical), 
        sum(Nk_empirical) / sum(Nk_colleague))
println()

# Calculate dN/dlogDp for plotting comparison
dNdlogDp_colleague = zeros(15)
dNdlogDp_empirical = zeros(15)

for idk in 1:15
    dlogDp = log10(Dp[idk+1]) - log10(Dp[idk])
    dNdlogDp_colleague[idk] = Nk_colleague[idk] / dlogDp
    dNdlogDp_empirical[idk] = Nk_empirical[idk] / dlogDp
end

println("="^70)
println("dN/dlog(Dp) for plotting [cm⁻³]")
println("="^70)
println()

@printf("%-6s %-12s %-15s %-15s\n", 
        "Bin", "Dp (μm)", "Colleague", "Empirical")
println("-"^70)

for idk in 1:15
    @printf("%-6d %-12.4f %-15.1f %-15.1f\n", 
            idk, Dpmid[idk], dNdlogDp_colleague[idk], dNdlogDp_empirical[idk])
end

println()
@printf("Peak values:\n")
@printf("  Colleague:  %.1f cm⁻³ (bin %d)\n", 
        maximum(dNdlogDp_colleague), argmax(dNdlogDp_colleague))
@printf("  Empirical:  %.1f cm⁻³ (bin %d)\n", 
        maximum(dNdlogDp_empirical), argmax(dNdlogDp_empirical))
println()

# Compare with expected values
println("="^70)
println("Physical Validation")
println("="^70)
println()

expected_range = (3000, 50000)  # #/cm³ for continental surface
N_colleague = sum(Nk_colleague)
N_empirical = sum(Nk_empirical)

println("Expected range for continental surface: $(expected_range) #/cm³")
println()

@printf("Colleague formula total: %.1f #/cm³\n", N_colleague)
if expected_range[1] <= N_colleague <= expected_range[2]
    println("  ✓ In expected range")
else
    println("  ✗ Outside expected range")
end
println()

@printf("Empirical formula total: %.1f #/cm³\n", N_empirical)
if expected_range[1] <= N_empirical <= expected_range[2]
    println("  ✓ In expected range")
else
    println("  ✗ Outside expected range")
end
println()

# The colleague's plot shows peak dN/dlogDp ~ 6000-6500 cm⁻³
println("From colleague's plot: peak dN/dlog(Dp) ≈ 6000-6500 cm⁻³")
@printf("Our colleague formula: %.1f cm⁻³\n", maximum(dNdlogDp_colleague))
@printf("Our empirical formula: %.1f cm⁻³\n", maximum(dNdlogDp_empirical))
println()

println("="^70)
println("CONCLUSION:")
println("="^70)
println()

if abs(maximum(dNdlogDp_empirical) - 6250) < abs(maximum(dNdlogDp_colleague) - 6250)
    println("Our EMPIRICAL formula matches the colleague's plot better!")
    println("This confirms: NK is in #/kg, use N = NK × ρ_air")
else
    println("The COLLEAGUE formula matches their plot better.")
    println("This suggests: NK interpretation needs revision.")
end

close(ds)
