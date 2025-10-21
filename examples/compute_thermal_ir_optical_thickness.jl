#!/usr/bin/env julia
"""
Compute Thermal IR Optical Thickness for a 10m Atmospheric Layer

This example demonstrates how to use vSmartMOM to compute the optical thickness
of a thin atmospheric layer (10m) near the ground with tropical water vapor
abundances in the thermal infrared.

The optical thickness τ is computed from the absorption cross-section σ:
    τ = σ × n × Δz

where:
    - σ = absorption cross-section [cm²/molecule]
    - n = number density [molecules/cm³]
    - Δz = layer thickness [cm]
"""

using vSmartMOM
using vSmartMOM.Absorption
using Plots
using Printf

println("="^70)
println("Thermal IR Optical Thickness Calculation")
println("="^70)
println()

# =============================================================================
# STEP 1: Define atmospheric conditions (tropical near-surface)
# =============================================================================

println("Setting up atmospheric parameters...")
println()

# Layer properties
Δz_m = 10.0                          # Layer thickness [m]
Δz_cm = Δz_m * 100.0                 # Layer thickness [cm]

# Tropical surface conditions
T = 300.0                            # Temperature [K]
P = 1013.25                          # Pressure [hPa] (sea level)

# Tropical water vapor mixing ratio
# Typical tropical boundary layer: ~20 g/kg = 0.02 kg H2O / kg dry air
# Convert to volume mixing ratio (VMR):
# VMR = mass mixing ratio × (M_air / M_H2O)
mass_mix_ratio = 0.02                # kg H2O / kg dry air
M_air = 28.97                        # g/mol (dry air)
M_H2O = 18.015                       # g/mol (water)
vmr_H2O = mass_mix_ratio * (M_air / M_H2O)  # dimensionless

@printf("Atmospheric conditions:\n")
@printf("  Temperature:    %.1f K\n", T)
@printf("  Pressure:       %.2f hPa\n", P)
@printf("  Layer thickness: %.1f m (%.0f cm)\n", Δz_m, Δz_cm)
@printf("  H2O mass mixing ratio: %.3f kg/kg\n", mass_mix_ratio)
@printf("  H2O VMR:        %.4f (%.2f%%)\n", vmr_H2O, vmr_H2O * 100)
println()

# Calculate number density from ideal gas law
# n = P / (k_B × T)
# where P is in Pa, k_B = 1.380649e-23 J/K
k_B = 1.380649e-23                   # Boltzmann constant [J/K]
P_Pa = P * 100.0                     # Convert hPa to Pa
n_total = P_Pa / (k_B * T)           # Total number density [molecules/m³]
n_total_cm3 = n_total * 1e-6         # Convert to [molecules/cm³]
n_H2O = n_total_cm3 * vmr_H2O        # H2O number density [molecules/cm³]

@printf("Number densities:\n")
@printf("  Total air: %.3e molecules/cm³\n", n_total_cm3)
@printf("  H2O:       %.3e molecules/cm³\n", n_H2O)
println()

# =============================================================================
# STEP 2: Download/read HITRAN data for H2O in thermal IR
# =============================================================================

println("Loading HITRAN data for H2O in thermal infrared...")
println()

# Thermal IR atmospheric window: 8-12 μm
# Convert to wavenumber: ν [cm⁻¹] = 10000 / λ [μm]
λ_min = 8.0                          # μm
λ_max = 12.0                         # μm
ν_max = 10000.0 / λ_min              # cm⁻¹ (1250)
ν_min = 10000.0 / λ_max              # cm⁻¹ (833)

@printf("Spectral range:\n")
@printf("  Wavelength: %.1f - %.1f μm\n", λ_min, λ_max)
@printf("  Wavenumber: %.1f - %.1f cm⁻¹\n", ν_min, ν_max)
println()

# Read HITRAN data for H2O (molecule 1, isotopologue 1)
# This will download from the HITRAN database if not already cached
h2o_data = read_hitran(artifact("H2O"), mol=1, iso=1, ν_min=ν_min, ν_max=ν_max)

@printf("HITRAN data loaded:\n")
@printf("  Number of transitions: %d\n", length(h2o_data.νᵢ))
@printf("  Wavenumber range: %.2f - %.2f cm⁻¹\n", 
        minimum(h2o_data.νᵢ), maximum(h2o_data.νᵢ))
println()

# =============================================================================
# STEP 3: Create absorption model
# =============================================================================

println("Creating Voigt line-shape model...")
println()

# Use Voigt profile (includes both Doppler and pressure broadening)
# Important for thermal IR where both effects are significant
h2o_model = make_hitran_model(
    h2o_data, 
    Voigt(),
    wing_cutoff=25.0,                # cm⁻¹ (line wing cutoff)
    vmr=vmr_H2O,                     # Volume mixing ratio
    architecture=CPU()
)

println("✓ Model created successfully")
println()

# =============================================================================
# STEP 4: Compute absorption cross-section
# =============================================================================

println("Computing absorption cross-sections...")
println()

# Define high-resolution wavenumber grid
ν_grid = ν_min:0.01:ν_max            # cm⁻¹ (0.01 cm⁻¹ resolution)

# Compute cross-section at specified P, T
# Returns σ in units of [cm²/molecule]
σ_H2O = absorption_cross_section(h2o_model, ν_grid, P, T)

@printf("Cross-section statistics:\n")
@printf("  Grid points:    %d\n", length(ν_grid))
@printf("  Mean σ:         %.3e cm²/molecule\n", mean(σ_H2O))
@printf("  Max σ:          %.3e cm²/molecule\n", maximum(σ_H2O))
@printf("  Min σ:          %.3e cm²/molecule\n", minimum(σ_H2O[σ_H2O .> 0]))
println()

# =============================================================================
# STEP 5: Calculate optical thickness
# =============================================================================

println("Computing optical thickness...")
println()

# Optical thickness: τ = σ × n × Δz
# where:
#   σ [cm²/molecule] × n [molecules/cm³] × Δz [cm] = τ [dimensionless]
τ_H2O = σ_H2O .* n_H2O .* Δz_cm

@printf("Optical thickness statistics:\n")
@printf("  Mean τ:         %.4f\n", mean(τ_H2O))
@printf("  Max τ:          %.4f\n", maximum(τ_H2O))
@printf("  Min τ:          %.4e\n", minimum(τ_H2O[τ_H2O .> 0]))
@printf("  Median τ:       %.4f\n", median(τ_H2O))
println()

# Calculate transmittance: T = exp(-τ)
transmittance = exp.(-τ_H2O)

@printf("Transmittance statistics:\n")
@printf("  Mean T:         %.4f (%.2f%%)\n", mean(transmittance), mean(transmittance)*100)
@printf("  Min T:          %.4f (%.2f%%)\n", minimum(transmittance), minimum(transmittance)*100)
println()

# =============================================================================
# STEP 6: Visualize results
# =============================================================================

println("Creating visualizations...")
println()

# Convert wavenumber back to wavelength for plotting
λ_grid = 10000.0 ./ ν_grid           # μm

# Plot 1: Absorption cross-section vs wavelength
p1 = plot(λ_grid, σ_H2O,
          xlabel="Wavelength [μm]",
          ylabel="Cross-section [cm²/molecule]",
          title="H₂O Absorption Cross-Section (T=$T K, P=$P hPa)",
          yaxis=:log,
          label="",
          linewidth=1.5,
          color=:blue,
          size=(900, 500),
          margin=5Plots.mm)

# Plot 2: Optical thickness vs wavelength
p2 = plot(λ_grid, τ_H2O,
          xlabel="Wavelength [μm]",
          ylabel="Optical Thickness τ",
          title="Optical Thickness for $(Δz_m)m Layer (Tropical Conditions)",
          yaxis=:log,
          label="",
          linewidth=1.5,
          color=:red,
          size=(900, 500),
          margin=5Plots.mm)

# Add reference line for τ=1
hline!(p2, [1.0], linestyle=:dash, color=:black, label="τ = 1", linewidth=2)

# Plot 3: Transmittance vs wavelength
p3 = plot(λ_grid, transmittance,
          xlabel="Wavelength [μm]",
          ylabel="Transmittance",
          title="Layer Transmittance exp(-τ) for $(Δz_m)m",
          ylims=(0, 1),
          label="",
          linewidth=1.5,
          color=:green,
          size=(900, 500),
          margin=5Plots.mm)

# Combine all plots
p_combined = plot(p1, p2, p3, layout=(3, 1), size=(900, 1200))

# Save plots
output_dir = "examples/output"
mkpath(output_dir)

savefig(p_combined, joinpath(output_dir, "thermal_ir_optical_thickness.png"))
println("✓ Plots saved to: $output_dir/thermal_ir_optical_thickness.png")
println()

# =============================================================================
# STEP 7: Summary and key wavelengths
# =============================================================================

println("="^70)
println("Summary of Key Results")
println("="^70)
println()

# Find wavelengths of maximum absorption
idx_max = argmax(τ_H2O)
@printf("Maximum optical thickness:\n")
@printf("  τ_max = %.4f at λ = %.3f μm (ν = %.2f cm⁻¹)\n", 
        τ_H2O[idx_max], λ_grid[idx_max], ν_grid[idx_max])
println()

# Find atmospheric window regions (where τ < 0.1)
window_mask = τ_H2O .< 0.1
if any(window_mask)
    window_λ = λ_grid[window_mask]
    @printf("Atmospheric window regions (τ < 0.1):\n")
    @printf("  Wavelength ranges: %.3f - %.3f μm\n", 
            minimum(window_λ), maximum(window_λ))
    @printf("  Average transmittance: %.2f%%\n", 
            mean(transmittance[window_mask]) * 100)
else
    println("No clear atmospheric windows in this spectral range")
end
println()

# Integration over spectral range
mean_τ = mean(τ_H2O)
@printf("Spectrally-averaged quantities:\n")
@printf("  <τ> = %.4f\n", mean_τ)
@printf("  <T> = %.4f (%.2f%%)\n", mean(transmittance), mean(transmittance)*100)
println()

println("="^70)
println("✓ Calculation complete!")
println("="^70)
println()
println("Physical interpretation:")
println("  • Optical thickness τ > 1: Layer is optically thick (significant absorption)")
println("  • Optical thickness τ < 1: Layer is optically thin (weak absorption)")
println("  • Transmittance = exp(-τ): Fraction of radiation transmitted through layer")
println()
println("For a $(Δz_m)m tropical layer:")
@printf("  • Average transmittance: %.2f%%\n", mean(transmittance)*100)
@printf("  • Strong absorption lines have τ > 1 (%.1f%% of spectrum)\n", 
        sum(τ_H2O .> 1.0) / length(τ_H2O) * 100)
println()
