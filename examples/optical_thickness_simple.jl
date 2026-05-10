#!/usr/bin/env julia
"""
Simple Example: Optical Thickness Computation

Quick demonstration of computing optical thickness for a 10m layer
with tropical water vapor in the thermal IR.
"""

using vSmartMOM
using vSmartMOM.Absorption
using Printf, Statistics

println("\n" * "="^70)
println("Simple Optical Thickness Calculation")
println("="^70 * "\n")

# =============================================================================
# Step 1: Define atmospheric layer
# =============================================================================

T = 300.0          # Temperature [K]
P = 1013.25        # Pressure [hPa]
Δz = 30.0          # Layer thickness [m]

# Water vapor: tropical conditions (2% mass mixing ratio)
vmr_H2O = 0.032    # Volume mixing ratio (~2% by mass)

println("Layer properties:")
@printf("  T = %.1f K, P = %.1f hPa, Δz = %.1f m\n", T, P, Δz)
@printf("  H₂O VMR = %.3f (%.1f%%)\n\n", vmr_H2O, vmr_H2O*100)

# =============================================================================
# Step 2: Calculate number density
# =============================================================================

# From ideal gas law: n = P/(k_B × T)
k_B = 1.380649e-23                    # Boltzmann constant [J/K]
n_air = (P * 100) / (k_B * T) * 1e-6  # Total [molecules/cm³]
n_H2O = n_air * vmr_H2O               # H₂O [molecules/cm³]

@printf("Number densities:\n")
@printf("  n_air = %.3e molecules/cm³\n", n_air)
@printf("  n_H2O = %.3e molecules/cm³\n\n", n_H2O)

# =============================================================================
# Step 3: Load HITRAN and compute cross-section
# =============================================================================

println("Loading HITRAN data for thermal IR (8-12 μm)...")

# Thermal IR window: 8-12 μm = 833-1250 cm⁻¹
h2o_data = read_hitran(artifact("H2O"), mol=1, iso=1, ν_min=250.0, ν_max=3350.0)
@printf("  Loaded %d transitions\n\n", length(h2o_data.νᵢ))

println("Creating absorption model (Voigt profile)...")
h2o_model = make_hitran_model(h2o_data, Voigt(), vmr=vmr_H2O, architecture=CPU())
println("  ✓ Model created\n")

println("Computing absorption cross-section...")
ν_grid = 250.0:0.1:3350.0                # Wavenumber grid [cm⁻¹]
σ = absorption_cross_section(h2o_model, ν_grid, P, T)  # [cm²/molecule]
@printf("  ✓ Computed %d points\n", length(σ))
@printf("  Max σ = %.3e cm²/molecule\n\n", maximum(σ))

# =============================================================================
# Step 4: Compute optical thickness
# =============================================================================

println("Computing optical thickness...")
Δz_cm = Δz * 100                     # Convert to [cm]
τ = σ .* n_H2O .* Δz_cm              # Optical thickness [dimensionless]

println("\n" * "="^70)
println("RESULTS:")
println("="^70)
@printf("\nOptical thickness τ for %.1f m layer:\n", Δz)
@printf("  Mean:   %.4f\n", mean(τ))
@printf("  Median: %.4f\n", median(τ))
@printf("  Max:    %.4f\n", maximum(τ))
@printf("  Min:    %.4e\n", minimum(τ[τ .> 0]))

# Transmittance
T_trans = exp.(-τ)
@printf("\nTransmittance exp(-τ):\n")
@printf("  Mean:   %.4f (%.1f%% of radiation transmitted)\n", 
        mean(T_trans), mean(T_trans)*100)
@printf("  Min:    %.4f (%.1f%% at strongest absorption)\n", 
        minimum(T_trans), minimum(T_trans)*100)

# Summary
println("\n" * "="^70)
println("INTERPRETATION:")
println("="^70)
println("\nFor a 10m tropical boundary layer:")
println("  • Optical thickness varies strongly with wavelength")
println("  • Mean τ ≈ $(round(mean(τ), digits=3))")
if mean(τ) < 1.0
    println("  • On average, the layer is OPTICALLY THIN (τ < 1)")
    println("  • Most radiation passes through with minimal absorption")
else
    println("  • On average, the layer is OPTICALLY THICK (τ > 1)")
    println("  • Significant absorption occurs in this layer")
end

# Count spectral regions
frac_thick = sum(τ .> 1.0) / length(τ)
@printf("\n  • %.1f%% of spectrum has τ > 1 (strong absorption lines)\n", frac_thick*100)
@printf("  • %.1f%% of spectrum has τ < 0.1 (atmospheric windows)\n", 
        sum(τ .< 0.1) / length(τ) * 100)

println("\n" * "="^70)
println("✓ Calculation complete!")
println("="^70 * "\n")

# =============================================================================
# Formula used
# =============================================================================

println("Formula:")
println("  τ = σ × n × Δz")
println("  where:")
println("    σ = absorption cross-section [cm²/molecule]")
println("    n = number density [molecules/cm³]")
println("    Δz = layer thickness [cm]")
println("")
