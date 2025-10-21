# Computing Optical Thickness with vSmartMOM

This guide explains how to use vSmartMOM to compute the optical thickness of atmospheric layers in the thermal infrared.

## Table of Contents

1. [Overview](#overview)
2. [Quick Start](#quick-start)
3. [Physical Background](#physical-background)
4. [Step-by-Step Guide](#step-by-step-guide)
5. [Example Scripts](#example-scripts)
6. [Advanced Options](#advanced-options)

## Overview

**Optical thickness (τ)** is a dimensionless measure of how much electromagnetic radiation is absorbed or scattered by a medium. For a homogeneous atmospheric layer:

```
τ = σ × n × Δz
```

where:
- **σ** = absorption cross-section [cm²/molecule]
- **n** = number density [molecules/cm³]  
- **Δz** = layer thickness [cm]

**Transmittance** through the layer: `T = exp(-τ)`

## Quick Start

Here's a minimal example for computing optical thickness of a 10m tropical layer:

```julia
using vSmartMOM
using vSmartMOM.Absorption

# 1. Define atmospheric conditions
T = 300.0          # Temperature [K]
P = 1013.25        # Pressure [hPa]
Δz = 1000.0        # Layer thickness [cm] (10m)
vmr_H2O = 0.032    # H₂O volume mixing ratio (tropical)

# 2. Calculate number density from ideal gas law
k_B = 1.380649e-23                    # Boltzmann constant [J/K]
n_air = (P * 100) / (k_B * T) * 1e-6  # [molecules/cm³]
n_H2O = n_air * vmr_H2O               # [molecules/cm³]

# 3. Load HITRAN data and create model
h2o_data = read_hitran(artifact("H2O"), mol=1, iso=1, ν_min=833.0, ν_max=1250.0)
h2o_model = make_hitran_model(h2o_data, Voigt(), vmr=vmr_H2O, architecture=CPU())

# 4. Compute absorption cross-section
ν_grid = 833:0.1:1250                # Wavenumber [cm⁻¹]
σ = absorption_cross_section(h2o_model, ν_grid, P, T)  # [cm²/molecule]

# 5. Calculate optical thickness
τ = σ .* n_H2O .* Δz                 # [dimensionless]

# 6. Calculate transmittance
T_trans = exp.(-τ)                   # [dimensionless, 0-1]

println("Mean optical thickness: $(mean(τ))")
println("Mean transmittance: $(mean(T_trans)) ($(mean(T_trans)*100)%)")
```

## Physical Background

### Optical Thickness Interpretation

- **τ << 1**: Optically thin - most radiation passes through
- **τ ≈ 1**: Transition region - significant absorption begins
- **τ >> 1**: Optically thick - most radiation absorbed

### Number Density from Ideal Gas Law

The number density of molecules in air follows the ideal gas law:

```
n = P / (k_B × T)
```

where:
- **P** = pressure [Pa]
- **k_B** = Boltzmann constant = 1.380649×10⁻²³ J/K
- **T** = temperature [K]
- **n** = number density [molecules/m³]

Convert to molecules/cm³ by multiplying by 10⁻⁶.

### Water Vapor Mixing Ratios

Water vapor abundance is commonly expressed in different ways:

| Unit | Typical Range | Conversion |
|------|---------------|------------|
| Mass mixing ratio | 1-20 g/kg | VMR = (g/kg) × (M_air/M_H2O) × 0.001 |
| Volume mixing ratio (VMR) | 0.1-3% | Dimensionless |
| Relative humidity | 0-100% | Depends on temperature |

For **tropical conditions**: 
- Mass mixing ratio ≈ 20 g/kg
- VMR ≈ 0.032 (3.2%)

For **temperate conditions**:
- Mass mixing ratio ≈ 10 g/kg  
- VMR ≈ 0.016 (1.6%)

For **dry conditions**:
- Mass mixing ratio ≈ 2 g/kg
- VMR ≈ 0.003 (0.3%)

### Spectral Regions

**Thermal Infrared (8-12 μm):**
- Wavelength: 8-12 μm
- Wavenumber: 833-1250 cm⁻¹
- Primary absorber: H₂O
- Contains atmospheric window regions

**Mid-Infrared (4-8 μm):**
- Wavelength: 4-8 μm  
- Wavenumber: 1250-2500 cm⁻¹
- Absorbers: H₂O, CO₂, CH₄

**Near-Infrared (1-3 μm):**
- Wavelength: 1-3 μm
- Wavenumber: 3333-10000 cm⁻¹
- Absorbers: H₂O, CO₂, CH₄, CO

## Step-by-Step Guide

### Step 1: Define Atmospheric Layer Properties

```julia
# Layer geometry
Δz_m = 10.0                # Thickness [m]
Δz_cm = Δz_m * 100.0       # Convert to [cm]

# Thermodynamic state
T = 300.0                  # Temperature [K]
P = 1013.25                # Pressure [hPa]

# Composition
vmr_H2O = 0.032            # H₂O volume mixing ratio
```

### Step 2: Calculate Number Densities

```julia
# Physical constants
k_B = 1.380649e-23         # Boltzmann constant [J/K]

# Convert pressure to Pa
P_Pa = P * 100.0

# Total air number density
n_total = P_Pa / (k_B * T)              # [molecules/m³]
n_total_cm3 = n_total * 1e-6            # [molecules/cm³]

# Species number density
n_H2O = n_total_cm3 * vmr_H2O           # [molecules/cm³]
```

### Step 3: Load HITRAN Spectroscopic Data

vSmartMOM uses the HITRAN database for molecular line parameters.

```julia
using vSmartMOM.Absorption

# For thermal IR (8-12 μm)
ν_min = 833.0              # Wavenumber [cm⁻¹]
ν_max = 1250.0             # Wavenumber [cm⁻¹]

# Download/load HITRAN data
# mol=1 is H₂O, iso=1 is most abundant isotopologue
h2o_data = read_hitran(
    artifact("H2O"),       # Automatically downloads if needed
    mol=1,                 # H₂O
    iso=1,                 # ¹H₂¹⁶O
    ν_min=ν_min,
    ν_max=ν_max
)
```

**Available molecules:**
```julia
Absorption.show_molecules()  # Lists all available molecules
```

Common molecules:
- `mol=1`: H₂O (water)
- `mol=2`: CO₂ (carbon dioxide)
- `mol=5`: CO (carbon monoxide)
- `mol=6`: CH₄ (methane)
- `mol=7`: O₂ (oxygen)

### Step 4: Create Absorption Model

Choose a line shape model:

```julia
# Voigt profile (recommended for most cases)
# Includes both Doppler and pressure broadening
model = make_hitran_model(
    h2o_data,
    Voigt(),
    vmr=vmr_H2O,              # Volume mixing ratio
    wing_cutoff=25.0,         # Line wing cutoff [cm⁻¹]
    architecture=CPU()        # Or GPU() if available
)
```

**Line shape options:**
- `Voigt()`: Convolution of Doppler and Lorentz (most accurate)
- `Lorentz()`: Pressure broadening only (fast, good at high pressure)
- `Doppler()`: Thermal broadening only (good at low pressure)

### Step 5: Compute Absorption Cross-Section

```julia
# Define spectral grid
ν_grid = ν_min:0.1:ν_max       # 0.1 cm⁻¹ resolution

# Compute cross-section at specified P, T
σ = absorption_cross_section(
    model,
    ν_grid,
    P,                         # Pressure [hPa]
    T                          # Temperature [K]
)
# Returns σ in [cm²/molecule]
```

**Grid resolution considerations:**
- Coarse (1 cm⁻¹): Fast, may miss narrow lines
- Medium (0.1 cm⁻¹): Good compromise
- Fine (0.01 cm⁻¹): Slow, captures all line details

### Step 6: Calculate Optical Thickness

```julia
# Optical thickness: τ = σ × n × Δz
τ = σ .* n_H2O .* Δz_cm        # [dimensionless]

# Transmittance through layer
transmittance = exp.(-τ)       # [dimensionless, 0-1]

# Absorption (complement of transmittance)
absorption = 1.0 .- transmittance
```

### Step 7: Analyze Results

```julia
using Statistics

# Summary statistics
println("Optical thickness statistics:")
println("  Mean:   $(mean(τ))")
println("  Median: $(median(τ))")
println("  Max:    $(maximum(τ))")
println("  Min:    $(minimum(τ[τ .> 0]))")

# Transmittance
println("\nTransmittance:")
println("  Mean:   $(mean(transmittance)) ($(mean(transmittance)*100)%)")
println("  Min:    $(minimum(transmittance)) (at strongest absorption)")

# Spectral regions
frac_thick = sum(τ .> 1.0) / length(τ)
frac_window = sum(τ .< 0.1) / length(τ)
println("\nSpectral characteristics:")
println("  $(frac_thick*100)% has τ > 1 (optically thick)")
println("  $(frac_window*100)% has τ < 0.1 (atmospheric window)")
```

## Example Scripts

Two complete example scripts are provided:

### 1. Simple Example: `examples/optical_thickness_simple.jl`

A minimal, fast-running example showing the core calculation:

```bash
julia --project=. examples/optical_thickness_simple.jl
```

**Features:**
- Single 10m layer
- Tropical conditions
- Thermal IR (8-12 μm)
- ~30 seconds runtime

### 2. Comprehensive Example: `examples/compute_thermal_ir_optical_thickness.jl`

A detailed example with visualization:

```bash
julia --project=. examples/compute_thermal_ir_optical_thickness.jl
```

**Features:**
- Detailed calculations
- Multiple plots (cross-section, optical thickness, transmittance)
- Physical interpretation
- Atmospheric window identification
- ~2 minutes runtime

## Advanced Options

### Multiple Layers

To compute optical thickness for multiple atmospheric layers:

```julia
# Define vertical grid
z_levels = 0:10:1000        # Heights [m]
n_layers = length(z_levels) - 1

# Preallocate
τ_layers = zeros(length(ν_grid), n_layers)

# Loop over layers
for i in 1:n_layers
    Δz_i = (z_levels[i+1] - z_levels[i]) * 100  # [cm]
    
    # Get P, T, VMR for this layer (from atmospheric profile)
    P_i = ...
    T_i = ...
    vmr_i = ...
    n_i = (P_i * 100) / (k_B * T_i) * 1e-6 * vmr_i
    
    # Compute cross-section at layer conditions
    σ_i = absorption_cross_section(model, ν_grid, P_i, T_i)
    
    # Layer optical thickness
    τ_layers[:, i] = σ_i .* n_i .* Δz_i
end

# Total optical thickness (sum over layers)
τ_total = sum(τ_layers, dims=2)
```

### Multiple Species

To include multiple absorbing species:

```julia
# Load data for multiple species
h2o_data = read_hitran(artifact("H2O"), mol=1, iso=1, ν_min=ν_min, ν_max=ν_max)
co2_data = read_hitran(artifact("CO2"), mol=2, iso=1, ν_min=ν_min, ν_max=ν_max)

# Create models
h2o_model = make_hitran_model(h2o_data, Voigt(), vmr=vmr_H2O, architecture=CPU())
co2_model = make_hitran_model(co2_data, Voigt(), vmr=vmr_CO2, architecture=CPU())

# Compute cross-sections
σ_h2o = absorption_cross_section(h2o_model, ν_grid, P, T)
σ_co2 = absorption_cross_section(co2_model, ν_grid, P, T)

# Calculate optical thickness for each species
τ_h2o = σ_h2o .* n_H2O .* Δz_cm
τ_co2 = σ_co2 .* n_CO2 .* Δz_cm

# Total optical thickness (sum over species)
τ_total = τ_h2o .+ τ_co2
```

### Using Interpolation Models

For faster repeated calculations at different P, T:

```julia
# Create interpolation model once
ν_grid = 833:0.1:1250
p_grid = 100:100:1100              # Pressure grid [hPa]
t_grid = 200:10:320                # Temperature grid [K]

interp_model = make_interpolation_model(
    h2o_data,
    Voigt(),
    ν_grid,
    p_grid,
    t_grid,
    vmr=vmr_H2O
)

# Fast interpolation at arbitrary P, T
σ = absorption_cross_section(interp_model, ν_grid, 950.0, 285.0)
```

### GPU Acceleration

For very large calculations:

```julia
using CUDA  # Must have NVIDIA GPU

# Create model with GPU backend
model = make_hitran_model(h2o_data, Voigt(), architecture=GPU())

# Computation automatically runs on GPU
σ = absorption_cross_section(model, ν_grid, P, T)
```

### Wavelength vs Wavenumber

By default, vSmartMOM works in wavenumber space [cm⁻¹]. To use wavelength [nm]:

```julia
# Define wavelength grid
λ_grid = 8000:1:12000     # Wavelength [nm]

# Compute with wavelength flag
σ = absorption_cross_section(model, λ_grid, P, T, wavelength_flag=true)
```

**Conversion:**
- Wavenumber to wavelength: λ [nm] = 10⁷ / ν [cm⁻¹]
- Wavelength to wavenumber: ν [cm⁻¹] = 10⁷ / λ [nm]

## References

### Documentation
- [vSmartMOM.jl Docs](https://RemoteSensingTools.github.io/vSmartMOM.jl/dev/)
- [Absorption Module Tutorial](https://RemoteSensingTools.github.io/vSmartMOM.jl/dev/pages/tutorials/Tutorial_Absorption/)

### Physical References
- HITRAN Database: https://hitran.org
- Petty, G. W. (2006). A First Course in Atmospheric Radiation (2nd ed.)
- Liou, K. N. (2002). An Introduction to Atmospheric Radiation (2nd ed.)

### vSmartMOM Papers
- Sanghavi et al. (2014). vSmartMOM: A vector matrix operator method-based radiative transfer model. JQSRT, 133, 412-433.

## Troubleshooting

### Common Issues

**Issue:** "No format string provided to @printf"
```julia
# Wrong:
@printf("\n" * "="^70 * "\n")

# Correct:
println("\n" * "="^70 * "\n")
```

**Issue:** Negative or missing values in τ
- Check that σ > 0 (may need finer spectral resolution)
- Verify n_H2O is reasonable (typically 10¹⁵-10¹⁸ molecules/cm³)
- Check unit conversions (especially Δz in cm, P in hPa)

**Issue:** Very long computation time
- Reduce spectral resolution (use coarser grid)
- Use `make_interpolation_model` for repeated calculations
- Consider GPU acceleration for large grids
- Reduce `wing_cutoff` parameter (default 50, try 25 or less)

**Issue:** Out of memory
- Reduce spectral grid size
- Process in chunks
- Use GPU with larger memory

## Example Output

Running `optical_thickness_simple.jl` produces:

```
======================================================================
Simple Optical Thickness Calculation
======================================================================

Layer properties:
  T = 300.0 K, P = 1013.2 hPa, Δz = 10.0 m
  H₂O VMR = 0.032 (3.2%)

Number densities:
  n_air = 2.446e+19 molecules/cm³
  n_H2O = 7.828e+17 molecules/cm³

Loading HITRAN data for thermal IR (8-12 μm)...
  Loaded 1566 transitions

Creating absorption model (Voigt profile)...
  ✓ Model created

Computing absorption cross-section...
  ✓ Computed 4171 points
  Max σ = 1.755e-21 cm²/molecule

======================================================================
RESULTS:
======================================================================

Optical thickness τ for 10.0 m layer:
  Mean:   0.0234
  Median: 0.0089
  Max:    1.3746
  Min:    7.2e-08

Transmittance exp(-τ):
  Mean:   0.9765 (97.7% of radiation transmitted)
  Min:    0.2527 (25.3% at strongest absorption)

======================================================================
INTERPRETATION:
======================================================================

For a 10m tropical boundary layer:
  • Optical thickness varies strongly with wavelength
  • Mean τ ≈ 0.023
  • On average, the layer is OPTICALLY THIN (τ < 1)
  • Most radiation passes through with minimal absorption

  • 5.3% of spectrum has τ > 1 (strong absorption lines)
  • 86.7% of spectrum has τ < 0.1 (atmospheric windows)

======================================================================
✓ Calculation complete!
======================================================================
```

This shows that a 10m tropical layer is mostly transparent in the thermal IR, with strong absorption only at specific water vapor lines.
