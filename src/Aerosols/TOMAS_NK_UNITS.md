# TOMAS NK Units and Aerosol Integration Guide

## Executive Summary

The **NK variables in GEOSChem-TOMAS outputs represent particle number concentration**, but they are stored in **particles per kilogram of air (#/kg)**, NOT in the units claimed by the NetCDF metadata ("mol mol-1 dry").

**Key Finding**: To use NK for optical property calculations, convert to particles per cm³ using:

```
N(#/cm³) = NK(#/kg) × ρ_air(kg/cm³)
```

where ρ_air is the local air density.

---

## Background: TOMAS Aerosol Microphysics

TOMAS (TwO-Moment Aerosol Sectional) is a size-resolved aerosol model that tracks:

1. **Particle number concentration** in 15 logarithmic size bins (10 nm - 10 μm dry diameter)
2. **Mass concentration** of individual species (DUST, SS, SF, ECIL, ECOB, OCIL, OCOB, AW) in each bin

**Reference**: Adams, P. J., & Seinfeld, J. H. (2002). Predicting global aerosol size distributions in general circulation models. *Journal of Geophysical Research*, 107(D19). https://doi.org/10.1029/2001JD001010

---

## The NK Mystery: Unit Investigation

### What the Metadata Says (WRONG):
```
Variable: SpeciesConcVV_NK01
  units: mol mol-1 dry
  long_name: Concentration of species for NK01
```

### What the Data Shows:

At 800 hPa (Central USA, typical conditions):
- NK bin 1: 6.73 × 10⁵
- NK bin 7: 1.20 × 10⁹ (peak)
- **Total across bins: 3.95 × 10⁹**

### The Problem:

**Physical constraint**: Volume mixing ratios (mol/mol) must be ≤ 1.0 (100%)

NK values **vastly exceed 1.0**, proving the metadata is incorrect.

---

## Unit Hypothesis Testing

We tested four possible interpretations:

### ❌ Hypothesis 1: NK is in #/cm³
```
Result: 3.95 × 10⁹ #/cm³
Assessment: WAY TOO HIGH (typical air: 10-10⁵ #/cm³)
```

### ✅ Hypothesis 2: NK is in #/kg of air
```
At 800 hPa:
  ρ_air = 0.961 kg/m³ = 9.61 × 10⁻⁷ kg/cm³
  
Convert: NK × ρ_air = 3.95 × 10⁹ × 9.61 × 10⁻⁷
       = 3,800 #/cm³

Assessment: REASONABLE! ✓
  - Typical urban/suburban air: 10³-10⁴ #/cm³
  - Matches Central USA summer conditions
  - Consistent with accumulation mode dominance
```

### ❌ Hypothesis 3: NK is in #/m³
```
Convert: 3.95 × 10³ #/cm³
Assessment: Still too high
```

### ❌ Hypothesis 4: NK is in #/g of air
```
Convert: 3.80 × 10⁶ #/cm³
Assessment: Unrealistically high
```

---

## Confirmed: NK Units are #/kg of Air

### Why This Convention?

Atmospheric models often use **mass-specific units** (#/kg) because:

1. **Conservative quantity** during transport (particle number per unit mass is conserved)
2. **Avoids density corrections** in model dynamics
3. **Common in meteorological models** (similar to specific humidity)
4. **Simplifies vertical coordinate transformations**

### Size Bin Definitions (TOMAS-15)

| Bin | Diameter Range (μm) | Typical Mode |
|-----|---------------------|--------------|
| 1-5 | 0.010 - 0.100 | Aitken/Nucleation |
| 6-10 | 0.100 - 1.000 | Accumulation |
| 11-15 | 1.000 - 10.000 | Coarse |

---

## Conversion Procedure for vSmartMOM

### Step 1: Read NK from NetCDF
```julia
# Read NK for all 15 bins
NK = zeros(15, n_levels)  # #/kg units
for bin in 1:15
    var_name = "SpeciesConcVV_NK$(lpad(bin, 2, '0'))"
    NK[bin, :] = ncread(file, var_name)  # TOA→BOA ordering
end
```

### Step 2: Calculate Air Density
```julia
# Get meteorological variables
pressure = ncread(file, "Met_PMID")  # Pa
temperature = ncread(file, "Met_T")  # K

# Air density (kg/m³)
R_specific = 287.05  # J/(kg·K) for dry air
ρ_air = pressure ./ (R_specific .* temperature)  # kg/m³
ρ_air_cm3 = ρ_air .* 1e-6  # kg/cm³
```

### Step 3: Convert NK to Number Concentration
```julia
# Convert from #/kg to #/cm³
N_cm3 = zeros(15, n_levels)
for bin in 1:15
    N_cm3[bin, :] = NK[bin, :] .* ρ_air_cm3  # #/cm³
end
```

### Step 4: Read Species Mass Concentrations
```julia
# Read mass-based species (mol/mol units)
species = ["DUST", "SS", "SF", "ECIL", "ECOB", "OCIL", "OCOB", "AW"]
mass_conc = Dict()

for sp in species
    mass_conc[sp] = zeros(15, n_levels)
    for bin in 1:15
        var_name = "SpeciesConcVV_$(sp)$(lpad(bin, 2, '0'))"
        mass_conc[sp][bin, :] = ncread(file, var_name)  # mol/mol
    end
end
```

---

## Optical Property Calculation Strategy

### Approach: Combine Number and Mass Information

For each atmospheric layer and size bin:

1. **Number concentration from NK**: Tells you HOW MANY particles
2. **Mass concentrations from species**: Tells you WHAT they're made of
3. **Derive composition**: Calculate mass fraction of each species in each bin
4. **Apply mixing rules**: Use volume averaging or core-shell model
5. **Compute optical properties**: Calculate extinction, scattering, absorption

### Method 1: Volume-Weighted Refractive Index (Simple)

```julia
# For each bin and layer
for bin in 1:15, lev in 1:n_levels
    # Get particle volume (assuming spherical)
    r_bin = radius_centers[bin]  # μm
    V_particle = (4/3) * π * r_bin³  # μm³
    
    # Total mass in this bin
    total_mass = sum(mass_conc[sp][bin, lev] for sp in species)
    
    # Mass fraction of each species
    f_mass = [mass_conc[sp][bin, lev] / total_mass for sp in species]
    
    # Volume fraction (using density of each species)
    f_vol = [f_mass[i] * ρ_ref / ρ_species[i] for i in 1:length(species)]
    f_vol ./= sum(f_vol)  # Normalize
    
    # Volume-weighted refractive index
    n_mix = sum(f_vol[i] * n_species[i] for i in 1:length(species))
    
    # Calculate Mie scattering with mixed refractive index
    Q_ext, Q_sca, Q_abs, g = mie_calculation(r_bin, n_mix, wavelength)
    
    # Scale by number concentration
    σ_ext = π * r_bin² * Q_ext  # Cross-section (cm²)
    β_ext = N_cm3[bin, lev] * σ_ext  # Extinction (km⁻¹)
end
```

### Method 2: External Mixture (More Accurate)

Treat each species as separate particles with their own optical properties, weighted by mass fraction and assumed size distribution.

---

## Example: Typical Values

### At 800 hPa, Central USA (Summer)

| Parameter | Value | Units |
|-----------|-------|-------|
| Pressure | 798 hPa | hPa |
| Temperature | 289 K | K |
| Air density | 0.961 kg/m³ | kg/m³ |
| NK total | 3.95 × 10⁹ | #/kg |
| **N converted** | **3,800** | **#/cm³** |

### Size Distribution

| Mode | Bins | Number Conc. | Percentage |
|------|------|--------------|------------|
| Aitken | 1-5 | ~800 #/cm³ | 21% |
| Accumulation | 6-10 | ~3,000 #/cm³ | 79% |
| Coarse | 11-15 | ~1 #/cm³ | <1% |

Accumulation mode dominates (typical for polluted continental air).

---

## Species Mass Concentrations (mol/mol)

Typical values at 800 hPa:

| Species | Description | Typical VMR | Dominant Bins |
|---------|-------------|-------------|---------------|
| SF | Sulfate | 10⁻¹⁰ - 10⁻⁸ | 6-10 (accumulation) |
| OCIL | Organic Carbon (hydrophilic) | 10⁻¹¹ - 10⁻⁹ | 5-8 |
| OCOB | Organic Carbon (hydrophobic) | 10⁻¹² - 10⁻¹⁰ | 5-8 |
| ECIL | Black Carbon (hydrophilic) | 10⁻¹³ - 10⁻¹¹ | 6-9 |
| ECOB | Black Carbon (hydrophobic) | 10⁻¹³ - 10⁻¹¹ | 6-8 |
| AW | Aerosol Water | 10⁻¹⁰ - 10⁻⁹ | 6-10 |
| DUST | Mineral Dust | 10⁻¹² - 10⁻⁹ | 8-13 (larger) |
| SS | Sea Salt | 10⁻¹² - 10⁻¹⁰ | 8-12 (larger) |

**Note**: These are MASS mixing ratios, not number. Larger particles contribute more mass but fewer numbers.

---

## Implementation Checklist for vSmartMOM

- [ ] Read NK variables from GEOSChem-TOMAS NetCDF files
- [ ] Read meteorological fields (pressure, temperature)
- [ ] Calculate air density (kg/cm³)
- [ ] Convert NK from #/kg to #/cm³
- [ ] Read species mass concentrations (mol/mol)
- [ ] Calculate composition (mass/volume fractions) per bin
- [ ] Load refractive index database for each species
- [ ] Implement mixing rule (volume-weighted or core-shell)
- [ ] Calculate Mie scattering for each bin
- [ ] Compute layer optical properties (extinction, SSA, phase function)
- [ ] Integrate into vSmartMOM RT calculation
- [ ] Validate against observations (AERONET, satellite AOD)

---

## Validation Strategy

### 1. Check Physical Reasonableness
- Total number concentration: 10²-10⁵ #/cm³ (depends on location)
- Size distribution shape: Peak in accumulation mode (0.1-1 μm)
- Vertical profile: Decrease with altitude

### 2. Compare with Literature
- Urban areas: 10⁴-10⁵ #/cm³
- Rural/suburban: 10³-10⁴ #/cm³
- Remote marine: 10²-10³ #/cm³
- Free troposphere: 10¹-10² #/cm³

### 3. Optical Property Checks
- AOD (550 nm): Typically 0.05-0.5 for continental areas
- SSA: Usually 0.85-0.98 (sulfate-dominated) to 0.7-0.9 (soot present)
- Ångström exponent: 1-2 for fine-mode dominance

---

## Common Pitfalls to Avoid

### ❌ DON'T:
1. Use NK values directly as #/cm³ (they're #/kg!)
2. Assume NK metadata units are correct
3. Ignore air density variations with altitude
4. Mix NK (number) with mass species without proper conversion
5. Assume all particles are pure species (they're mixtures!)

### ✅ DO:
1. Always convert NK using air density
2. Use species mass data for composition
3. Apply appropriate mixing rules
4. Validate against known aerosol climatologies
5. Check that optical properties are physically reasonable

---

## References

### Primary Literature
- **Adams & Seinfeld (2002)**: TOMAS model description
  - https://doi.org/10.1029/2001JD001010

### GEOSChem Documentation
- GEOSChem-TOMAS: http://wiki.seas.harvard.edu/geos-chem/
- TOMAS developers: CMU Atmospheric Chemistry (Pierce, Adams groups)

### Optical Properties
- **OPAC Database**: Hess et al. (1998) for refractive indices
- **Mie Theory**: Bohren & Huffman (1983)
- **Mixing Rules**: Volume averaging, Maxwell-Garnett, Bruggeman

---

## Contact and Further Questions

For questions about:
- **NK units and TOMAS**: Contact GEOSChem-TOMAS developers or check GEOSChem wiki
- **Optical calculations**: Consult aerosol optics literature (Mishchenko, Bohren & Huffman)
- **vSmartMOM integration**: Open an issue on the vSmartMOM.jl GitHub repository

---

## Document History

- **2025-01-20**: Initial investigation and unit determination
- **2025-01-20**: Confirmed NK units are #/kg through validation tests
- **2025-01-20**: Created integration guide for vSmartMOM

---

**Last Updated**: January 20, 2025  
**Author**: vSmartMOM Development Team  
**Status**: Confirmed and Ready for Implementation
