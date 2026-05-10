# TOMAS NK Units: Definitive Guide

## Executive Summary

**CORRECT INTERPRETATION (VALIDATED):**

The NK variables in GEOSChem-TOMAS output represent **aerosol particle number** in units of:

```
1000 × (particles per mole of air)
```

Despite metadata saying "mol mol-1 dry", NK uses a special convention where particle numbers are normalized to **moles of air** (not kg of air).

## Conversion Formula (CORRECT)

```julia
N (#/cm³) = (NK / 1000) × (n_air / volume)
```

where:
```julia
n_air = Met_AD / M_air        # moles of air in grid cell [mol]
M_air = 28.9644e-3           # kg/mol (molar mass of dry air)  
volume = Met_AIRVOL × 1e6    # cm³
```

### Full Expansion

```julia
N (#/cm³) = (NK / 1000) × (Met_AD / M_air) / (Met_AIRVOL × 1e6)
```

## WARNING: Common Mistakes

### ❌ WRONG Formula (gives 29× too large):
```julia
N = NK × ρ_air  # DON'T USE THIS!
```

This formula is **mathematically incorrect** because:
- NK is per **mole** of air, not per **kg** of air
- Missing the factor of 1/1000
- Results in values 28.96× (= M_air) too large

### ❌ WRONG Interpretation:
"NK is in #/kg" - This is NOT correct despite giving reasonable-looking numbers

## Mathematical Proof

### Test Case
Location: GEOSChem.Custom.20190702_0000z.nc4, face=2, x=6, y=19, level=1

**Data:**
- Met_AD = 2.244×10¹³ kg
- Met_AIRVOL = 2.027×10¹³ m³
- NK01 = 7.472×10⁵
- Sum(NK01-NK15) = 1.754×10¹¹

**Correct Formula Result:**
```julia
n_air = 2.244e13 / 28.9644e-3 = 7.748e14 mol
vol = 2.027e13 × 1e6 = 2.027e19 cm³
N_total = (1.754e11 / 1000) × (7.748e14 / 2.027e19)
        = 6,704.52 #/cm³  ✅ PHYSICALLY REASONABLE
```

**Wrong Formula Result:**
```julia
ρ_air = 2.244e13 / 2.027e13 = 1.107 kg/m³ = 1.107e-6 kg/cm³
N_total = 1.754e11 × 1.107e-6
        = 194,192 #/cm³  ❌ 29× TOO LARGE
```

### Validation: Comparison with Mass Species

When comparing NK with aerosol mass species (SO4, DUST, AW, etc.) converted to particle number:

```
N_NK / N_species = 29.31 ≈ M_air
```

This factor of M_air arises because:
- NK uses: particles per **mole** of air
- Species use: mass per mass (mol/mol, then converted via molecular weight and density)
- The 1/M_air factor in NK formula creates this difference

**When properly scaled by M_air, the size distributions have identical shapes**, confirming they represent the same physical quantity with different unit conventions.

## Full Example Code

```julia
using NCDatasets

function convert_NK_to_concentration(ds, idx, idy, idf, idl)
    """
    Convert NK values to particle number concentration
    
    Args:
        ds: NCDataset object
        idx, idy, idf: Spatial indices (1-indexed)
        idl: Vertical level (1-indexed, where 1=surface after reversal)
    
    Returns:
        Array of concentrations [#/cm³] for 15 bins
    """
    
    # Convert to 0-indexed for Python/NetCDF
    i, j, f, l = idx-1, idy-1, idf-1, idl-1
    
    # Read meteorology
    Met_AD = ds["Met_AD"][1, l, f, j, i]      # kg
    Met_AIRVOL = ds["Met_AIRVOL"][1, l, f, j, i]  # m³
    
    # Constants
    M_air = 28.9644e-3  # kg/mol
    
    # Conversion factors
    n_air = Met_AD / M_air      # mol
    vol_cm3 = Met_AIRVOL * 1e6  # cm³
    
    # Read and convert all NK values
    N = Float64[]
    for k in 1:15
        NK = ds["SpeciesConcVV_NK$(lpad(k,2,'0'))"][1, l, f, j, i]
        push!(N, (NK / 1000) * (n_air / vol_cm3))
    end
    
    return N
end
```

## Size Distribution: dN/dlogD

TOMAS-15 has logarithmically-spaced bins. To get dN/dlogD:

```julia
# Bin boundaries (diameter in nm)
D_bounds = [10.0, 15.8, 25.1, 39.8, 63.1, 100.0, 158.5, 251.2, 
            398.1, 631.0, 1000.0, 1584.9, 2511.9, 3981.1, 6309.6, 10000.0]

# Calculate Δlog₁₀(D) for each bin
delta_logD = log10.(D_bounds[2:end]) .- log10.(D_bounds[1:end-1])

# Convert N to dN/dlogD
dN_dlogD = N ./ delta_logD  # cm⁻³
```

## Typical Values

| Location | Pressure | N_total (#/cm³) | Notes |
|----------|----------|----------------|-------|
| Central USA | 800 hPa | ~35,000 | Continental, polluted |
| South Pacific | 800 hPa | ~3,800 | Marine, clean |
| TOA (test case) | ~50 hPa | ~6,700 | Upper atmosphere |

## Historical Note: Why the Confusion?

### What Went Wrong

Initial analysis suggested NK was in "#/kg" because the formula `N = NK × ρ_air` gave seemingly reasonable results (~3,800-5,000 #/cm³ at 800 hPa).

### Why This Was Wrong

1. The results were only "reasonable" by coincidence
2. When tested at different levels or compared with mass species, the discrepancy appeared
3. The ratio NK/species was always ~29 ≈ M_air, which is the smoking gun
4. The colleague's formula `(NK/1000) × (n_air/vol)` gives correct values matching their plots

### Lesson Learned

Always validate against multiple test cases and independent methods. A formula that works at one location may be fundamentally wrong!

## Unit Convention Explanation

### Why GEOSChem Uses This Convention

GEOSChem stores all concentrations as "mol/mol dry air" for consistency. For number concentrations:

1. Can't use molecular weight (particles don't have molar mass)
2. Solution: Define artificial "molar mass" = 1×10⁻³ N/mol
3. Then: mol_NK/mol_air = (N/1000) particles/mol_air
4. Stored as: NK = 1000 × N/mol_air

### Advantages
- Consistent with all other GEOSChem species
- Preserves fractional representation
- Works with GEOSChem's transport and chemistry operators

### Disadvantages  
- Confusing for users expecting #/kg or #/cm³
- Metadata doesn't explain the convention
- Easy to misinterpret

## Related Variables

### Mass Species

Variables like `SpeciesConcVV_SO4`, `SpeciesConcVV_DUST`, etc. use **true mol/mol dry air**:

```julia
# For mass species (e.g., SO4)
mol_species = species_value × n_air         # mol
mass = mol_species × MW_species            # g
# Then convert to particle number using density and size
```

These are fundamentally different from NK and should not be converted the same way.

### Meteorology

- `Met_AD`: Air mass per grid box [kg]
- `Met_AIRVOL`: Air volume per grid box [m³]
- `Met_DELP`: Pressure thickness [Pa]
- `Met_PS2WET`: Surface pressure [Pa]

## References

1. **TOMAS Model:** Adams & Seinfeld (2002), JGR, doi:10.1029/2001JD001010
2. **GEOSChem TOMAS:** http://wiki.seas.harvard.edu/geos-chem/index.php/TOMAS_aerosol_microphysics
3. **Validation:** NK_SPECIES_COMPARISON.md (this repository)

## See Also

- `NK_SPECIES_COMPARISON.md` - Detailed comparison showing NK and species match when scaled by M_air
- `data_exploration/NK_vs_Species_comparison.png` - Visualization of distributions
- `examples/geoschem_integration.jl` - Working example code

## Summary

✅ **Correct formula:** `N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)`

✅ **Unit interpretation:** NK is in `1000 × (#/mol_air)`

✅ **Validated by:**
   - Colleague's working code with plots
   - Comparison with mass species (ratio = M_air)
   - Multiple test locations

❌ **Don't use:** `N = NK × ρ_air` (gives 29× wrong answer)

❌ **NK is NOT in #/kg** (despite this seeming to work at some locations)
