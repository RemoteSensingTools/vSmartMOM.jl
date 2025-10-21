# GEOSChem-TOMAS NK Units: Resolution Summary

**Date:** October 21, 2025  
**Issue:** Confusion about NK variable units in GEOSChem-TOMAS output  
**Resolution:** NK is in `1000 × (particles/mol_air)`, NOT `#/kg`

---

## The Debate

### Initial Hypothesis (WRONG)
NK is in `#/kg` (particles per kg of air)

**Conversion:** `N = NK × ρ_air`

**Why it seemed right:**
- Gave reasonable values (~3,800-5,000 #/cm³)
- Simple and intuitive
- Worked at some test locations

### Correct Interpretation (VALIDATED)
NK is in `1000 × (#/mol_air)` (1000 times particles per mole of air)

**Conversion:** `N = (NK/1000) × (n_air/vol)` where `n_air = Met_AD/M_air`

**Why it's correct:**
- Colleague's working code uses this formula
- Gives correct values matching their plots (~6,700 #/cm³)
- Consistent with mass species (ratio = M_air)
- Size distributions match when properly scaled

---

## The Smoking Gun

### Ratio Test

When comparing NK with sum of mass species (converted to particle number):

```
N_NK / N_species = 29.31 ≈ M_air (28.96 g/mol)
```

This ratio is **constant across all size bins**, proving:
1. The shapes are identical (same physical distribution)
2. The difference is purely a unit convention
3. NK includes a 1/M_air factor that mass species don't have

### Why the M_air Factor Appears

**NK formula:**
```julia
N = (NK/1000) × (Met_AD / M_air) / (Met_AIRVOL × 1e6)
                 └─────┬─────┘
                   divides by M_air
```

**Species formula:**
```julia
M = MW × species × (Met_AD / M_air) / Met_AIRVOL
N = M / m_particle
  = M / (ρ × V_particle)
# The M_air factor cancels out in mass calculations
```

---

## Validation Evidence

### 1. Direct Comparison (Test Location)
- **File:** GEOSChem.Custom.20190702_0000z.nc4
- **Location:** Face 2, x=6, y=19, level 1 (TOA)
- **Met_AD:** 2.244×10¹³ kg
- **Met_AIRVOL:** 2.027×10¹³ m³

**Results:**
| Method | Total N (#/cm³) | Assessment |
|--------|----------------|------------|
| Colleague formula | 6,704.52 | ✅ Matches their plot (~6,000) |
| Wrong formula (NK×ρ) | 194,192 | ❌ 29× too large |

### 2. Mathematical Proof

Working backwards from correct answer:
```julia
N_correct = 6,704.52 #/cm³
NK_total = 1.754×10¹¹

# Colleague's formula
n_air = 2.244e13 / 28.9644e-3 = 7.748e14 mol
vol = 2.027e13 × 1e6 = 2.027e19 cm³
N = (NK/1000) × (n_air/vol) = 6,704.52 ✅

# Wrong formula  
ρ_air = 2.244e13 / 2.027e13 = 1.107 kg/m³
N = NK × ρ_air × 1e-6 = 194,192 ❌
```

### 3. Size Distribution Match

See `NK_vs_Species_comparison.png`:
- Unscaled: NK is 29× larger than species
- Scaled by M_air: Curves overlap perfectly
- Ratio plot: Constant at ~29 across all bins

---

## Common Mistakes to Avoid

### ❌ DON'T: Use `N = NK × ρ_air`
- Gives 29× wrong answer
- NK is per MOLE not per KG

### ❌ DON'T: Forget the /1000 factor
- NK is scaled by 1000 in the file
- Must divide by 1000 in conversion

### ❌ DON'T: Mix up volume units
- NK: use cm³ for concentration
- Mass species: use m³ for concentration

### ✅ DO: Use the validated formula
```julia
N (#/cm³) = (NK / 1000) × (Met_AD / M_air) / (Met_AIRVOL × 1e6)
```

### ✅ DO: Check consistency
- NK total vs species total should differ by ~M_air
- Size distribution shapes should match when scaled

---

## Why GEOSChem Uses This Convention

### Internal Representation
GEOSChem stores all concentrations as "mol/mol dry air" for consistency:
- **Mass species:** True mol/mol (mol_species/mol_air)
- **Number species (NK):** Artificial molar mass of 1×10⁻³ N/mol

### The Math
```
If molar_mass_NK = 1×10⁻³ N/mol, then:
  
  mol_NK/mol_air = N/1000 particles/mol_air
  
  Therefore: NK_stored = 1000 × (N/mol_air)
```

### Advantages
- Consistent with GEOSChem's transport/chemistry operators
- All species use same fractional representation
- Integrates seamlessly with model physics

### Disadvantages
- Confusing for users
- Metadata doesn't explain the convention
- Easy to misinterpret as #/kg or other units

---

## Final Recommendations

### For Users

1. **Always use:** `N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)`

2. **Validate your results:**
   - Typical surface: 1,000-50,000 #/cm³
   - Clean marine: 100-1,000 #/cm³  
   - Polluted continental: 10,000-100,000 #/cm³

3. **Cross-check with mass species:**
   - Calculate N from species mass
   - Ratio should be ~M_air (28.96)

### For vSmartMOM Integration

Use the validated conversion in aerosol readers:

```julia
function read_TOMAS_number(ds, idx, idy, idf, idl)
    M_air = 28.9644e-3  # kg/mol
    
    Met_AD = ds["Met_AD"][1, idl, idf, idy, idx]
    Met_AIRVOL = ds["Met_AIRVOL"][1, idl, idf, idy, idx]
    
    n_air = Met_AD / M_air
    vol_cm3 = Met_AIRVOL * 1e6
    
    N = Float64[]
    for k in 1:15
        NK = ds["SpeciesConcVV_NK$(lpad(k,2,'0'))"][1, idl, idf, idy, idx]
        push!(N, (NK / 1000) * (n_air / vol_cm3))
    end
    
    return N  # #/cm³
end
```

---

## Documentation Files

1. **TOMAS_NK_UNITS_CORRECTED.md** - Complete unit guide
2. **NK_SPECIES_COMPARISON.md** - Detailed validation analysis
3. **COLLEAGUE_CODE_ANALYSIS.md** - Analysis of working code
4. **This file** - Executive summary

---

## Acknowledgments

Thanks to Christian's colleague for:
- Sharing their working Julia code
- Providing plots showing correct values
- Clarifying the mol/mol convention with M_NK = 1×10⁻³ N/mol

This resolved weeks of confusion and validated the correct interpretation!

---

## Key Takeaway

**NK is fundamentally particles-per-mole-of-air (with 1000× scaling), NOT particles-per-kg.**

The confusion arose because:
1. `N = NK × ρ_air` gives reasonable values (by coincidence)
2. Metadata says "mol mol-1" without explaining the artificial molar mass
3. The factor-of-29 discrepancy only becomes obvious when comparing methods

**Trust the colleague's code. Trust the M_air ratio. Trust the mathematics.**

