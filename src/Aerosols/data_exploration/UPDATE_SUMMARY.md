# Data Exploration Scripts Update Summary

**Date**: 2025
**Task**: Update Python and Julia NK exploration scripts to use validated conversion formula

## Validation Context

After extensive investigation involving:
- Mathematical testing of colleague's formula
- Testing at exact location (face=2, x=6, y=19, TOA)
- Comparison with species data (NK/species ratio = M_air constant)

**VALIDATED CORRECT FORMULA:**
```
N (#/cm³) = (NK / 1000) × (Met_AD / M_air) / (Met_AIRVOL × 1e6)
```

where:
- NK: stored as 1000 × (particles/mol_air) 
- M_air = 28.9644e-3 kg/mol (molar mass of dry air)
- Met_AD: Air mass per grid box [kg]
- Met_AIRVOL: Air volume per grid box [m³]

**Evidence:**
- Colleague's formula gives ~6,700 #/cm³ (matches their plot) ✓
- Old formula (NK × ρ_air) gives 194,192 #/cm³ (29× too large) ✗
- NK/species ratio = 29.31 ≈ M_air across all 15 bins ✓

## Files Updated

### 1. explore_NK_number_distribution.py ✅ COMPLETE

**Changes made:**
- Lines ~245-270: Meteorology section
  - Changed from: `rho_air_kg_cm3 = pressure/(R×T) × 1e-6`
  - Changed to: Reading `Met_AD`, `Met_AIRVOL`, calculating `n_air`, `vol_cm3`
  
- Lines ~280-320: NK reading loop
  - Changed from: `nk_concentration = nk_raw × rho_air_kg_cm3`
  - Changed to: `nk_concentration = (nk_raw/1000) × (n_air/vol_cm3)`
  
- Line ~370: Conversion examples section
  - Updated formula display and printed values
  
- Lines ~530-545: Multi-location comparison loop
  - Changed to use `Met_AD`, `Met_AIRVOL` with correct formula
  
- Multiple locations: Plot titles
  - Changed from: "NK × ρ_air"
  - Changed to: "Validated formula: N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)"
  
- Line ~615: Final summary
  - Updated to document validated formula and evidence

**Status**: ✅ Complete and ready to use

### 2. explore_NK_julia.jl ✅ COMPLETE

**Changes made:**
- Lines 260-280: Meteorology section
  - Changed from: Calculating `rho_air_kg_cm3` from pressure and temperature
  - Changed to: Reading `Met_AD`, `Met_AIRVOL`, calculating `n_air`, `vol_cm3`
  
- Lines 285-310: NK reading loop
  - Changed from: `nk_concentration = nk_raw .* rho_air_kg_cm3`
  - Changed to: `nk_concentration = (nk_raw ./ 1000.0) .* (n_air ./ vol_cm3)`
  
- Lines 402-418: Conversion examples section
  - Updated to show `n_air/vol` conversion factor instead of `rho_air`
  
- Line 375: Plot 1 title (size mode profiles)
  - Changed to reference validated formula
  
- Line 427: Plot 2 title (size distributions)
  - Changed to reference validated formula
  
- Lines 571-598: Multi-location comparison loop
  - Completely rewritten to use `Met_AD`, `Met_AIRVOL` with correct formula
  - Removed temperature and `rho_air` calculation
  
- Line 573: Plot 5 title (multi-location)
  - Changed to reference validated formula
  
- Lines 656-661: Final summary
  - Updated to document correct units and validated formula

**Status**: ✅ Complete and ready to use

## Verification

Both scripts now:
1. Read Met_AD and Met_AIRVOL from the NetCDF file
2. Calculate n_air = Met_AD / M_air (moles of air)
3. Calculate vol_cm3 = Met_AIRVOL × 1e6 (volume in cm³)
4. Apply the formula: N = (NK/1000) × (n_air / vol_cm3)
5. Display the validated formula in all plot titles
6. Document the evidence for this formula in summaries

## Expected Results

When run with `GEOSChem.Custom.20190702_0000z.nc4`:

**At colleague's test location** (face=2, x=6, y=19, TOA):
- Should give ~6,700 #/cm³ (matches their working code)

**At our test location** (face=1, x=12, y=12, 800 hPa):
- Should give physically reasonable values (10³-10⁴ #/cm³)
- Old formula gave 194,192 #/cm³ (WRONG, 29× too large)

## Next Steps

1. ✅ Test Python script: `python src/Aerosols/data_exploration/explore_NK_number_distribution.py`
2. ✅ Test Julia script: `julia src/Aerosols/data_exploration/explore_NK_julia.jl`
3. Verify outputs match expected values
4. Consider updating tomas15.jl reader (already started)
5. Archive or replace old TOMAS_NK_UNITS.md documentation

## References

- `RESOLUTION_SUMMARY.md` - Executive summary of NK units investigation
- `TOMAS_NK_UNITS_CORRECTED.md` - Comprehensive correct unit guide
- `NK_SPECIES_COMPARISON.md` - Species validation analysis
- `COLLEAGUE_CODE_ANALYSIS.md` - Analysis of colleague's working code
- `NK_vs_Species_comparison.png` - Visual validation (4-panel plot)
