# Data Exploration Tools

This subfolder contains scripts for exploring and validating TOMAS-15 aerosol data from GEOSChem output files.

## Contents

### Main Exploration Scripts
- **`explore_NK_number_distribution.py`** - Python script for comprehensive NK analysis with bimodal fitting
- **`explore_NK_julia.jl`** - Julia version of the NK exploration script (identical functionality)
- **`README.md`** - Detailed documentation for the exploration scripts

### Validation Scripts
- **`explore_tomas_aerosols.py`** - Initial TOMAS species exploration
- **`validate_NK_vs_mass_species.py`** - Cross-validation of NK vs reconstructed number from mass
- **`investigate_NK_units.py`** - Unit testing and validation for NK variable
- **`quick_aerosol_check.jl`** - Quick Julia-based data checks

## Purpose

These scripts were used to:
1. Understand the TOMAS-15 aerosol size distribution data structure
2. **Determine true units of NK** (found to be #/kg, not mol/mol as metadata claims)
3. **Validate unit conversion** using air density: N(#/cm³) = NK(#/kg) × ρ_air(kg/cm³)
4. Develop bimodal log-normal fitting for size distributions
5. Compare aerosol loading across different global locations

## Key Findings

- **NK units**: Particles per kilogram of air (#/kg), despite metadata claiming "mol mol-1 dry"
- **Conversion validated**: Converted values (10³-10⁴ #/cm³) match typical atmospheric measurements
- **Size modes**: Successfully separated Aitken and accumulation modes via bimodal fitting
- **Geographic variability**: Clean marine (South Pacific): ~3,800 #/cm³, Continental: ~35,000-90,000 #/cm³

See the main README.md file in this directory for complete documentation.

## Integration Status

The findings from these exploration scripts have been documented in:
- `../TOMAS_NK_UNITS.md` - Comprehensive unit analysis and implementation guide
- Main vSmartMOM aerosol module (integration in progress)

## Usage

Always run from the repository root directory:
```bash
cd /path/to/vSmartMOM.jl
python3 src/Aerosols/data_exploration/explore_NK_number_distribution.py
julia src/Aerosols/data_exploration/explore_NK_julia.jl
```
