# NK Number Size Distribution Exploration Scripts

This directory contains scripts to analyze TOMAS-15 aerosol number concentration (NK) data from GEOSChem output files.

**Location**: `src/Aerosols/data_exploration/`

## Available Scripts

### Python Version: `explore_NK_number_distribution.py`

**Requirements:**
```bash
pip install netCDF4 numpy matplotlib scipy
```

**Usage:**
```bash
python3 src/Aerosols/data_exploration/explore_NK_number_distribution.py
```

### Julia Version: `explore_NK_julia.jl`

**Requirements:**
```julia
using Pkg
Pkg.add(["NCDatasets", "Plots", "LsqFit"])
```

**Usage:**
```bash
julia src/Aerosols/data_exploration/explore_NK_julia.jl
```

## Additional Validation Scripts

- **`explore_tomas_aerosols.py`** - Initial exploration of TOMAS aerosol species
- **`validate_NK_vs_mass_species.py`** - Validation comparing NK to mass species
- **`investigate_NK_units.py`** - Comprehensive unit testing for NK variable
- **`quick_aerosol_check.jl`** - Quick Julia checks for aerosol data

## What These Scripts Do

Both scripts perform identical analysis and produce the same 5 plots:

1. **NK_01_vertical_profile.png** - Vertical profile by size mode (Aitken/Accumulation/Coarse)
2. **NK_02_size_distributions_multi_altitude.png** - Size distributions at multiple altitudes (900, 800, 500, 300, 200 hPa)
3. **NK_03_heatmap.png** - 2D heatmap showing altitude vs size distribution
4. **NK_04_boundary_layer_detail.png** - Boundary layer (~900 hPa) with **bimodal log-normal fit**
5. **NK_05_multi_location_comparison.png** - Size distributions across 6 global locations

## Analyzed Locations

1. **Central USA** (36.8°N, 262.1°E) - Continental/polluted
2. **Amazon Basin** (14.0°N, 94.4°E) - Biomass burning
3. **Sahara Desert** (61.1°N, 287.7°E) - Dust-dominated
4. **India/South Asia** (-2.0°N, 155.6°E) - Polluted
5. **Northern China** (-55.3°N, 215.0°E) - Polluted
6. **South Pacific** (-21.6°N, 184.4°E) - **Clean marine background**

## Key Features

### Unit Conversion
- **Input**: NK stored as **#/kg** (particles per kilogram of air)
- **Metadata claim**: "mol mol-1 dry" (INCORRECT!)
- **Conversion**: N(#/cm³) = NK(#/kg) × ρ_air(kg/cm³)
- **Air density**: ρ_air = P/(R_specific × T), where R_specific = 287.05 J/(kg·K)

### Size Modes
- **Aitken mode**: < 0.1 μm (bins 1-5)
- **Accumulation mode**: 0.1-1.0 μm (bins 6-10)
- **Coarse mode**: > 1.0 μm (bins 11-15)

### Bimodal Fitting
Plot 4 (boundary layer) includes a sophisticated **bimodal log-normal fit** that separates:
- **Mode 1 (Aitken)**: Smaller particles (~0.2 μm median diameter)
- **Mode 2 (Accumulation)**: Larger particles (~0.8 μm median diameter)

Example results for Central USA at 897 hPa:
- Aitken: N = 6.75×10⁴ #/cm³, d_med = 0.200 μm, σ_g = 1.789 (75% of total)
- Accumulation: N = 2.24×10⁴ #/cm³, d_med = 0.825 μm, σ_g = 1.171 (25% of total)

## Output

All plots are saved to `test/aerosol_exploration_output/` in PNG format (150 DPI).

**Note**: Make sure to run these scripts from the repository root directory so that relative paths work correctly:
```bash
# From repository root:
python3 src/Aerosols/data_exploration/explore_NK_number_distribution.py
julia src/Aerosols/data_exploration/explore_NK_julia.jl
```

## Validation

Results have been validated against:
- **Typical atmospheric values**: 10³-10⁴ #/cm³ for continental air ✓
- **Clean marine background**: South Pacific shows ~3,770 #/cm³ (23× lower than continental) ✓
- **Physical constraints**: All converted values are physically reasonable ✓

## Technical Details

### TOMAS-15 Aerosol Model
- **Bins**: 15 logarithmic size sections
- **Size range**: 10 nm - 10 μm dry diameter
- **Spacing**: Even in log-space (factor of √10 per 3 bins)
- **Reference**: Adams & Seinfeld (2002) JGR

### Data File
- **Format**: NetCDF4 (cubed-sphere GEOS-Chem output)
- **Grid**: 24×24×6 cubed-sphere faces
- **Levels**: 72 vertical levels
- **Variables**: NK01-NK15 (number concentration by size bin)

## Performance

- **Python version**: ~10-15 seconds
- **Julia version**: ~15-20 seconds (first run with compilation), ~5-10 seconds (subsequent runs)

Both produce identical results within numerical precision.

## Citation

If you use this analysis in publications, please cite:
- Adams, P. J., & Seinfeld, J. H. (2002). Predicting global aerosol size distributions in general circulation models. *Journal of Geophysical Research*, 107(D19), 4370. doi:10.1029/2001JD001010
