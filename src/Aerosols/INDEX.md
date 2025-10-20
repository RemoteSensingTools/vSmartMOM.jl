# Aerosols Module Documentation

This directory contains the aerosol module implementation for vSmartMOM.jl, including comprehensive documentation of the development process.

## Core Implementation Files

- **`Aerosols.jl`** - Main module file
- **`types.jl`** - Type definitions for aerosol structures
- **`readers.jl`** - Data readers for various aerosol formats
- **`optical_properties.jl`** - Aerosol optical property calculations
- **`refractive_index.jl`** - Refractive index handling

## Aerosol Schemes

See the `schemes/` directory:
- **`tomas15.jl`** - TOMAS-15 sectional aerosol model
- **`two_moment.jl`** - Two-moment aerosol scheme

## Documentation Files

### Implementation Documentation
- **`README.md`** - Main module overview and usage guide
- **`AEROSOL_DESIGN.md`** - Detailed design document for aerosol framework
- **`AEROSOL_FRAMEWORK_DESIGN.md`** - Framework architecture and interfaces
- **`AEROSOL_IMPLEMENTATION_SUMMARY.md`** - Summary of implementation status
- **`AEROSOL_NEXT_STEPS.md`** - Roadmap and future development tasks

### Integration Documentation
- **`GEOSCHEM_INTEGRATION_SUMMARY.md`** - GEOSChem data integration guide
- **`TEST_RESULTS_GEOSCHEM.md`** - Test results for GEOSChem integration
- **`TOMAS_NK_UNITS.md`** - **Critical**: NK variable units analysis and conversion procedures

### Data Exploration
- **`data_exploration/`** - Scripts and tools for analyzing aerosol data
  - Python and Julia exploration scripts
  - Validation and unit testing tools
  - Comprehensive plotting and analysis

## Quick Start

### Understanding NK Units (IMPORTANT!)

Before working with TOMAS aerosol data, read **`TOMAS_NK_UNITS.md`**. Key finding:
- **NK is stored as #/kg** (particles per kilogram of air)
- **NOT mol/mol** as metadata incorrectly claims
- **Conversion**: N(#/cm³) = NK(#/kg) × ρ_air(kg/cm³)

### Running Data Exploration Scripts

From repository root:
```bash
# Python
python3 src/Aerosols/data_exploration/explore_NK_number_distribution.py

# Julia
julia src/Aerosols/data_exploration/explore_NK_julia.jl
```

See `data_exploration/README.md` for details.

## Development Status

✅ **Completed**:
- Type system design
- TOMAS-15 scheme structure
- NK unit analysis and validation
- Data exploration tools (Python & Julia)
- Bimodal size distribution fitting
- GEOSChem reader infrastructure

🚧 **In Progress**:
- Optical property calculations
- Refractive index database integration
- Volume-weighted mixing rules

📋 **Planned**:
- Full vSmartMOM radiative transfer integration
- AERONET validation
- Additional aerosol schemes

## Key Findings from Data Exploration

1. **NK Units**: Definitively determined to be #/kg, validated against atmospheric observations
2. **Size Distributions**: Successfully fitted bimodal log-normal distributions (Aitken + Accumulation modes)
3. **Geographic Variability**: 
   - Clean marine (South Pacific): ~3,800 #/cm³
   - Continental/polluted: ~35,000-90,000 #/cm³
4. **Vertical Structure**: Aerosol concentration peaks in boundary layer, decreases with altitude

## File Organization History

- **October 20, 2025**: Reorganized structure
  - Moved exploration scripts to `data_exploration/`
  - Moved documentation files from root to `src/Aerosols/`
  - Created comprehensive READMEs
  - See `../../FILE_REORGANIZATION_SUMMARY.md` for details

## Contributing

When adding new aerosol schemes or functionality:
1. Follow the type structure in `types.jl`
2. Document any unit conversions carefully (see TOMAS_NK_UNITS.md as example)
3. Add exploration/validation scripts to `data_exploration/`
4. Update this index when adding new files

## References

- Adams, P. J., & Seinfeld, J. H. (2002). Predicting global aerosol size distributions in general circulation models. *JGR*, 107(D19), 4370.
- GEOSChem TOMAS documentation: http://wiki.seas.harvard.edu/geos-chem/
