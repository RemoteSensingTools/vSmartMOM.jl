# Aerosols Module - Quick Reference

## 📍 Location
`src/Aerosols/`

## 📚 Documentation Files (9 total)

### Core Documentation
- **INDEX.md** - 🌟 **START HERE** - Complete navigation and quick start guide
- **README.md** - Main module overview and usage
- **TOMAS_NK_UNITS.md** - ⚠️ **CRITICAL** - NK units analysis (#/kg, not mol/mol!)

### Design & Implementation
- **AEROSOL_DESIGN.md** - Detailed design document
- **AEROSOL_FRAMEWORK_DESIGN.md** - Framework architecture
- **AEROSOL_IMPLEMENTATION_SUMMARY.md** - Implementation status
- **AEROSOL_NEXT_STEPS.md** - Roadmap and future tasks

### Integration & Testing
- **GEOSCHEM_INTEGRATION_SUMMARY.md** - GEOSChem data integration
- **TEST_RESULTS_GEOSCHEM.md** - Test results

## 🛠️ Code Files (6 total)
- `Aerosols.jl` - Main module
- `types.jl` - Type definitions
- `readers.jl` - Data readers
- `optical_properties.jl` - Optical calculations
- `refractive_index.jl` - Refractive index handling
- `schemes/tomas15.jl` - TOMAS-15 implementation
- `schemes/two_moment.jl` - Two-moment scheme

## 🔬 Data Exploration (8 files)
Location: `src/Aerosols/data_exploration/`

### Main Scripts
- **explore_NK_number_distribution.py** - Comprehensive Python analysis with bimodal fitting
- **explore_NK_julia.jl** - Identical Julia version

### Validation Scripts
- `explore_tomas_aerosols.py` - Initial exploration
- `validate_NK_vs_mass_species.py` - Cross-validation
- `investigate_NK_units.py` - Unit testing
- `quick_aerosol_check.jl` - Quick checks

### Documentation
- `README.md` - Full documentation for exploration scripts
- `OVERVIEW.md` - Summary of findings

## 🚀 Quick Start

### Read Documentation
```bash
# Start here for complete overview
cat src/Aerosols/INDEX.md

# Critical: Understand NK units before using TOMAS data
cat src/Aerosols/TOMAS_NK_UNITS.md
```

### Run Exploration Scripts
```bash
# From repository root:
python3 src/Aerosols/data_exploration/explore_NK_number_distribution.py
julia src/Aerosols/data_exploration/explore_NK_julia.jl
```

### Key Outputs
- 5 plots generated in `test/aerosol_exploration_output/`
- NK units validated: **#/kg** (not mol/mol!)
- Conversion: `N(#/cm³) = NK(#/kg) × ρ_air(kg/cm³)`

## 📊 Key Findings

1. **NK Units**: Definitively #/kg, validated against observations
2. **Size Distributions**: Bimodal (Aitken + Accumulation)
3. **Geographic Range**: 
   - Clean marine: ~3,800 #/cm³
   - Continental/polluted: ~35,000-90,000 #/cm³
4. **Boundary Layer**: 75% Aitken mode, 25% Accumulation mode

## 🗂️ File Organization

```
src/Aerosols/
├── [9 .md documentation files]
├── [6 .jl code files]
├── data_exploration/
│   ├── [6 exploration scripts]
│   └── [2 .md documentation files]
└── schemes/
    ├── tomas15.jl
    └── two_moment.jl
```

**Total**: 24 files organized in 2 subdirectories

## 🔗 Related Files
- Root: `FILE_REORGANIZATION_SUMMARY.md` - Details of October 2025 reorganization
- Root: `GEOSChem.Custom.20190702_0000z.nc4` - Example data file
- Test: `aerosol_exploration_output/` - Plot outputs

## 💡 Tips

1. **Always** read `TOMAS_NK_UNITS.md` before working with NK data
2. **Run from repo root** to ensure relative paths work
3. **Check INDEX.md** for complete navigation
4. **Explore data first** using the data_exploration scripts before implementing

## 🆘 Need Help?

1. Start with `INDEX.md` for navigation
2. Check `README.md` for module overview
3. See `data_exploration/README.md` for script usage
4. Review specific design docs for implementation details

---
Last updated: October 20, 2025
