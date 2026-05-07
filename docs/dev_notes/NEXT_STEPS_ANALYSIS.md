# vSmartMOM.jl: Comprehensive Next Steps Analysis

**Date**: October 14, 2025  
**Branch**: io-update  
**Inspired by**: Oceananigans.jl and CliMA best practices

---

## Executive Summary

vSmartMOM.jl is a mature, scientifically robust radiative transfer package with excellent modularity. The recent IO refactoring (io-update branch) successfully centralizes input handling and removes unsafe eval usage. This analysis identifies strategic next steps across seven key areas to enhance maintainability, performance, usability, and community adoption.

---

## 1. Architecture & Module Organization ⭐⭐⭐

### Current State (Strengths)
- ✅ **Excellent modular structure**: Absorption, Scattering, InelasticScattering, CoreRT, SolarModel, IO
- ✅ **Clear separation of concerns**: Physics modules independent from RT solver
- ✅ **Good GPU abstraction** via Architectures module (borrowed from Oceananigans pattern)
- ✅ **Comprehensive type system** with abstract types and concrete implementations

### Areas for Improvement

#### 1.1 Export Management (Priority: MEDIUM)
**Current Issue**: Exports scattered across multiple module files
```julia
# CoreRT/CoreRT.jl - Line 96+
export model_from_parameters, rt_run, default_parameters
export GaussQuadFullSphere, LambertianSurfaceScalar, LambertianSurfaceSpectrum

# Absorption/Absorption.jl - Line 36+
export AbstractCrossSectionModel, HitranModel, InterpolationModel
export compute_absorption_cross_section, absorption_cross_section
# ... many more
```

**Recommendation** (from CliMA pattern):
```julia
# Create src/Exports.jl
module Exports

# Public API - these are the user-facing functions
export rt_run, model_from_parameters, default_parameters
export parameters_from_yaml, read_parameters
export make_mie_model, compute_aerosol_optical_properties
export read_hitran, make_hitran_model

# Re-export from submodules for convenience
export CPU, GPU, default_architecture
export Stokes_I, Stokes_IQ, Stokes_IQU, Stokes_IQUV
export RadauQuad, GaussLegQuad, GaussQuadFullSphere

end # module Exports
```

**Benefits**:
- Single source of truth for public API
- Easier to manage breaking changes
- Clear documentation of intended user interface
- Matches Oceananigans.jl pattern

#### 1.2 Module File Organization (Priority: LOW)
**Current**: Good hierarchical structure already in place
**Enhancement**: Consider adding `src/Public.jl` following Julia 1.11+ conventions
```julia
# src/Public.jl
public rt_run, model_from_parameters, default_parameters
public parameters_from_yaml, read_parameters
# ... explicit public declarations
```

---

## 2. Testing Infrastructure ⭐⭐⭐⭐⭐

### Current State
- ✅ Comprehensive test coverage across modules
- ✅ Integration with test suite (test/runtests.jl)
- ✅ Benchmark comparisons (6SV1, Natraj)
- ✅ RAMI validation framework

### Critical Gaps & Recommendations

#### 2.1 **Missing: Continuous Integration** (Priority: CRITICAL)
**Current**: No CI/CD visible in repo

**Recommendation**: Add GitHub Actions workflow
```yaml
# .github/workflows/CI.yml
name: CI
on:
  push:
    branches: [main, io-update]
  pull_request:
    branches: [main]
jobs:
  test:
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        julia-version: ['1.9', '1.10', '1.11']
        os: [ubuntu-latest]
        arch: [x64]
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v2
        with:
          version: ${{ matrix.julia-version }}
      - uses: julia-actions/cache@v2
      - uses: julia-actions/julia-buildpkg@v1
      - uses: julia-actions/julia-runtest@v1
      - uses: julia-actions/julia-processcoverage@v1
      - uses: codecov/codecov-action@v4
```

**Also add**:
- `.github/workflows/Documentation.yml` for automated docs builds
- `.github/workflows/CompatHelper.yml` for dependency updates
- Coverage reporting (Codecov/Coveralls)

#### 2.2 **Add Unit Tests for New IO System** (Priority: HIGH)
**Missing tests for**:
```julia
# test/test_IO.jl (NEW FILE NEEDED)
@testset "IO Parsing" begin
    @testset "Surface string parsing" begin
        # Test parse_surface_str with various inputs
        @test parse_surface_str("LambertianSurfaceScalar(0.3)", Float64) isa CoreRT.LambertianSurfaceScalar
        @test parse_surface_str("rpvSurfaceScalar(0.1, 0.5, -0.2, 30.0)", Float32) isa CoreRT.rpvSurfaceScalar
        # Edge cases
        @test_throws AssertionError parse_surface_str("InvalidSurface(0.1)", Float64)
    end
    
    @testset "Mapping tables" begin
        # Verify all map entries resolve correctly
        for key in keys(IO.POLARIZATION_MAP)
            @test IO.POLARIZATION_MAP[key](Float64) isa Scattering.AbstractPolarizationType
        end
    end
    
    @testset "spec_bands parsing" begin
        # Test with and without units
        # Test error handling
    end
end
```

#### 2.3 **Performance Regression Testing** (Priority: MEDIUM)
**Current**: Timing macros present but no regression tracking

**Recommendation**: Add `test/benchmarks/` with PkgBenchmark.jl
```julia
# test/benchmarks/benchmarks.jl
using BenchmarkTools, PkgBenchmark

SUITE = BenchmarkGroup()

SUITE["rt_run"] = BenchmarkGroup(["radiative_transfer"])
params = default_parameters()
model = model_from_parameters(params)
SUITE["rt_run"]["small"] = @benchmarkable rt_run($model)

# Run with: using PkgBenchmark; benchmarkpkg("vSmartMOM")
```

#### 2.4 **Aqua.jl Quality Checks** (Priority: MEDIUM)
Following Oceananigans pattern:
```julia
# test/test_quality.jl
using Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(
        vSmartMOM;
        ambiguities = false,  # Many expected from submodule dispatch
        piracy = true,
        unbound_args = true,
    )
end
```

---

## 3. Documentation Enhancements ⭐⭐⭐⭐

### Current State
- ✅ Excellent module-level documentation
- ✅ New IO docs (Overview, Schema, Examples)
- ✅ Principles page (NEW, great addition!)
- ✅ Literate.jl tutorials
- ✅ Comprehensive docstrings

### Recommendations

#### 3.1 **API Reference Organization** (Priority: HIGH)
**Current**: Docs show methods but could be better organized

**Add structured API reference pages**:
```markdown
# docs/src/api/index.md
# API Reference

## High-Level Interface
- [`rt_run`](@ref)
- [`model_from_parameters`](@ref)
- [`default_parameters`](@ref)

## Input/Output
- [`parameters_from_yaml`](@ref)
- [`read_parameters`](@ref)
- [`read_atmos_profile`](@ref)

## Absorption Module
[Link to detailed page]

## Scattering Module
[Link to detailed page]

## CoreRT Module
[Link to detailed page]
```

#### 3.2 **Add "Quick Start" / "Getting Started" Tutorial** (Priority: HIGH)
**Current**: Example page exists but could be friendlier

**Recommendation**: Create interactive notebook-style tutorial
```julia
# docs/src/pages/tutorials/Tutorial_QuickStart.jl
# # Quick Start: Your First RT Simulation
#
# This 5-minute tutorial shows you how to run your first
# radiative transfer calculation with vSmartMOM.jl

using vSmartMOM

## Step 1: Load default parameters (or from file)
params = default_parameters()

## Step 2: Customize (optional)
# params.sza = 45.0  # Change solar zenith angle

## Step 3: Build model
model = model_from_parameters(params)

## Step 4: Run RT
R, T = rt_run(model)

## Step 5: Visualize
using Plots
plot(model.params.spec_bands[1], R[1, 1, :], 
     xlabel="Wavenumber (cm⁻¹)", ylabel="Reflectance",
     title="Top-of-Atmosphere Reflectance")
```

#### 3.3 **Performance Tips Section** (Priority: MEDIUM)
Document GPU usage, Float32 vs Float64, parallelization

#### 3.4 **Fix Remaining Doc Warnings** (Priority: LOW)
```julia
# From last docs build:
# - Invalid link in Scattering/Types.md ✅ FIXED
# - Missing docs for 191 internal functions (expected for internals)
# - Cross-ref warnings ✅ PARTIALLY FIXED

# Remaining: Add docstrings for commonly-referenced internals
```

---

## 4. Code Quality & Maintainability ⭐⭐⭐

### Current Strengths
- ✅ Consistent coding style
- ✅ Good use of Parameters.jl for structs
- ✅ DocStringExtensions for consistent docstrings
- ✅ Type stability considerations

### Improvements

#### 4.1 **Add `.JuliaFormatter.toml`** (Priority: HIGH)
Following Oceananigans/CliMA pattern:
```toml
# .JuliaFormatter.toml
style = "blue"
indent = 4
margin = 100
always_for_in = true
whitespace_typedefs = true
whitespace_ops_in_indices = true
remove_extra_newlines = true
```

Add formatting CI check:
```yaml
# .github/workflows/Format.yml
name: Format Check
on: [push, pull_request]
jobs:
  format:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v2
      - run: |
          julia -e 'using Pkg; Pkg.add("JuliaFormatter")'
          julia -e 'using JuliaFormatter; format(".", verbose=true) || error("Formatting required")'
```

#### 4.2 **Address TODO Comments** (Priority: MEDIUM)
**Found 50+ TODO comments**. Prioritize:

**HIGH PRIORITY**:
```julia
# src/CoreRT/tools/atmo_prof.jl:148
# TODO: This needs a VCD_dry weighted average!

# src/CoreRT/rt_run.jl:119
# TODO: if RS_type!=noRS, create inelastic arrays

# src/CoreRT/CoreKernel/interaction_hdrf.jl:38
# TODO: Use Radau quadrature and include insolation
```

**Create issues** for each significant TODO and link them:
```julia
# TODO(#123): This needs a VCD_dry weighted average!
```

#### 4.3 **Type Annotations for Public Functions** (Priority: MEDIUM)
Enhance type stability documentation:
```julia
# Before:
function rt_run(model; i_band = 1)

# After:
function rt_run(model::vSmartMOM_Model; i_band::Integer = 1)::Tuple{Array, Array}
```

---

## 5. Performance Optimization ⭐⭐⭐⭐

### Current State
- ✅ GPU support via CUDA.jl and KernelAbstractions
- ✅ Batched operations
- ✅ TimerOutputs for profiling
- ✅ Float32/Float64 support

### Recommendations

#### 5.1 **Add Performance Benchmarks Section in Docs** (Priority: HIGH)
```markdown
# docs/src/pages/performance.md
# Performance Guide

## Hardware Configurations
- CPU: Single core, multi-core, server
- GPU: NVIDIA, AMD (via AMDGPU.jl future support)

## Float Precision Trade-offs
- Float64: ~2x memory, slower, more accurate
- Float32: ~2x faster, GPU-friendly, sufficient for most RT

## Batch Processing
...
```

#### 5.2 **Profile-Guided Optimization** (Priority: MEDIUM)
Add profiling utilities:
```julia
# src/Profiling.jl (NEW)
module Profiling

using Profile, ProfileView

function profile_rt_run(params=default_parameters(); nruns=10)
    model = model_from_parameters(params)
    
    # Warm-up
    rt_run(model)
    
    # Profile
    Profile.clear()
    @profile for _ in 1:nruns
        rt_run(model)
    end
    
    return Profile.print(), ProfileView.view()
end

export profile_rt_run
end
```

#### 5.3 **Memory Allocation Audit** (Priority: LOW)
Run allocation profiler on critical paths:
```julia
using AllocCheck
@check_allocs rt_run(model)
```

---

## 6. Extensibility & User Experience ⭐⭐⭐⭐

### Recommendations

#### 6.1 **Extension Points via Package Extensions** (Priority: HIGH)
Julia 1.9+ feature - create optional integrations:
```julia
# ext/vSmartMOMPlotsExt.jl
module vSmartMOMPlotsExt

using vSmartMOM
using Plots

# Add convenience plotting recipes
@recipe function f(model::vSmartMOM.CoreRT.vSmartMOM_Model)
    xguide --> "Wavenumber (cm⁻¹)"
    yguide --> "Reflectance"
    label --> "TOA Reflectance"
    # ...
end

end # module
```

In Project.toml:
```toml
[extensions]
vSmartMOMPlotsExt = "Plots"

[compat]
Plots = "1"
```

#### 6.2 **Interactive Examples with Pluto.jl** (Priority: MEDIUM)
Create interactive notebooks:
```julia
# docs/notebooks/interactive_rt.jl
### A Pluto.jl notebook ###
# Interactive RT parameter explorer
using Pluto, vSmartMOM, PlutoUI

@bind sza Slider(0:90, default=30)
@bind albedo Slider(0:0.01:1, default=0.3)

# Real-time RT calculation and plot
```

#### 6.3 **Parameter Validation Layer** (Priority: MEDIUM)
Add comprehensive validation:
```julia
# src/Validation.jl (NEW)
module Validation

function validate_parameters(params::vSmartMOM_Parameters)
    errors = String[]
    
    # Physical constraints
    0 ≤ params.sza ≤ 90 || push!(errors, "SZA must be in [0, 90]")
    all(0 .≤ params.vza .≤ 90) || push!(errors, "VZA must be in [0, 90]")
    
    # Consistency checks
    length(params.brdf) == length(params.spec_bands) || 
        push!(errors, "Number of surfaces must match number of bands")
    
    !isempty(errors) && throw(ArgumentError(join(errors, "\n")))
    return true
end

export validate_parameters
end
```

#### 6.4 **Examples Gallery** (Priority: LOW)
Create `examples/` directory with common use cases:
```
examples/
├── satellite_retrieval.jl
├── multi_band_processing.jl
├── aerosol_sensitivity.jl
├── gpu_acceleration.jl
└── custom_surface_brdf.jl
```

---

## 7. Community & Ecosystem ⭐⭐⭐⭐⭐

### Critical for Adoption

#### 7.1 **GitHub Repository Health** (Priority: CRITICAL)
**Add missing files**:
```
.github/
├── ISSUE_TEMPLATE/
│   ├── bug_report.md ✅ EXISTS
│   ├── feature_request.md (NEW)
│   └── question.md (NEW)
├── PULL_REQUEST_TEMPLATE.md (NEW)
├── CONTRIBUTING.md (NEW)
└── workflows/
    ├── CI.yml (NEW)
    ├── Documentation.yml (NEW)
    ├── CompatHelper.yml (NEW)
    └── TagBot.yml (NEW)
```

#### 7.2 **CONTRIBUTING.md** (Priority: HIGH)
```markdown
# Contributing to vSmartMOM.jl

## Getting Started
1. Fork the repository
2. Create a branch (`git checkout -b feature/amazing-feature`)
3. Make changes
4. Run tests (`julia --project -e 'using Pkg; Pkg.test()'`)
5. Format code (`julia -e 'using JuliaFormatter; format(".")'`)
6. Submit PR

## Development Setup
...

## Coding Standards
- Follow Julia style guide
- Add docstrings for public functions
- Include tests for new features
- Keep type-stable

## Where to Contribute
- See [open issues](link)
- Check TODO comments in code
- Improve documentation
- Add examples
```

#### 7.3 **Community Engagement** (Priority: HIGH)
- **Set up Discussions**: Enable GitHub Discussions for Q&A
- **Zulip/Slack**: Join Julia Earth/Climate community channels
- **JuliaHub registration**: Ensure package is properly registered
- **CITATION.bib**: Add citation information
```bibtex
@software{vSmartMOM,
  author = {Sanghavi, Suniti and Frankenberg, Christian and ...},
  title = {vSmartMOM.jl: Vector Smart Matrix-Operator Method},
  year = {2024},
  url = {https://github.com/RemoteSensingTools/vSmartMOM.jl},
  doi = {...}
}
```

#### 7.4 **Interoperability** (Priority: MEDIUM)
**Connect with Julia ecosystem**:
- **SciML**: Integrate with DifferentialEquations.jl for time-evolution
- **GeoStats.jl**: Spatial processing of RT output
- **JuliaClimate**: Collaborate with other atmospheric packages
- **ModelingToolkit.jl**: Symbolic modeling capabilities

---

## 8. Specific Technical Debt Items

### From Code Analysis

#### 8.1 **Inelastic Scattering** (Priority: MEDIUM)
Multiple TODOs for inelastic/Raman scattering:
```julia
# Many files reference: "TODO: if RS_type!=noRS, create..."
# Indicates incomplete Raman implementation
```
**Action**: Create tracking issue, break into subtasks

#### 8.2 **Observer Altitude** (Priority: LOW)
```julia
# Currently documented as "not used"
# Only TOA observations supported
```
**Action**: Either implement or remove from parameters

#### 8.3 **Canopy RT** (Priority: LOW)
`rt_run_canopy.jl` appears experimental
**Action**: Document status, add tests or mark as experimental

#### 8.4 **Eval Usage in spec_bands** (Priority: LOW)
```julia
# Currently documented in IO docs ✅
# Consider future: parse ranges without eval
```
**Action**: Track as enhancement for future

---

## 9. Immediate Action Plan (Next 2-4 weeks)

### Week 1: Infrastructure (CRITICAL)
1. ✅ **Merge io-update branch** (assuming tests pass)
2. **Set up CI/CD** (GitHub Actions)
   - CI.yml for testing
   - Documentation.yml for docs
   - CompatHelper.yml
   - TagBot.yml
3. **Add CONTRIBUTING.md**
4. **Enable GitHub Discussions**

### Week 2: Testing & Quality
5. **Add IO unit tests** (test/test_IO.jl)
6. **Add Aqua.jl quality checks**
7. **Set up JuliaFormatter**
8. **Fix critical TODO items** (create GitHub issues for rest)

### Week 3: Documentation
9. **Add Quick Start tutorial**
10. **Restructure API reference**
11. **Add Performance Guide**
12. **Fix remaining doc warnings**

### Week 4: Polish
13. **Create examples gallery**
14. **Add CITATION.bib**
15. **Announce on Julia Discourse**
16. **Submit to JuliaHub showcase**

---

## 10. Long-Term Vision (6-12 months)

### Scientific Enhancements
- **Complete Raman scattering** implementation
- **Add VLIDORT comparison** benchmarks
- **3D RT capabilities** (stretch goal)
- **Bayesian inversion** examples with Turing.jl

### Software Engineering
- **Package extensions** for plotting, sensitivity analysis
- **GPU performance optimization** (profile-guided)
- **ARM/Apple Silicon** native support
- **Python wrapper** (PyCall/PythonCall) for broader adoption

### Community
- **Workshop/tutorial** at Julia conference
- **Paper in JOSS** (Journal of Open Source Software)
- **Integration** with CliMA ecosystem
- **Active contributor community** (5-10 regular contributors)

---

## Comparison with Oceananigans.jl & CliMA Patterns

### What vSmartMOM Does Well (matches best practices)
✅ Modular architecture  
✅ GPU abstraction layer  
✅ Comprehensive type system  
✅ Good documentation structure  
✅ Physics separated from infrastructure  

### What to Adopt from Oceananigans/CliMA
📋 Structured exports (single source of truth)  
📋 Comprehensive CI/CD pipeline  
📋 JuliaFormatter integration  
📋 Aqua.jl quality checks  
📋 Package extensions (Julia 1.9+)  
📋 Interactive examples (Pluto.jl)  
📋 Community engagement infrastructure  

---

## Conclusion

vSmartMOM.jl is a **scientifically mature and architecturally sound** package. The io-update branch represents excellent progress toward modern Julia best practices. The highest-leverage improvements are:

1. **Infrastructure**: CI/CD, formatting, testing framework
2. **Documentation**: Quick start guide, performance tips, API reference
3. **Community**: CONTRIBUTING.md, discussions, examples
4. **Quality**: Unit tests for new code, address critical TODOs

With these enhancements, vSmartMOM.jl can become a **flagship package** in the Julia atmospheric science ecosystem, on par with Oceananigans.jl in ocean modeling.

---

**Next immediate action**: Review this document with team, prioritize items, create GitHub project board to track progress.
