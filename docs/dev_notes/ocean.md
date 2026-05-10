# Ocean Radiative Transfer Implementation Plan

## Status: Phase 1 — Start with OceanOptics.jl; remaining phases deferred

---

## 1. Motivation and Context

vSmartMOM.jl currently handles atmospheric RT with various surface BRDFs (Lambertian, RPV, RossLi, CoxMunk, Canopy). The `CoxMunkSurface` models wind-roughened ocean *reflection* but treats the ocean as opaque — no light penetrates the water column. This plan adds full coupled atmosphere-ocean RT:

- Fresnel transmission into/out of the water column
- Multi-layer ocean volume scattering/absorption by water constituents
- Chlorophyll fluorescence emission
- APAR (Absorbed Photosynthetically Active Radiation) computation per layer
- Ocean floor lower boundary condition

### References
- **Fell (1997)**: "Validierung eines Modells zur Simulation des Strahlungstransportes in Atmosphäre und Ozean" — PhD thesis, FU Berlin. Located at `docs/papers/1997_Thesis_Fell.pdf`. Uses the Matrix Operator Method (identical to vSmartMOM's adding-doubling) for coupled atmosphere-ocean systems. Chapter 3 covers the ocean floor and water surface treatment. Section 2.2.9 covers fluorescence. Chapter 5 covers water constituent optical properties.
- **Ocean Optics Web Book**: https://www.oceanopticsbook.info/ — comprehensive online reference for ocean optics theory and measurements
- **Gordon (1992)**: Bio-optical model relating phytoplankton/detritus to IOPs
- **Smith & Baker (1981)**: Pure water absorption coefficients
- **Morel (1974)**: Pure water scattering coefficients

### Design Pattern
The implementation mirrors the existing `CanopySurface` pattern: the ocean (surface interface + water column + ocean floor) is encapsulated as an `OceanSurface <: AbstractSurfaceType` that internally runs adding-doubling through ocean sub-layers before presenting an effective surface reflectance/source to the atmospheric RT solver.

---

## 2. Architecture Overview

```
TOA ─────────────────────────────────── τ = 0
  │  Atmosphere (existing vSmartMOM RT)
  │  Layers: Rayleigh + aerosol + gas absorption
  │  Quadrature: N_atm points (standard Gauss-Lobatto)
BOA ─────────────────────────────────── τ_atm
  │
  ├─ Fresnel Interface Layer ←──────── NEW: maps between grids
  │    R_AA [N_atm × N_atm]           atmosphere-side Fresnel reflection
  │    R_OO [N_ocean × N_ocean]        ocean-side reflection (includes TIR)
  │    T_AO [N_ocean × N_atm]          atmosphere→ocean transmission
  │    T_OA [N_atm × N_ocean]          ocean→atmosphere transmission
  │    Quadrature: N_ocean = N_atm + N_tir points
  │
  ├─ Ocean Sub-Layers ←─────────────── NEW: adding-doubling on ocean grid
  │    Layer 1: z=0..5m
  │    Layer 2: z=5..10m
  │    ...
  │    Each layer has: τ, ω₀, phase function from constituent IOPs
  │    Optional: fluorescence source terms (isotropic, m=0 only)
  │    Optional: APAR accumulation
  │
  └─ Ocean Floor ←──────────────────── NEW: Lambertian reflector (lower BC)
       R_B = 2ρ·1·M·C for m=0          (Fell Eq. 3.5)
       R_B = 0 for m>0
```

### Key Technical Challenge: Quadrature Grid Mismatch

The ocean operates on an **extended quadrature grid** due to Snell's law:
- Atmospheric grid: N_atm points in μ ∈ [0, 1]
- Ocean grid: N_atm refracted points + N_tir TIR-zone points

The Fresnel interface matrices are **rectangular**, mapping between grids. The final `AddedLayer` returned to the atmosphere is on the atmospheric grid. This is the most mathematically subtle part of the implementation.

---

## 3. Part 1: OceanOptics.jl Package (START HERE)

Create a new Julia package (like `CanopyOptics.jl`) for Inherent Optical Properties (IOPs). This keeps ocean physics models reusable outside vSmartMOM.

### 3.1 Package Setup

```bash
# Create package skeleton
cd ~/code/gitHub
julia -e 'using Pkg; Pkg.generate("OceanOptics")'
cd OceanOptics
# Add to git, set up CI, etc.
```

### 3.2 Package Structure

```
OceanOptics.jl/
├── Project.toml
├── src/
│   ├── OceanOptics.jl              # Module entry point
│   ├── types.jl                    # All type definitions
│   ├── pure_water.jl               # Pure seawater absorption & scattering
│   ├── cdom.jl                     # CDOM/Gelbstoff (colored dissolved organic matter)
│   ├── phytoplankton.jl            # Phytoplankton absorption & scattering (Gordon 1992)
│   ├── detritus.jl                 # Detrital absorption & scattering
│   ├── inorganic_matter.jl         # Inorganic suspended particles
│   ├── pigments.jl                 # Additional pigment absorption spectra
│   ├── phase_functions.jl          # Oceanic phase functions (Petzold, Fournier-Forand)
│   ├── fluorescence.jl             # Chlorophyll fluorescence model
│   ├── iop_mixing.jl               # Combine constituents → total IOPs
│   └── water_refraction.jl         # Refractive index (Segelstein + T/S corrections)
├── data/
│   ├── smith_baker_1981.csv        # Pure water absorption [λ(nm), a(1/m)]
│   ├── pope_fry_1997.csv           # Updated pure water absorption (alternative)
│   ├── petzold_sdh.csv             # Petzold "San Diego Harbor" phase function
│   ├── phyto_absorption_normalized.csv   # a_PHY(λ)/a_PHY(440) from Gordon (1992)
│   └── detritus_absorption_normalized.csv # a_DET(λ)/a_DET(440) from Gordon (1992)
└── test/
    ├── runtests.jl
    ├── test_pure_water.jl
    ├── test_cdom.jl
    ├── test_phytoplankton.jl
    ├── test_phase_functions.jl
    └── test_fluorescence.jl
```

### 3.3 Module Entry Point: `src/OceanOptics.jl`

```julia
module OceanOptics

using LinearAlgebra
using Interpolations  # For spectral data lookup tables

# Type definitions
include("types.jl")

# Constituent IOPs
include("pure_water.jl")
include("cdom.jl")
include("phytoplankton.jl")
include("detritus.jl")
include("inorganic_matter.jl")
include("pigments.jl")

# Phase functions
include("phase_functions.jl")

# Fluorescence
include("fluorescence.jl")

# IOP mixing
include("iop_mixing.jl")

# Water refractive index
include("water_refraction.jl")

# Exports
export AbstractOceanConstituent, PureWater, CDOM, Phytoplankton, InorganicMatter
export ChlorophyllFluorescence, PhycocyaninFluorescence, PhycoerythrinFluorescence
export OceanLayer, OceanProfile
export PetzoldPhaseFunction, FournierForandPhaseFunction
export compute_absorption, compute_scattering, compute_iops
export compute_phase_function, phase_function_moments
export fluorescence_emission, fluorescence_activation
export water_refractive_index

end # module
```

### 3.4 Type Definitions: `src/types.jl`

```julia
# ============================================================================
# Abstract base type for ocean water constituents
# ============================================================================
abstract type AbstractOceanConstituent end

# ============================================================================
# Pure seawater
# ============================================================================
"""
    PureWater{FT}

Pure seawater optical properties.
- Absorption: Smith & Baker (1981) or Pope & Fry (1997) tabulated data
- Scattering: Einstein-Smoluchowski theory (Morel 1974)
    b_H2O(λ) = 0.00193 * (550/λ)^4.32  [1/m, λ in nm]
    p(cosθ) = 0.06225 * (1 + 0.835 cos²θ)
- Scattering is nearly Rayleigh-like with depolarization ratio ρ=0.09
"""
struct PureWater{FT} <: AbstractOceanConstituent
    "Use Pope & Fry (1997) instead of Smith & Baker (1981)"
    use_pope_fry::Bool
    "Temperature correction factor [optional, Pegau & Zanefeld 1993]"
    temperature::Union{Nothing, FT}    # °C
    "Salinity correction factor [optional, Buiteveld et al. 1994]"
    salinity::Union{Nothing, FT}       # PSU

    function PureWater{FT}(; use_pope_fry=false, temperature=nothing, salinity=nothing) where FT
        new{FT}(use_pope_fry, temperature, salinity)
    end
end
PureWater(; kwargs...) = PureWater{Float64}(; kwargs...)

# ============================================================================
# CDOM (Colored Dissolved Organic Matter) / Gelbstoff / Yellow Substance
# ============================================================================
"""
    CDOM{FT}

CDOM is a pure absorber (no scattering). Spectral absorption follows an
exponential decay (Bricaud et al. 1981):

    a_GS(λ) = a_ref * exp(-slope * (λ - λ_ref))

Typical values:
- Open ocean: a_GS(440) ~ 0.01-0.1 m⁻¹
- Coastal: a_GS(440) ~ 0.1-1.0 m⁻¹
- Rivers/lakes: a_GS(440) ~ 1-20 m⁻¹
- Slope q ~ 0.014 nm⁻¹ (range 0.01-0.02)
"""
struct CDOM{FT} <: AbstractOceanConstituent
    "Absorption coefficient at reference wavelength [1/m]"
    a_ref::FT
    "Reference wavelength [nm] (typically 440)"
    lambda_ref::FT
    "Exponential slope [1/nm] (typically 0.014)"
    slope::FT

    function CDOM{FT}(; a_ref, lambda_ref=FT(440), slope=FT(0.014)) where FT
        new{FT}(a_ref, lambda_ref, slope)
    end
end

# ============================================================================
# Phytoplankton
# ============================================================================
"""
    Phytoplankton{FT}

Phytoplankton optical properties using the Gordon (1992) bio-optical model.

The model separates organic marine particles into phytoplankton and detritus,
with the ratio depending on chlorophyll concentration:

    f(<CHL>) = 0.5 + 0.25 * log10(<CHL>)    # Eq. 5.7

Absorption (Gordon 1992):
    a_P(λ) = [a*_PHY(440) * a_PHY(λ) + a*_DET(440) * a_DET(λ) * (1-f)/f] * <CHL>
    a*_PHY(440) = 0.034 m²/mg,  a*_DET(440) = 0.030 m²/mg

Scattering (Gordon 1992):
    b_P(λ) = b_P(550) * <CHL>^0.62 * {f + (1-f) * 550/λ}
    b_P(550) = 0.30 m²/mg  (for <CHL> > 0.1 mg/m³)

Phase function: Petzold "San Diego Harbor" with modified forward peak.
Backscatter ratios: b_b,PHY = 0.0025, b_b,DET = 0.018
"""
struct Phytoplankton{FT} <: AbstractOceanConstituent
    "Chlorophyll-a concentration [mg/m³]"
    Chl::FT
    "Specific absorption coefficient at 440nm [m²/mg] (default 0.034)"
    a_star_phy_440::FT
    "Specific detrital absorption at 440nm [m²/mg] (default 0.030)"
    a_star_det_440::FT
    "Scattering coefficient at 550nm [m⁻¹/(mg/m³)^0.62] (default 0.30)"
    b_P_550::FT
    "Phytoplankton backscatter ratio (default 0.0025)"
    bb_ratio_phy::FT
    "Detritus backscatter ratio (default 0.018)"
    bb_ratio_det::FT

    function Phytoplankton{FT}(; Chl,
                                 a_star_phy_440=FT(0.034),
                                 a_star_det_440=FT(0.030),
                                 b_P_550=FT(0.30),
                                 bb_ratio_phy=FT(0.0025),
                                 bb_ratio_det=FT(0.018)) where FT
        new{FT}(Chl, a_star_phy_440, a_star_det_440, b_P_550, bb_ratio_phy, bb_ratio_det)
    end
end

# ============================================================================
# Inorganic suspended matter
# ============================================================================
"""
    InorganicMatter{FT}

Inorganic (mineral) suspended particles. Optical properties are poorly known;
typical parameterization from Fell (1997):

    c_AS(λ) = c_AS(440) * 440/λ     # Extinction [1/m]
    ω₀ = 0.95                        # Single-scattering albedo (weakly absorbing)

Phase function: modified Petzold with b_b = 0.036
"""
struct InorganicMatter{FT} <: AbstractOceanConstituent
    "Extinction coefficient at 440 nm [1/m]"
    c_ref::FT
    "Single-scattering albedo (default 0.95)"
    ssa::FT
    "Backscatter ratio (default 0.036)"
    bb_ratio::FT

    function InorganicMatter{FT}(; c_ref, ssa=FT(0.95), bb_ratio=FT(0.036)) where FT
        new{FT}(c_ref, ssa, bb_ratio)
    end
end

# ============================================================================
# Additional pigments (beyond Chl-a)
# ============================================================================
"""
    AbstractPigment

Pigments other than chlorophyll-a that contribute to absorption and
possibly fluorescence in the water column.
"""
abstract type AbstractPigment end

"""
    PhycocyaninPigment{FT} <: AbstractPigment

Phycocyanin — blue pigment found in cyanobacteria.
- Absorption peak: ~620 nm (broad, ~590-640 nm)
- Fluorescence: ~650 nm (if enabled)
- Typical concentrations: 0.1-100 μg/L in cyanobacteria blooms
"""
struct PhycocyaninPigment{FT} <: AbstractPigment
    "Concentration [μg/L]"
    concentration::FT
    "Specific absorption at 620 nm [m²/mg] (default ~0.007)"
    a_star_620::FT
end

"""
    PhycoerythrinPigment{FT} <: AbstractPigment

Phycoerythrin — red pigment found in some cyanobacteria and red algae.
- Absorption peak: ~565 nm (broad, ~490-580 nm)
- Fluorescence: ~580 nm (if enabled)
"""
struct PhycoerythrinPigment{FT} <: AbstractPigment
    "Concentration [μg/L]"
    concentration::FT
    "Specific absorption at 565 nm [m²/mg]"
    a_star_565::FT
end

"""
    FucoxanthinPigment{FT} <: AbstractPigment

Fucoxanthin — carotenoid pigment dominant in diatoms and brown algae.
- Absorption: broad 400-550 nm
- No fluorescence (energy transferred to Chl-a)
- Its presence indicates diatom-dominated phytoplankton community
"""
struct FucoxanthinPigment{FT} <: AbstractPigment
    "Concentration [μg/L]"
    concentration::FT
end

# ============================================================================
# Fluorescence types
# ============================================================================
"""
    AbstractFluorophore

Base type for fluorescent materials in the water column.
"""
abstract type AbstractFluorophore end

"""
    ChlorophyllFluorescence{FT} <: AbstractFluorophore

Chlorophyll-a fluorescence model following Gordon (1979) / Fell (1997).

Spectral emission is a Gaussian centered at 685 nm:
    g^C(λ) = (1/√(2π)σ) * exp(-(λ - λ₀)²/(2σ²))

The fluorescence quantum yield φ^C represents the fraction of absorbed
photons re-emitted as fluorescence. It varies with phytoplankton physiology:
- Healthy, low-light adapted: φ ~ 0.01-0.05
- High-light, nutrient replete: φ ~ 0.005-0.01
- Stressed (nutrient-limited): φ ~ 0.02-0.05
- Typical average: φ ~ 0.003 (Fischer et al. 1986)

The activation function (Fell Eq. 2.96) integrates absorbed PAR:
    A^C(τ) = ∫_Λ' a^C(τ; λ') · E°(τ; λ') · λ' dλ'
where E° is scalar irradiance and the integral is over excitation wavelengths.
"""
struct ChlorophyllFluorescence{FT} <: AbstractFluorophore
    "Fluorescence quantum yield [-] (default 0.003)"
    quantum_yield::FT
    "Emission center wavelength [nm] (default 685)"
    center_wavelength::FT
    "Emission Gaussian width σ [nm] (default 10.6, giving FWHM ≈ 25 nm)"
    sigma::FT
    "Excitation wavelength range [nm] (default (400, 700))"
    excitation_range::Tuple{FT, FT}

    function ChlorophyllFluorescence{FT}(;
            quantum_yield=FT(0.003),
            center_wavelength=FT(685),
            sigma=FT(10.6),
            excitation_range=(FT(400), FT(700))) where FT
        new{FT}(quantum_yield, center_wavelength, sigma, excitation_range)
    end
end

"""
    PhycocyaninFluorescence{FT} <: AbstractFluorophore

Phycocyanin fluorescence: emission at ~650 nm, excited by ~600-630 nm.
"""
struct PhycocyaninFluorescence{FT} <: AbstractFluorophore
    quantum_yield::FT
    center_wavelength::FT     # ~650 nm
    sigma::FT                 # ~15 nm
    excitation_range::Tuple{FT, FT}  # (590, 640)

    function PhycocyaninFluorescence{FT}(;
            quantum_yield=FT(0.01),
            center_wavelength=FT(650),
            sigma=FT(15),
            excitation_range=(FT(590), FT(640))) where FT
        new{FT}(quantum_yield, center_wavelength, sigma, excitation_range)
    end
end

"""
    PhycoerythrinFluorescence{FT} <: AbstractFluorophore

Phycoerythrin fluorescence: emission at ~580 nm, excited by ~490-570 nm.
"""
struct PhycoerythrinFluorescence{FT} <: AbstractFluorophore
    quantum_yield::FT
    center_wavelength::FT     # ~580 nm
    sigma::FT                 # ~12 nm
    excitation_range::Tuple{FT, FT}  # (490, 570)

    function PhycoerythrinFluorescence{FT}(;
            quantum_yield=FT(0.01),
            center_wavelength=FT(580),
            sigma=FT(12),
            excitation_range=(FT(490), FT(570))) where FT
        new{FT}(quantum_yield, center_wavelength, sigma, excitation_range)
    end
end

# ============================================================================
# Phase function types
# ============================================================================
"""
    AbstractOceanPhaseFunction

Base type for oceanic phase functions.
"""
abstract type AbstractOceanPhaseFunction end

"""
    PetzoldPhaseFunction{FT} <: AbstractOceanPhaseFunction

Petzold (1972) "San Diego Harbor" volume scattering function.
The most commonly used oceanic particle phase function. Highly forward-peaked
(backscatter ratio b_b/b ~ 0.018 for the average particle mixture).

The forward peak (θ < 0.1°) is modified following Fell (1997) to achieve
a specified backscatter ratio for different particle types.
"""
struct PetzoldPhaseFunction{FT} <: AbstractOceanPhaseFunction
    "Target backscatter ratio b_b/b"
    backscatter_ratio::FT
    "Maximum angle for forward peak modification [degrees] (default 0.1)"
    theta_forward_max::FT

    function PetzoldPhaseFunction{FT}(; backscatter_ratio=FT(0.018),
                                        theta_forward_max=FT(0.1)) where FT
        new{FT}(backscatter_ratio, theta_forward_max)
    end
end

"""
    FournierForandPhaseFunction{FT} <: AbstractOceanPhaseFunction

Fournier-Forand analytical phase function (Fournier & Forand 1994,
Fournier & Jonasz 1999). Based on anomalous diffraction by a
Junge-type particle size distribution. Widely used in modern ocean
optics as it has an analytical form parameterized by the backscatter ratio.

This is preferred over Petzold for most applications because the backscatter
ratio directly determines the full phase function shape.
"""
struct FournierForandPhaseFunction{FT} <: AbstractOceanPhaseFunction
    "Backscatter ratio b_b/b"
    backscatter_ratio::FT

    function FournierForandPhaseFunction{FT}(; backscatter_ratio=FT(0.018)) where FT
        new{FT}(backscatter_ratio)
    end
end

# ============================================================================
# Ocean profile
# ============================================================================
"""
    OceanLayer{FT}

A single vertically-homogeneous ocean layer with specified constituents.
"""
struct OceanLayer{FT}
    "Depth of layer top [m] (0 = surface)"
    depth_top::FT
    "Depth of layer bottom [m]"
    depth_bottom::FT
    "Water constituents in this layer"
    constituents::Vector{AbstractOceanConstituent}
    "Additional pigments (beyond those in constituents)"
    pigments::Vector{AbstractPigment}
    "Fluorophores active in this layer"
    fluorophores::Vector{AbstractFluorophore}
    "Phase function model for particles (default: Petzold)"
    phase_function::AbstractOceanPhaseFunction
end

"""
    OceanProfile{FT}

Vertical profile of the ocean, consisting of multiple layers from surface
to a specified depth. The profile is specified from surface downward (z=0
at surface, increasing with depth).

Recommended vertical resolution (Fell 1997 Section 5.2.1):
- 0-10m: 1m layers (fine resolution near surface)
- 10-20m: 2m layers
- 20-60m: 5m layers
- Below 60m: can be coarser or omitted (negligible signal at surface)
"""
struct OceanProfile{FT}
    "Ocean layers, ordered from surface downward"
    layers::Vector{OceanLayer{FT}}

    function OceanProfile(layers::Vector{OceanLayer{FT}}) where FT
        # Validate: layers should be contiguous, non-overlapping, surface-downward
        for i in 2:length(layers)
            @assert layers[i].depth_top ≈ layers[i-1].depth_bottom "Layer $i top must match layer $(i-1) bottom"
        end
        @assert layers[1].depth_top ≈ zero(FT) "First layer must start at surface (depth_top=0)"
        new{FT}(layers)
    end
end

"""Total depth of the ocean profile [m]"""
total_depth(p::OceanProfile) = p.layers[end].depth_bottom

"""Number of layers"""
n_layers(p::OceanProfile) = length(p.layers)

"""Geometric thickness of a layer [m]"""
layer_thickness(layer::OceanLayer) = layer.depth_bottom - layer.depth_top
```

### 3.5 Pure Water: `src/pure_water.jl`

```julia
# ============================================================================
# Pure seawater absorption
# ============================================================================

"""
Built-in Smith & Baker (1981) pure water absorption coefficients.
Table spans 200-800 nm in 10 nm steps.
"""
const SMITH_BAKER_1981 = Dict{Int, Float64}(
    # λ(nm) => a(1/m)
    200 => 3.07, 210 => 1.99, 220 => 1.31, 230 => 0.927, 240 => 0.720,
    250 => 0.559, 260 => 0.457, 270 => 0.373, 280 => 0.288, 290 => 0.215,
    300 => 0.141, 310 => 0.105, 320 => 0.0844, 330 => 0.0678, 340 => 0.0561,
    350 => 0.0463, 360 => 0.0379, 370 => 0.0300, 380 => 0.0220, 390 => 0.0191,
    400 => 0.0171, 410 => 0.0162, 420 => 0.0153, 430 => 0.0144, 440 => 0.0145,
    450 => 0.0145, 460 => 0.0156, 470 => 0.0156, 480 => 0.0176, 490 => 0.0196,
    500 => 0.0257, 510 => 0.0357, 520 => 0.0477, 530 => 0.0507, 540 => 0.0558,
    550 => 0.0638, 560 => 0.0708, 570 => 0.0799, 580 => 0.108, 590 => 0.157,
    600 => 0.244, 610 => 0.289, 620 => 0.309, 630 => 0.319, 640 => 0.329,
    650 => 0.349, 660 => 0.400, 670 => 0.430, 680 => 0.450, 690 => 0.500,
    700 => 0.650, 710 => 0.839, 720 => 1.169, 730 => 1.799, 740 => 2.38,
    750 => 2.47, 760 => 2.55, 770 => 2.51, 780 => 2.36, 790 => 2.16, 800 => 2.07
)

"""
    absorption(w::PureWater, lambda_nm)

Pure water absorption coefficient [1/m] at wavelength lambda_nm [nm].
Interpolates linearly between tabulated values.
"""
function compute_absorption(w::PureWater{FT}, lambda_nm::FT) where FT
    # Get sorted table keys for interpolation
    lambdas = sort(collect(keys(SMITH_BAKER_1981)))
    a_vals = [FT(SMITH_BAKER_1981[l]) for l in lambdas]

    # Linear interpolation
    if lambda_nm <= FT(lambdas[1])
        return a_vals[1]
    elseif lambda_nm >= FT(lambdas[end])
        return a_vals[end]
    end

    # Find bracketing indices
    idx = searchsortedlast(lambdas, lambda_nm)
    idx = clamp(idx, 1, length(lambdas)-1)
    λ₁, λ₂ = FT(lambdas[idx]), FT(lambdas[idx+1])
    a₁, a₂ = a_vals[idx], a_vals[idx+1]

    # Linear interpolation
    frac = (lambda_nm - λ₁) / (λ₂ - λ₁)
    return a₁ + frac * (a₂ - a₁)
end

# ============================================================================
# Pure seawater scattering (Einstein-Smoluchowski / Morel 1974)
# ============================================================================

"""
    compute_scattering(w::PureWater, lambda_nm)

Pure water scattering coefficient [1/m] at wavelength lambda_nm [nm].
Morel (1974):
    b_H2O(λ) = 0.00193 * (550/λ)^4.32

Note the exponent 4.32 (not exactly 4) accounts for the wavelength
dependence of the depolarization ratio.
"""
function compute_scattering(w::PureWater{FT}, lambda_nm::FT) where FT
    FT(0.00193) * (FT(550) / lambda_nm)^FT(4.32)
end

"""
    pure_water_phase_function(cos_theta)

Pure water phase function (Morel 1974):
    p(cosθ) = 0.06225 * (1 + 0.835 * cos²θ)

This is a modified Rayleigh phase function with depolarization ratio ρ=0.09.
The factor 0.06225 ensures ∫₄π p dΩ = 1.

Note: for polarized RT, use the full depolarized Rayleigh matrix with ρ=0.09.
"""
function pure_water_phase_function(::Type{FT}, cos_theta::FT) where FT
    FT(0.06225) * (one(FT) + FT(0.835) * cos_theta^2)
end

"""
    compute_iops(w::PureWater, lambda_nm)

Total IOPs for pure water: (absorption, scattering, extinction, ssa)
"""
function compute_iops(w::PureWater{FT}, lambda_nm::FT) where FT
    a = compute_absorption(w, lambda_nm)
    b = compute_scattering(w, lambda_nm)
    c = a + b
    ω = b / c
    return (; a, b, c, ω)
end
```

### 3.6 CDOM: `src/cdom.jl`

```julia
"""
    compute_absorption(cdom::CDOM, lambda_nm)

CDOM absorption [1/m] following Bricaud et al. (1981):
    a_GS(λ) = a_ref * exp(-slope * (λ - λ_ref))

CDOM is a pure absorber: scattering is negligible.
"""
function compute_absorption(cdom::CDOM{FT}, lambda_nm::FT) where FT
    cdom.a_ref * exp(-cdom.slope * (lambda_nm - cdom.lambda_ref))
end

compute_scattering(::CDOM{FT}, ::FT) where FT = zero(FT)

function compute_iops(cdom::CDOM{FT}, lambda_nm::FT) where FT
    a = compute_absorption(cdom, lambda_nm)
    b = zero(FT)
    c = a
    ω = zero(FT)
    return (; a, b, c, ω)
end
```

### 3.7 Phytoplankton: `src/phytoplankton.jl`

```julia
"""
Normalized phytoplankton absorption spectrum a_PHY(λ)/a_PHY(440).
From Gordon (1992), digitized from measurements.
"""
const PHYTO_ABSORPTION_NORM = Dict{Int, Float64}(
    # λ(nm) => a_PHY(λ)/a_PHY(440)
    400 => 0.87, 410 => 0.93, 420 => 0.97, 430 => 1.00, 440 => 1.00,
    450 => 0.94, 460 => 0.87, 470 => 0.76, 480 => 0.64, 490 => 0.52,
    500 => 0.41, 510 => 0.33, 520 => 0.27, 530 => 0.24, 540 => 0.22,
    550 => 0.21, 560 => 0.20, 570 => 0.19, 580 => 0.18, 590 => 0.17,
    600 => 0.16, 610 => 0.15, 620 => 0.14, 630 => 0.13, 640 => 0.12,
    650 => 0.11, 660 => 0.13, 670 => 0.33, 680 => 0.26, 690 => 0.09,
    700 => 0.03
)

"""
Normalized detrital absorption spectrum a_DET(λ)/a_DET(440).
From Gordon (1992).
"""
const DETRITUS_ABSORPTION_NORM = Dict{Int, Float64}(
    400 => 1.50, 410 => 1.24, 420 => 1.11, 430 => 1.04, 440 => 1.00,
    450 => 0.96, 460 => 0.91, 470 => 0.86, 480 => 0.82, 490 => 0.77,
    500 => 0.72, 510 => 0.67, 520 => 0.63, 530 => 0.59, 540 => 0.55,
    550 => 0.52, 560 => 0.49, 570 => 0.46, 580 => 0.44, 590 => 0.41,
    600 => 0.39, 610 => 0.37, 620 => 0.35, 630 => 0.34, 640 => 0.32,
    650 => 0.31, 660 => 0.30, 670 => 0.29, 680 => 0.28, 690 => 0.27,
    700 => 0.26
)

"""
    phytoplankton_fraction(Chl)

Fraction of total organic particles that are phytoplankton (Gordon 1992):
    f(<CHL>) = 0.5 + 0.25 * log10(<CHL>)

Defined for 0.01 ≤ <CHL> ≤ 100 mg/m³. Clamped to [0.01, 0.99].
"""
function phytoplankton_fraction(Chl::FT) where FT
    f = FT(0.5) + FT(0.25) * log10(max(Chl, FT(0.01)))
    return clamp(f, FT(0.01), FT(0.99))
end

"""
    compute_absorption(p::Phytoplankton, lambda_nm)

Particle absorption coefficient [1/m] from Gordon (1992) bio-optical model:
    a_P(λ) = [a*_PHY(440)·a_PHY(λ) + a*_DET(440)·a_DET(λ)·(1-f)/f] · <CHL>

where f = phytoplankton_fraction(<CHL>).
"""
function compute_absorption(p::Phytoplankton{FT}, lambda_nm::FT) where FT
    f = phytoplankton_fraction(p.Chl)

    # Interpolate normalized absorption spectra
    a_phy_norm = _interp_spectrum(PHYTO_ABSORPTION_NORM, lambda_nm, FT)
    a_det_norm = _interp_spectrum(DETRITUS_ABSORPTION_NORM, lambda_nm, FT)

    a_P = (p.a_star_phy_440 * a_phy_norm + p.a_star_det_440 * a_det_norm * (one(FT) - f) / f) * p.Chl
    return a_P
end

"""
    compute_scattering(p::Phytoplankton, lambda_nm)

Particle scattering coefficient [1/m] from Gordon (1992):
    b_P(λ) = b_P(550) · <CHL>^0.62 · {f·<CHL> + [1 - f·<CHL>]·550/λ}

Phytoplankton particles have nearly spectrally flat scattering,
while smaller detritus particles scatter proportional to 1/λ.
"""
function compute_scattering(p::Phytoplankton{FT}, lambda_nm::FT) where FT
    f = phytoplankton_fraction(p.Chl)
    b = p.b_P_550 * p.Chl^FT(0.62) * (f + (one(FT) - f) * FT(550) / lambda_nm)
    return b
end

function compute_iops(p::Phytoplankton{FT}, lambda_nm::FT) where FT
    a = compute_absorption(p, lambda_nm)
    b = compute_scattering(p, lambda_nm)
    c = a + b
    ω = c > zero(FT) ? b / c : zero(FT)
    return (; a, b, c, ω)
end

# Helper: linear interpolation of Dict{Int,Float64} spectrum
function _interp_spectrum(table::Dict{Int,Float64}, lambda_nm::FT, ::Type{FT}) where FT
    lambdas = sort(collect(keys(table)))
    if lambda_nm <= FT(lambdas[1])
        return FT(table[lambdas[1]])
    elseif lambda_nm >= FT(lambdas[end])
        return FT(table[lambdas[end]])
    end
    idx = searchsortedlast(lambdas, lambda_nm)
    idx = clamp(idx, 1, length(lambdas)-1)
    λ₁, λ₂ = FT(lambdas[idx]), FT(lambdas[idx+1])
    v₁, v₂ = FT(table[lambdas[idx]]), FT(table[lambdas[idx+1]])
    frac = (lambda_nm - λ₁) / (λ₂ - λ₁)
    return v₁ + frac * (v₂ - v₁)
end
```

### 3.8 Inorganic Matter: `src/inorganic_matter.jl`

```julia
"""
    compute_absorption(im::InorganicMatter, lambda_nm)
    compute_scattering(im::InorganicMatter, lambda_nm)

Inorganic suspended matter IOPs:
    c(λ) = c_ref * (440/λ)
    a(λ) = (1 - ω₀) * c(λ)
    b(λ) = ω₀ * c(λ)
"""
function compute_absorption(im::InorganicMatter{FT}, lambda_nm::FT) where FT
    c = im.c_ref * (FT(440) / lambda_nm)
    return (one(FT) - im.ssa) * c
end

function compute_scattering(im::InorganicMatter{FT}, lambda_nm::FT) where FT
    c = im.c_ref * (FT(440) / lambda_nm)
    return im.ssa * c
end

function compute_iops(im::InorganicMatter{FT}, lambda_nm::FT) where FT
    c = im.c_ref * (FT(440) / lambda_nm)
    a = (one(FT) - im.ssa) * c
    b = im.ssa * c
    ω = im.ssa
    return (; a, b, c, ω)
end
```

### 3.9 Phase Functions: `src/phase_functions.jl`

```julia
"""
    fournier_forand(cos_theta, backscatter_ratio)

Fournier-Forand phase function. Analytical form parameterized by
the backscatter ratio b_b/b. This is the preferred phase function
for ocean RT as it has a closed form that matches measured data well.

Reference: Fournier & Forand (1994), Fournier & Jonasz (1999).
"""
function fournier_forand(::Type{FT}, cos_theta::FT, bb_ratio::FT) where FT
    # Compute refractive index and Junge slope from bb_ratio
    # Using the Fournier-Forand parameterization
    # ... (implementation details from Mobley et al. 2002)
    # This is a multi-line formula; see Ocean Optics Web Book for details
    error("Not yet implemented — see oceanopticsbook.info for the full formula")
end

"""
    petzold_sdh(cos_theta)

Petzold (1972) "San Diego Harbor" average particle VSF, normalized as a phase function.
Tabulated values from the original measurements. Highly forward-peaked.
"""
function petzold_sdh(::Type{FT}, cos_theta::FT) where FT
    # Tabulated Petzold phase function
    # ... (from data/petzold_sdh.csv or built-in)
    error("Not yet implemented — load from tabulated data")
end

"""
    phase_function_legendre_moments(pf::AbstractOceanPhaseFunction, N_moments)

Compute Legendre polynomial expansion coefficients ωₗ of a phase function:
    p(cosθ) = Σₗ (2l+1) ωₗ Pₗ(cosθ)

These are used to build Greek coefficients / Z-matrices for the RT solver.
"""
function phase_function_legendre_moments(pf::AbstractOceanPhaseFunction, N_moments::Int)
    # Numerical integration of ωₗ = ∫₋₁¹ p(cosθ) Pₗ(cosθ) d(cosθ)
    # using Gauss-Legendre quadrature
    error("Not yet implemented")
end
```

### 3.10 Fluorescence: `src/fluorescence.jl`

```julia
"""
    fluorescence_emission(f::ChlorophyllFluorescence, lambda_nm)

Chlorophyll-a fluorescence spectral emission function g^C(λ) [1/nm].
Gaussian centered at λ₀=685nm with σ=10.6nm (FWHM ≈ 25nm).

Normalized: ∫ g^C(λ) dλ = 1 over the fluorescence emission band.
"""
function fluorescence_emission(f::ChlorophyllFluorescence{FT}, lambda_nm::FT) where FT
    g = exp(-(lambda_nm - f.center_wavelength)^2 / (2 * f.sigma^2)) / (sqrt(2 * FT(π)) * f.sigma)
    return g
end

# Same for other fluorophores
function fluorescence_emission(f::PhycocyaninFluorescence{FT}, lambda_nm::FT) where FT
    g = exp(-(lambda_nm - f.center_wavelength)^2 / (2 * f.sigma^2)) / (sqrt(2 * FT(π)) * f.sigma)
    return g
end

function fluorescence_emission(f::PhycoerythrinFluorescence{FT}, lambda_nm::FT) where FT
    g = exp(-(lambda_nm - f.center_wavelength)^2 / (2 * f.sigma^2)) / (sqrt(2 * FT(π)) * f.sigma)
    return g
end

"""
    fluorescence_activation(f::ChlorophyllFluorescence, a_chl_spectrum, E_scalar_spectrum,
                            lambda_grid)

Compute the chlorophyll fluorescence activation function A^C (Fell Eq. 2.96):

    A^C = ∫_{Λ'} a^C(λ') · E°(λ') · λ' dλ'

Arguments:
- `a_chl_spectrum`: Chlorophyll absorption [1/m] at each wavelength in `lambda_grid`
- `E_scalar_spectrum`: Scalar irradiance [W/m²/nm] at each wavelength
- `lambda_grid`: Wavelength grid [nm] covering excitation range

Returns: A^C in [W·nm/m]
"""
function fluorescence_activation(f::ChlorophyllFluorescence{FT},
                                  a_chl_spectrum::Vector{FT},
                                  E_scalar_spectrum::Vector{FT},
                                  lambda_grid::Vector{FT}) where FT
    # Filter to excitation range
    λ_min, λ_max = f.excitation_range
    activation = zero(FT)

    for i in 1:length(lambda_grid)-1
        λ = lambda_grid[i]
        if λ_min ≤ λ ≤ λ_max
            dλ = lambda_grid[i+1] - lambda_grid[i]
            activation += a_chl_spectrum[i] * E_scalar_spectrum[i] * λ * dλ
        end
    end

    return activation
end

"""
    fluorescence_source(f::ChlorophyllFluorescence, lambda_nm, activation, c_total)

Fluorescence source function J^C(λ) [W/m²/sr/nm] (Fell Eq. 2.95):

    J^C(λ) = g^C(λ) · φ^C / (4π · c_total · λ) · A^C

This is isotropic: emitted equally in all directions.
Only contributes to Fourier moment m=0.
"""
function fluorescence_source(f::ChlorophyllFluorescence{FT},
                              lambda_nm::FT,
                              activation::FT,
                              c_total::FT) where FT
    g = fluorescence_emission(f, lambda_nm)
    J = g * f.quantum_yield / (4 * FT(π) * c_total * lambda_nm) * activation
    return J
end
```

### 3.11 IOP Mixing: `src/iop_mixing.jl`

```julia
"""
    compute_total_iops(constituents::Vector{AbstractOceanConstituent}, lambda_nm)

Combine all constituent IOPs into total layer optical properties:
    a_total = Σ a_i
    b_total = Σ b_i
    c_total = a_total + b_total
    ω₀ = b_total / c_total
"""
function compute_total_iops(constituents::Vector{<:AbstractOceanConstituent},
                            lambda_nm::FT) where FT
    a_total = zero(FT)
    b_total = zero(FT)

    for constituent in constituents
        iops = compute_iops(constituent, lambda_nm)
        a_total += iops.a
        b_total += iops.b
    end

    c_total = a_total + b_total
    ω₀ = c_total > zero(FT) ? b_total / c_total : zero(FT)

    return (; a=a_total, b=b_total, c=c_total, ω=ω₀)
end

"""
    compute_layer_optical_depth(layer::OceanLayer, lambda_nm)

Optical depth of an ocean layer: τ = c_total · Δz
where Δz is the geometric thickness [m].
"""
function compute_layer_optical_depth(layer::OceanLayer{FT}, lambda_nm::FT) where FT
    iops = compute_total_iops(layer.constituents, lambda_nm)
    dz = layer_thickness(layer)
    τ = iops.c * dz
    return (; τ, ω=iops.ω, a=iops.a, b=iops.b, c=iops.c)
end
```

### 3.12 Water Refractive Index: `src/water_refraction.jl`

```julia
"""
    water_refractive_index(lambda_nm; temperature=20.0, salinity=35.0)

Real refractive index of seawater as a function of wavelength,
temperature, and salinity (Quan & Fry 1995).

For the imaginary part (absorption), use the absorption coefficient from
`compute_absorption(PureWater(), lambda_nm)` via:
    k = a · λ / (4π)

Default: T=20°C, S=35 PSU (open ocean average).
For air-water interface calculations, the relative index is n_water/n_air ≈ n_water.
"""
function water_refractive_index(::Type{FT}, lambda_nm::FT;
                                 temperature::FT=FT(20),
                                 salinity::FT=FT(35)) where FT
    # Quan & Fry (1995) empirical formula
    # Valid for 400-700 nm, 0-30°C, 0-35 PSU
    λ = lambda_nm  # nm
    T = temperature # °C
    S = salinity    # PSU

    n = FT(1.31405) +
        FT(1.779e-4) * T +
        FT(-1.05e-6) * T^2 +
        FT(1.6e-8) * T^3 +
        FT(-2.02e-6) * T * S +
        (FT(15.868) + FT(0.01155) * T + FT(-0.00423) * T^2 + FT(-4382) / λ^2) / λ +
        FT(-0.00423) * S +
        FT(0.01155) * S / λ

    # Simplified version for broader wavelength range:
    # n ≈ 1.333 at 550nm, varies ~1.329 (700nm) to ~1.345 (400nm)
    return n
end
```

---

## 4. Part 2: OceanSurface in vSmartMOM (DEFERRED)

### 4.1 Type Definition

Add to `src/CoreRT/types.jl` after `CanopySurface` (line ~571):

```julia
mutable struct OceanSurface{FT} <: AbstractSurfaceType
    # Surface interface
    surface_reflection::Union{CoxMunkSurface{FT}, Nothing}  # nothing = flat Fresnel
    n_water::Union{Nothing, Complex{FT}}  # override refractive index

    # Ocean body
    ocean_profile::Any      # OceanOptics.OceanProfile{FT}
    n_layers::Int
    depth::FT               # Total depth [m]

    # Ocean floor (lower BC, like canopy soil)
    floor::AbstractSurfaceType

    # Quadrature settings
    n_quad_tir::Int         # Extra TIR zone points (default 5)

    # Options
    rough_transmission::Bool    # Rough surface Fresnel transmission
    compute_apar::Bool          # Store APAR per layer

    # Cache (lazily initialized)
    _cache::Any             # Union{Nothing, OceanCache}
end
```

### 4.2 Quadrature Remapping

New file: `src/CoreRT/Surfaces/ocean_quadrature.jl`

Snell's law maps atmospheric μᵢ to ocean μ*ᵢ (Fell Eq. 3.32):
```
μ*ᵢ = √(1 - (n₁/n₂)²·(1 - μᵢ²))
```

Modified weights (Fell Eq. 3.33):
```
c*ᵢ = (μᵢ/μ*ᵢ)·(n₁/n₂)²·cᵢ
```

TIR zone: additional N' Gauss-Legendre points in [0, μ_C] where μ_C = √(1-(n₁/n₂)²).

Energy conservation: distribute weight residual ΔC* = 1 - μ_C - Σc*ᵢ (Fell Eq. 3.34).

### 4.3 Fresnel Interface Layer

New file: `src/CoreRT/Surfaces/fresnel_interface.jl`

Flat surface operators (all diagonal, Fell Eq. 3.35-3.38):
- `R^AA[i,i] = r_F(μᵢ, n_O/n_A)` — atmospheric Fresnel reflection
- `R^OO[i,i] = r_F(μ*ᵢ, n_A/n_O)` for refracted angles; `= 1` for TIR zone
- `T^AO[i,i] = t_F(μᵢ, n_O/n_A)·(n_O/n_A)²` — downward transmission
- `T^OA[i,i] = t_F(μ*ᵢ, n_A/n_O)·(n_A/n_O)²` — upward transmission

Note: T^AO and T^OA are **rectangular** matrices (N_ocean × N_atm and N_atm × N_ocean).

For polarized RT: extend existing `fresnel.jl` with transmission Mueller matrices.

### 4.4 create_surface_layer! for OceanSurface

File: `src/CoreRT/Surfaces/ocean_surface.jl`

```
1. Lazy cache init
   - Compute ocean quadrature (Snell remapping + TIR points)
   - Build Fresnel interface matrices
   - Allocate ocean-grid AddedLayer/CompositeLayer

2. Ocean sub-layer loop (on ocean quadrature grid)
   For each layer:
   a. Compute IOPs from OceanOptics
   b. Build CoreScatteringOpticalProperties(τ, ω, Z⁺⁺, Z⁻⁺)
   c. elemental! → doubling! → interaction!
   d. If fluorescence & m==0: add J_fluor to source terms
   e. If compute_apar: accumulate APAR

3. Ocean floor
   - Create floor AddedLayer (Lambertian: R_B = 2ρ·1·M·C for m=0)
   - interaction! with ocean composite

4. Fresnel interface coupling
   - R_eff = R_AA + T_OA · R_ocean · (I - R_OO · R_ocean)⁻¹ · T_AO
   - J_eff from interaction of Fresnel with ocean composite

5. Copy R_eff, J_eff into atmospheric-grid added_layer
```

### 4.5 YAML Configuration

```yaml
ocean:
  depth: 50.0
  n_layers: 10
  floor_albedo: 0.05
  n_quad_tir: 5
  flat_surface: false

  constituents:
    pure_water: true
    cdom:
      a_440: 0.1
      slope: 0.014
    phytoplankton:
      Chl: 0.5
    inorganic_matter:
      c_440: 0.0

  fluorescence:
    enabled: true
    quantum_yield: 0.003
```

### 4.6 Integration Points

| File | Change |
|------|--------|
| `src/CoreRT/types.jl` | Add `OceanSurface{FT}` struct |
| `src/CoreRT/CoreRT.jl` | Include new files, export OceanSurface |
| `src/CoreRT/rt_run.jl` | Pre-RT ocean cache init + surface dispatch (like canopy at lines 140-227) |
| `src/CoreRT/rt_run_lin.jl` | Linearized surface dispatch (Phase 4) |
| `src/IO/Parameters.jl` | `_parse_ocean_section()`, BRDF_MAP entry |
| `src/CoreRT/Surfaces/fresnel.jl` | Add Fresnel transmission Mueller matrices |
| `test/runtests.jl` | Add ocean test set |

---

## 5. Part 3: Fluorescence and APAR (DEFERRED)

### Fluorescence Implementation
- Two-pass approach (Fell Section 2.2.9):
  - Pass 1: Compute scalar irradiance E°(z, λ) at each depth
  - Pass 2: Compute activation A^C(z) and add fluorescence source J^C
- Simplified single-pass: approximate E° from Beer-Lambert attenuation of direct solar beam
- Isotropic emission (m=0 only): J adds to j₀⁺ and j₀⁻ equally

### APAR Computation
- Natural by-product: `APAR(z) = ∫₄₀₀⁷⁰⁰ a_phyto(z,λ)·E°(z,λ) dλ`
- Store per layer in OceanCache.apar

### Additional Fluorophores
- Phycocyanin: excitation ~620nm, emission ~650nm
- Phycoerythrin: excitation ~565nm, emission ~580nm
- Same mathematical framework, different (λ₀, σ, φ, excitation range)

---

## 6. Part 4: Linearization (DEFERRED)

### Retrieval Parameters
- Chlorophyll concentration (Chl)
- CDOM absorption (a_ref at 440nm)
- Ocean floor albedo
- Fluorescence quantum yield
- Wind speed (existing CoxMunk)

### Approach
- Finite-difference Jacobians (matching canopy pattern)
- `create_surface_layer!` linearized dispatch in `ocean_surface_lin.jl`
- ParameterLayout extension for ocean parameter indexing

---

## 7. Implementation Phases

### Phase 1: OceanOptics.jl package (START HERE)
- [x] Package skeleton
- [ ] types.jl — all type definitions
- [ ] pure_water.jl — Smith & Baker absorption + Morel scattering
- [ ] cdom.jl — exponential absorption model
- [ ] phytoplankton.jl — Gordon (1992) bio-optical model
- [ ] inorganic_matter.jl — mineral particle IOPs
- [ ] phase_functions.jl — Petzold and Fournier-Forand
- [ ] fluorescence.jl — Gaussian emission + activation
- [ ] iop_mixing.jl — combine constituents
- [ ] water_refraction.jl — refractive index
- [ ] Tests for all of the above

### Phase 2: OceanSurface in vSmartMOM (flat surface, pure water)
- [ ] OceanSurface type definition
- [ ] Snell's law quadrature remapping
- [ ] Flat Fresnel interface operators
- [ ] Single-layer ocean RT
- [ ] Lambertian ocean floor
- [ ] create_surface_layer! dispatch
- [ ] YAML parsing
- [ ] Energy conservation tests

### Phase 3: Realistic ocean + fluorescence
- [ ] Multi-layer ocean with constituent profiles
- [ ] Delta-M truncation for oceanic phase functions
- [ ] Chlorophyll fluorescence source terms
- [ ] APAR computation
- [ ] Rough-surface transmission (extend CoxMunk)
- [ ] Additional pigments

### Phase 4: Linearization
- [ ] Finite-difference Jacobians for ocean parameters
- [ ] ParameterLayout extension
- [ ] Jacobian validation tests

### Phase 5: Polish
- [ ] GPU support
- [ ] Performance optimization
- [ ] Validation against MOMO/HydroLight/OSOAA
- [ ] Documentation

---

## 8. Verification Plan

1. **Unit tests**: Fresnel R+T=1, Snell's law, TIR zone, IOP formulas vs published values
2. **Energy conservation**: R + T + A = 1 for full atmosphere-ocean system
3. **Analytical limits**: Pure absorbing ocean (R=0), conservative ocean (kD from diffusion theory)
4. **Fell thesis**: Reproduce MOMO validation problems 1-6 from Chapter 4
5. **Fluorescence**: Signal magnitude O(10⁻³) relative to reflected solar near 685nm
6. **Jacobians**: Analytic vs finite-difference with `rel_errors()` from test_helpers.jl
