"""
    GCHPAerosolOptics

Couple GEOS-Chem High Performance (GCHP / GOCART–TOMAS) aerosol output to vSmartMOM
radiative transfer. Reads a GCHP column, builds a per-layer, per-bin, speciated aerosol
state, and produces aerosol optical properties — from a quick column AOD up to the full
per-layer scattering optics (τ, ϖ, Greek coefficients) the RT solver consumes.

SCOPE: this is the I/O + optics-assembly layer. The actual Mie computation and RT live in
vSmartMOM. See `proposals/gchp_aerosol_optics/STUDENT_BRIEF.md`.

STATUS: scaffold. Most bodies are TODO; many already have a working reference on the
vSmartMOM `gchp-io` branch (cited inline) that should be ported here.
"""
module GCHPAerosolOptics

using vSmartMOM                # core: Scattering (Mie), CoreRT types
using NCDatasets
using Distributions           # LogNormal, fitting
using Interpolations          # RI(λ) interpolation, spline across bins
using StaticArrays
using LinearAlgebra           # SVD / EOF reduction
import YAML

# ── I/O: GCHP column → scene ─────────────────────────────────────────────────
include("io/gchp_scene.jl")

# ── Size-distribution representation (the Mode A / Mode B fork) ───────────────
include("representation/representation.jl")
include("representation/sectional.jl")
include("representation/lognormal_fit.jl")

# ── Speciation → effective refractive index ──────────────────────────────────
include("mixing/effective_ri.jl")

# ── Optics ───────────────────────────────────────────────────────────────────
include("optics/column_aod.jl")     # Stage 1 (exists on gchp-io): bins → AOD
include("optics/layer_optics.jl")   # Stage 3 (NEW): sectional → AerosolOptics

# ── ML training dataset assembly ─────────────────────────────────────────────
include("ml/scene_dataset.jl")
include("ml/reduction.jl")

export GCHPFile, GCHPScene, scene_at, scenes
export AerosolRepresentation, Sectional, LogNormalFit, SplineSectional
export effective_ri
export column_aod, layer_aerosol_optics
export SceneDataset, build_dataset, eof_reduce

end # module
