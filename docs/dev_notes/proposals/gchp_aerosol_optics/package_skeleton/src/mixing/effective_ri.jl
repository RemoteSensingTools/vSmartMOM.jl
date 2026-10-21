# Speciation → effective complex refractive index. This is what the aerosol *types* mostly
# feed into (brief §4, §6). PORT TARGET: `effective_ri` + the mixing-rule types in
# gchp-io:src/Aerosols/aod_diagnostics.jl and gchp-io:src/Aerosols/abstract_types.jl.

"""
    MixingRule

Effective-medium rule combining per-species complex RIs by volume fraction:
  - `:external`       — externally mixed (each species optically independent).
  - `:volume`         — volume-weighted average (Lorentz–Lorenz style).
  - `:maxwell_garnett`— matrix + inclusions (e.g. BC inclusions in a sulfate host).
  - `:bruggeman`      — symmetric effective medium.
Core/shell (coated BC) is a future extension; note it where relevant.
"""
const MixingRule = Symbol   # placeholder; gchp-io has concrete types — port those.

"""
    effective_ri(composition, ri_db, λ, rule, scheme) -> Complex

Per-bin (or per-mode) effective refractive index at wavelength `λ` (μm). `composition`
holds per-species volume/mass fractions; `ri_db` is the refractive-index database
(`data/refractive_indices_database.yaml`, interpolated in λ). Aerosol water (`AW`) enters
here as the hygroscopic-growth / RH effect on RI.

PORT: gchp-io already implements this for the AOD path; reuse verbatim.
"""
function effective_ri end   # TODO (port)

"""
    effective_ri_endpoints(composition, ri_db, λ_lo, λ_hi, rule, scheme) -> (n_lo, n_hi)

CORE-ENHANCEMENT HOOK (brief §6). vSmartMOM's band-endpoint Mie
(`model_from_parameters.jl:447-520`) currently uses a FIXED RI. To let RI vary across the
window, evaluate `effective_ri` at both band endpoints and thread (n_lo, n_hi) into the
endpoint Mie calls — mirroring the existing Sanghavi k(λ) endpoint interpolation. The
actual wiring lands in vSmartMOM core, not this package; this helper produces the inputs.
"""
function effective_ri_endpoints end   # TODO (new; pairs with a core change)
