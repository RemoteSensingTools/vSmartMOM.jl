# Stage 3 (NEW — the core missing piece): sectional/parametric state → per-layer
# AerosolOptics that the RT solver consumes. There is NO reference implementation; the AOD
# path explicitly does not produce phase matrices.
#
# What the RT wants per band, per layer (see vSmartMOM
# src/CoreRT/tools/model_from_parameters.jl): an `AerosolOptics` with extinction `k`,
# single-scattering albedo `ϖ`, and GreekCoefs (α, β, γ, δ, ϵ, ζ); plus the per-layer τ
# that feeds `τ_aer[iBand][iAer, nSpec, iZ]`.

"""
    layer_aerosol_optics(scene, ilev, λ; rep, mixing, ri_db, decomp=vSmartMOM.NAI2())
        -> vSmartMOM.Scattering.AerosolOptics

Build the full scattering optics for one layer at wavelength `λ`. Two strategies:

  STRATEGY 1 — per-bin Mie → summed Greek (pairs with `Sectional`/`SplineSectional`):
    for each bin: effective RI → `make_mie_model` / `compute_aerosol_optical_properties`,
    weight phase matrices by number density, sum → one layer `AerosolOptics`; normalize ϖ.

  STRATEGY 2 — lognormal-fit → core Mie (pairs with `LogNormalFit`):
    emit `vSmartMOM.Aerosol(LogNormal(μ,σ), n_eff)` and call the existing core Mie path
    directly. Less new code; inherits the fit approximation.

Reuse `vSmartMOM.Scattering` for the Mie + Greek machinery — do not reimplement Mie.
"""
function layer_aerosol_optics end   # TODO (new — Stage 3)

"""
    column_aerosol_optics(scene, band; rep, mixing, ri_db) -> Vector{AerosolOptics}

All layers for a band. Decision knob (brief §5): does the size-distribution SHAPE vary with
height (call `layer_aerosol_optics` per layer) or only the LOADING (compute one
AOD-weighted mean distribution, scale by per-layer concentration)? Expose both; default to
the faithful per-layer path and measure the radiance difference.
"""
function column_aerosol_optics end   # TODO (new — Stage 3)

"""
    scene_to_rt_parameters(scene, base_params; rep, mixing, ri_db) -> vSmartMOM parameters

Assemble the per-column "model setup" that drives an end-to-end RT run (Stage 5): inject the
per-layer aerosol optics + τ profile into a vSmartMOM parameter/model object so
`vSmartMOM.rt_run` can produce radiances. PARTIAL prior art: gchp-io
`parameters_from_scene` / `scene_to_dict` (atmosphere side; aerosol-optics side is new).
"""
function scene_to_rt_parameters end   # TODO (new — Stage 5)
