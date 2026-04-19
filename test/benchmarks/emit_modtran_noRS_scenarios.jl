## =============================================================================
## emit_modtran_noRS_scenarios.jl
##
## Drives noRS (Rayleigh + aerosol + surface) forward simulations across all
## scenarios defined in ~/data/EMIT_MODTRANcomp/emit_20250404.json. For each
## scenario this script computes:
##   - Total (multiple scattering) reflected radiance at TOA
##   - Total transmitted radiance at BOA
##   - Single-scattering reflected radiance at TOA
##   - Single-scattering transmitted radiance at BOA
##
## The base atmospheric profile, geometry and arbitrary aerosol microphysics
## are defined in test/test_parameters/ParamsEMIT_MODTRANcomp.yaml; this
## script only modifies the LUT-grid dimensions (H2OSTR, AOT550, GNDALT) per
## scenario. TSZ is iterated but its physical meaning is ambiguous from the
## JSON alone — extend `apply_scenario!` below if it should influence the
## forward model.
##
## Results are saved as a JLD2 file. To keep turnaround reasonable, the
## default iterates over a coarse subset of the LUT grid; change STRIDE = 1
## below (or remove the stride slicing) to run the full grid.
## =============================================================================

using CUDA
CUDA.device!(1)

using Revise
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.InelasticScattering
using vSmartMOM.Absorption: read_hitran, make_hitran_model,
                            Voigt, HumlicekWeidemann32SDErrorFunction,
                            absorption_cross_section
using JSON
using JLD2
using NCDatasets
using Dates
using Printf

# ---- Paths (all absolute) ---------------------------------------------------

const REPO_ROOT    = "/home/sanghavi/code/github/vSmartMOM.jl"
const YAML_PATH    = "/home/sanghavi/code/github/vSmartMOM.jl/test/test_parameters/ParamsEMIT_MODTRANcomp.yaml"
const JSON_PATH    = "/home/sanghavi/data/EMIT_MODTRANcomp/emit_20250404.json"
const DAT_DIR      = "/home/sanghavi/data/EMIT_MODTRANcomp/vSmartMOM_out"
const OUTPUT_PATH  = "/home/sanghavi/data/EMIT_MODTRANcomp/emit_modtran_noRS_results.jld2"
const NC_PATH      = "/home/sanghavi/data/EMIT_MODTRANcomp/vSmartMOM_out/emit_noRS_vSmartMOM.nc"

# ---- Coarse-stride controls -------------------------------------------------
# The full LUT grid contains ~800 points which is impractical for a single
# run. Bump these down to 1 to sweep the full grid.
const STRIDE_H2O    = 2
const STRIDE_AOT    = 1
const STRIDE_GNDALT = 2
const STRIDE_TSZ    = 1

# ---- LUT-grid helpers -------------------------------------------------------

"Build a LUT axis from a {min, max, spacing} dict, inclusive of max."
function lut_axis(d::AbstractDict)
    lo   = Float32(d["min"])
    hi   = Float32(d["max"])
    step = Float32(d["spacing"])
    # `range(lo, hi, step=step)` drops `hi` if it isn't a clean multiple;
    # collect and re-filter to be inclusive-ish.
    vals = collect(lo:step:hi)
    if isempty(vals) || vals[end] < hi - 1e-9
        # `hi` was e.g. 5.001 to force inclusion of 5.0 — nothing to add.
    end
    return vals
end

# ---- Scenario application ---------------------------------------------------

"Approximate surface pressure [hPa] from ground altitude GNDALT [km]."
p_surface_from_gndalt(gndalt_km::Real) = Float32(1013.25f0 * exp(-Float32(gndalt_km) / 8.0f0))

"""
    interp_at_p(profile, p_centers, p_target)

Log-pressure linear interpolation of a layer-center `profile` (defined at
pressures `p_centers`) to the target pressure `p_target`. Clamps to the
boundary values when `p_target` lies outside the profile range.
"""
function interp_at_p(profile, p_centers, p_target)
    k = searchsortedlast(p_centers, p_target)
    k < 1              && return profile[1]
    k >= length(profile) && return profile[end]
    lp1  = log(p_centers[k])
    lp2  = log(p_centers[k+1])
    frac = (log(p_target) - lp1) / (lp2 - lp1)
    return profile[k] + frac * (profile[k+1] - profile[k])
end

"""
    apply_scenario!(params, base_q, base_p, base_T, base_vmr; h2o, aot, gndalt, tsz)

Mutate `params` in place for a single LUT-grid point. The `base_*` arguments
are copies of the original profiles (as loaded from YAML) so that every
scenario starts from the same baseline rather than accumulating edits.

When GNDALT raises the surface above some half-levels, the pressure grid is
truncated: all half-levels ≥ p_surf are removed and replaced by p_surf, and
the layer-center profiles (T, q, VMRs) are interpolated to the new surface.
"""
function apply_scenario!(params, base_q, base_p, base_T, base_vmr;
                         h2o::Real, aot::Real,
                         gndalt::Real, tsz::Real)

    p_surf = p_surface_from_gndalt(gndalt)

    # --- AOT550: set the reference aerosol optical depth at λ_ref = 0.55 µm.
    params.scattering_params.rt_aerosols[1].τ_ref = aot

    # --- GNDALT + profiles: truncate the grid if the new surface is above
    #     any existing half-levels.
    # k = last index in base_p that is strictly below p_surf.
    k = searchsortedlast(base_p, p_surf - eps(Float32(p_surf)))

    if k >= length(base_p)
        # No truncation needed — new surface is at or below the original.
        params.p = copy(base_p)
        params.p[end] = p_surf
        params.T = copy(base_T)
        params.q = base_q .* h2o
        for (gas, val) in base_vmr
            params.absorption_params.vmr[gas] = val isa AbstractVector ? copy(val) : val
        end
    else
        # Truncate: keep half-levels 1:k, append p_surf.
        new_p = vcat(base_p[1:k], p_surf)
        n_layers = k   # one fewer than number of half-levels

        # Layer-center pressures on the original (full) grid.
        p_centers = [(base_p[i] + base_p[i+1]) / 2 for i in 1:length(base_p)-1]

        # Helper: truncate a layer-center profile to n_layers, interpolating
        # the last value at p_surf.
        function truncate_profile(prof)
            out = similar(prof, n_layers)
            out[1:k-1] .= prof[1:k-1]
            out[k] = interp_at_p(prof, p_centers, p_surf)
            return out
        end

        params.p = new_p
        params.T = truncate_profile(base_T)
        params.q = truncate_profile(base_q) .* h2o

        for (gas, val) in base_vmr
            params.absorption_params.vmr[gas] =
                val isa AbstractVector ? truncate_profile(val) : val
        end
    end

    # --- TSZ: MODTRAN observer zenith angle (180° = nadir).
    # Convert to vSmartMOM VZA (0° = nadir): VZA = 180 − TSZ.
    params.vza = [180 - tsz]
    nothing
end

# ---- Fast profile recomputation ---------------------------------------------

"""
    recompute_vcd(T, p_half, q)

Recompute `vcd_dry` and `vmr_h2o` from a new specific-humidity profile `q`,
using the same arithmetic as `compute_atmos_profile_fields` but without
touching HITRAN or building a full model.  Returns `(vcd_dry, vmr_h2o)`.
"""
function recompute_vcd(T, p_half, q; g₀=9.8032465)
    FT = eltype(T)
    Nₐ        = FT(6.02214179e+23)
    dry_mass   = FT(28.9644e-3)
    wet_mass   = FT(18.01534e-3)
    n_layers   = length(T)
    vmr_h2o    = zeros(FT, n_layers)
    vcd_dry    = zeros(FT, n_layers)
    for i in 1:n_layers
        Δp         = p_half[i+1] - p_half[i]
        vmr_h2o[i] = (dry_mass / wet_mass) * q[i] / (1 - q[i])
        vmr_dry    = 1 - vmr_h2o[i]
        M          = vmr_dry * dry_mass + vmr_h2o[i] * wet_mass
        vcd        = Nₐ * Δp / (M * g₀ * 100^2) * 100
        vcd_dry[i] = vmr_dry * vcd
    end
    return vcd_dry, vmr_h2o
end

# ---- GNDALT truncation helpers ----------------------------------------------

"""
    build_absorption_models(params, profile)

Build lightweight absorption-model objects for every molecule that contributes
to τ_abs.  These are cheap to construct (just struct wrappers around parsed
HITRAN data or LUT handles) and let us call `absorption_cross_section` for a
single (p, T) later without repeating the full line-by-line sweep.

Returns `(mol_models, h2o_q_model)` where
  • `mol_models::Vector{NamedTuple}` — one entry per molecule with fields
    `(:name, :model, :vmr, :is_fixed)`.  `vmr` is scalar or Vector on the
    *full* profile grid.
  • `h2o_q_model` — the HitranModel for the q-based H₂O (wing_cutoff = 150).
"""
function build_absorption_models(params, profile)
    mol_models = []

    # Fixed molecules (O2 via LUT, CO2 via HITRAN, etc.)
    for mol in params.absorption_params.fixed_molecules[1]
        if haskey(params.absorption_params.luts, mol)
            m = params.absorption_params.luts[mol]
        else
            hd = read_hitran(artifact(mol), iso=1)
            m  = make_hitran_model(hd,
                     params.absorption_params.broadening_function,
                     wing_cutoff = params.absorption_params.wing_cutoff,
                     CEF = params.absorption_params.CEF,
                     architecture = params.architecture, vmr = 0)
        end
        push!(mol_models, (name=mol, model=m,
              vmr=profile.vmr[mol], is_fixed=true))
    end

    # Variable molecules (H2O, CH4, O3, N2O, CO — all HITRAN or LUT)
    for mol in params.absorption_params.variable_molecules[1]
        if haskey(params.absorption_params.luts, mol)
            m = params.absorption_params.luts[mol]
        else
            hd = read_hitran(artifact(mol), iso=1)
            m  = make_hitran_model(hd,
                     params.absorption_params.broadening_function,
                     wing_cutoff = params.absorption_params.wing_cutoff,
                     CEF = params.absorption_params.CEF,
                     architecture = params.architecture, vmr = 0)
        end
        push!(mol_models, (name=mol, model=m,
              vmr=profile.vmr[mol], is_fixed=false))
    end

    # q-based H₂O (special wing_cutoff = 150)
    hd = read_hitran(artifact("H2O"), iso=1)
    h2o_q_model = make_hitran_model(hd, Voigt(),
                      wing_cutoff = 150,
                      CEF = HumlicekWeidemann32SDErrorFunction(),
                      architecture = params.architecture, vmr = 0)

    return mol_models, h2o_q_model
end

"""
    compute_bottom_layer_τ!(τ_layer, τ_qH2O_layer,
                            mol_models, h2o_q_model, grid,
                            p, T, vcd_dry, vmr_h2o, vmr_dict)

Compute absorption optical depth for a single layer at pressure `p` and
temperature `T`.  Writes the total into `τ_layer` and the q-based H₂O
component into `τ_qH2O_layer` (both pre-allocated vectors of length n_spec).
"""
function compute_bottom_layer_τ!(τ_layer, τ_qH2O_layer,
                                 mol_models, h2o_q_model, grid,
                                 p, T, vcd_dry_k, vmr_h2o_k, vmr_dict)
    fill!(τ_layer, 0)
    fill!(τ_qH2O_layer, 0)

    # q-based H₂O
    σ = Array(absorption_cross_section(h2o_q_model, grid, p, T))
    τ_qH2O_layer .= σ .* vcd_dry_k .* vmr_h2o_k
    τ_layer     .+= τ_qH2O_layer

    # All other molecules
    for entry in mol_models
        vmr_k = vmr_dict[entry.name]
        σ = Array(absorption_cross_section(entry.model, grid, p, T))
        τ_layer .+= σ .* vcd_dry_k .* vmr_k
    end
end

"""
    truncate_and_assemble(profile_full, τ_abs_full, τ_abs_qH2O_full,
                          τ_rayl_full, τ_aer_full,
                          mol_models, h2o_q_model, grid,
                          base_q_full, base_p_full, base_T_full, base_vmr_full,
                          p_surf, aer_params, depol_rayl, curr_band_λ)

Given the full-atmosphere cached arrays and a new surface pressure `p_surf`,
return truncated `(τ_abs, τ_abs_qH2O, τ_rayl, τ_aer, profile_trunc, q_trunc)`
by slicing the first k-1 layers unchanged and recomputing only layer k.
"""
function truncate_and_assemble(profile_full, τ_abs_full, τ_abs_qH2O_full,
                               τ_rayl_full, τ_aer_full,
                               mol_models, h2o_q_model, grid,
                               base_q_full, base_T_full,
                               p_surf, aer_params, depol_rayl, curr_band_λ)

    FT  = eltype(profile_full.T)
    p_half_full = profile_full.p_half
    p_surf_FT   = FT(p_surf)              # promote once; keep all arithmetic in FT

    # ---- Truncation index k ----
    k = searchsortedlast(p_half_full, p_surf_FT - eps(p_surf_FT))

    if k >= length(p_half_full)
        # No truncation — return copies of the full arrays with p_half[end] = p_surf
        p_half_new = copy(p_half_full)
        p_half_new[end] = p_surf_FT
        p_full_new = (p_half_new[1:end-1] .+ p_half_new[2:end]) ./ FT(2)
        # Recompute vcd for the adjusted bottom boundary
        q_new = copy(base_q_full)
        vcd_dry_new, vmr_h2o_new = recompute_vcd(base_T_full, p_half_new, q_new)
        prof = CoreRT.AtmosphericProfile(copy(base_T_full), p_full_new, q_new,
                   p_half_new, vmr_h2o_new, vcd_dry_new,
                   zeros(FT, length(base_T_full)), # vcd_h2o placeholder
                   Dict{String, Union{Real, Vector}}(
                       g => v isa AbstractVector ? copy(v) : v
                       for (g,v) in profile_full.vmr))
        # Scale τ_rayl and τ_aer for the tiny p_half[end] change
        s = reshape(vcd_dry_new ./ profile_full.vcd_dry, 1, :)
        return (τ_abs      = τ_abs_full .* s,
                τ_abs_qH2O = τ_abs_qH2O_full .* s,
                τ_rayl     = τ_rayl_full .* s,
                τ_aer      = copy(τ_aer_full),
                profile    = prof,
                q_base     = q_new)
    end

    # ---- Build truncated half-level grid ----
    n_layers = k
    p_half_new = vcat(p_half_full[1:k], p_surf_FT)
    p_full_new = (p_half_new[1:end-1] .+ p_half_new[2:end]) ./ FT(2)

    # ---- Truncated T, q, vmr profiles ----
    # Layer-center pressures on the original (full) grid for interpolation
    p_centers_full = (p_half_full[1:end-1] .+ p_half_full[2:end]) ./ FT(2)

    T_trunc     = zeros(FT, n_layers)
    q_trunc     = zeros(FT, n_layers)
    T_trunc[1:k-1]  .= base_T_full[1:k-1]
    q_trunc[1:k-1]  .= base_q_full[1:k-1]
    T_trunc[k]  = interp_at_p(base_T_full, p_centers_full, p_full_new[k])
    q_trunc[k]  = interp_at_p(base_q_full, p_centers_full, p_full_new[k])

    vmr_trunc = Dict{String, Union{Real, Vector}}()
    for (g, v) in profile_full.vmr
        if v isa AbstractVector
            vt = zeros(FT, n_layers)
            vt[1:k-1] .= v[1:k-1]
            vt[k] = interp_at_p(v, p_centers_full, p_full_new[k])
            vmr_trunc[g] = vt
        else
            vmr_trunc[g] = v
        end
    end

    vcd_dry_new, vmr_h2o_new = recompute_vcd(T_trunc, p_half_new, q_trunc)

    prof = CoreRT.AtmosphericProfile(T_trunc, p_full_new, q_trunc,
               p_half_new, vmr_h2o_new, vcd_dry_new,
               zeros(FT, n_layers),  # vcd_h2o placeholder
               vmr_trunc)

    n_spec = size(τ_abs_full, 1)

    # ---- τ_abs: reuse layers 1:k-1, recompute layer k ----
    τ_abs_trunc      = zeros(FT, n_spec, n_layers)
    τ_abs_qH2O_trunc = zeros(FT, n_spec, n_layers)
    τ_abs_trunc[:,      1:k-1] .= τ_abs_full[:,      1:k-1]
    τ_abs_qH2O_trunc[:, 1:k-1] .= τ_abs_qH2O_full[:, 1:k-1]

    # Bottom layer k: evaluate cross-sections at one (p, T) point
    vmr_bottom = Dict{String, Real}()
    for entry in mol_models
        vmr_bottom[entry.name] = entry.vmr isa AbstractVector ?
            interp_at_p(entry.vmr, p_centers_full, p_full_new[k]) :
            entry.vmr
    end

    τ_k      = zeros(FT, n_spec)
    τ_qH2O_k = zeros(FT, n_spec)
    compute_bottom_layer_τ!(τ_k, τ_qH2O_k,
        mol_models, h2o_q_model, grid,
        p_full_new[k], T_trunc[k], vcd_dry_new[k], vmr_h2o_new[k],
        vmr_bottom)
    τ_abs_trunc[:, k]      .= τ_k
    τ_abs_qH2O_trunc[:, k] .= τ_qH2O_k

    # ---- τ_rayl: recompute (fast arithmetic) ----
    τ_rayl_trunc = CoreRT.getRayleighLayerOptProp(
                       p_surf_FT, curr_band_λ, depol_rayl, vcd_dry_new)

    # ---- τ_aer: recompute vertical profile (fast LogNormal) ----
    n_aer  = length(aer_params)
    τ_aer_trunc = zeros(FT, n_aer, n_spec, n_layers)
    for (ia, ap) in enumerate(aer_params)
        vert = CoreRT.getAerosolLayerOptProp(1, ap.z₀, ap.σ₀, p_half_new, T_trunc)
        τ_aer_trunc[ia, :, :] .= (ap.τ_ref / ap.k_ref) .* ap.k_spec .* vert'
    end

    return (τ_abs = τ_abs_trunc, τ_abs_qH2O = τ_abs_qH2O_trunc,
            τ_rayl = τ_rayl_trunc, τ_aer = τ_aer_trunc,
            profile = prof, q_base = q_trunc)
end

# ---- Result extraction ------------------------------------------------------

"""
Return (R_TOA_I, T_BOA_I) as 1-D spectra (n_spec,) for the first VZA and the
I Stokes component. These are the reflected radiance at TOA and the
transmitted radiance at BOA as produced by `rt_run` / `rt_run_ss`.
"""
function extract_I_radiances(R_SFI::AbstractArray, T_SFI::AbstractArray)
    R_I = Array(R_SFI)[1, 1, :]
    T_I = Array(T_SFI)[1, 1, :]
    return R_I, T_I
end

# ---- Per-scenario .dat output -----------------------------------------------

"""
    write_scenario_dat(dir, λ_nm, ν, R_tot, T_tot, R_ss, T_ss,
                       hem_R, hem_T, hem_R_ss, hem_T_ss;
                       h2o, aot, gndalt, tsz)

Write a single-scenario text file with columns:
    wl_nm  wn_cm-1  R_tot  T_tot  R_ss  T_ss  hem_R  hem_T  hem_R_ss  hem_T_ss
The filename encodes the four LUT-grid coordinates.
"""
function write_scenario_dat(dir, λ_nm, ν, R_tot_I, T_tot_I, R_ss_I, T_ss_I,
                            hem_R, hem_T, hem_R_ss, hem_T_ss;
                            h2o, aot, gndalt, tsz)
    fname = @sprintf("vSmartMOM_H2O%.4f_AOT%.4f_GNDALT%.3f_TSZ%.1f.dat",
                     h2o, aot, gndalt, tsz)
    path = joinpath(dir, fname)
    open(path, "w") do io
        @printf(io, "# vSmartMOM noRS   H2OSTR=%.4f  AOT550=%.4f  GNDALT=%.3f km  TSZ=%.1f\n",
                h2o, aot, gndalt, tsz)
        @printf(io, "# %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n",
                "wl_nm", "wn_cm-1", "R_tot", "T_tot", "R_ss", "T_ss",
                "hem_R", "hem_T", "hem_R_ss", "hem_T_ss")
        for i in eachindex(λ_nm)
            @printf(io, "  %16.6f %16.6f %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
                    λ_nm[i], ν[i], R_tot_I[i], T_tot_I[i], R_ss_I[i], T_ss_I[i],
                    hem_R[i], hem_T[i], hem_R_ss[i], hem_T_ss[i])
        end
    end
    return path
end

# ---- NetCDF output ----------------------------------------------------------

"""
    write_results_nc(path, h2o_axis, aot_axis, gndalt_axis, tsz_axis,
                     λ_nm, ν_axis, R_tot, T_tot, R_ss, T_ss, metadata)

Write all simulation results to a NetCDF-4 file with dimensions
(H2OSTR, AOT550, GNDALT, TSZ, wl) mirroring the MODTRAN .nc layout.
"""
function write_results_nc(path, h2o_axis, aot_axis, gndalt_axis, tsz_axis,
                          λ_nm, ν_axis, R_tot, T_tot, R_ss, T_ss,
                          hemR_tot, hemT_tot, hemR_ss, hemT_ss, metadata)
    ds = NCDataset(path, "c")

    # Global attributes
    ds.attrib["title"]       = "vSmartMOM noRS radiance LUT"
    ds.attrib["source"]      = "vSmartMOM.jl"
    ds.attrib["json_source"] = metadata["json_source"]
    ds.attrib["yaml_source"] = metadata["yaml_source"]
    ds.attrib["timestamp"]   = metadata["timestamp"]
    ds.attrib["n_total"]     = metadata["n_total"]

    # Dimensions
    nH  = length(h2o_axis)
    nA  = length(aot_axis)
    nG  = length(gndalt_axis)
    nT  = length(tsz_axis)
    nW  = length(λ_nm)

    defDim(ds, "H2OSTR",              nH)
    defDim(ds, "AOT550",              nA)
    defDim(ds, "surface_elevation_km", nG)
    defDim(ds, "TSZ",                 nT)
    defDim(ds, "wl",                  nW)

    # Coordinate variables
    v = defVar(ds, "H2OSTR",              Float64, ("H2OSTR",))
    v[:] = Float64.(h2o_axis)

    v = defVar(ds, "AOT550",              Float64, ("AOT550",))
    v[:] = Float64.(aot_axis)

    v = defVar(ds, "surface_elevation_km", Float64, ("surface_elevation_km",))
    v[:] = Float64.(gndalt_axis)

    v = defVar(ds, "TSZ",                 Float64, ("TSZ",))
    v[:] = Float64.(tsz_axis)

    v = defVar(ds, "wl",                  Float64, ("wl",))
    v.attrib["units"] = "nm"
    v[:] = Float64.(λ_nm)

    v = defVar(ds, "wn",                  Float64, ("wl",))
    v.attrib["units"] = "cm-1"
    v[:] = Float64.(ν_axis)

    # Data variables — dimensions match (H2OSTR, AOT550, GNDALT, TSZ, wl)
    data_dims = ("H2OSTR", "AOT550", "surface_elevation_km", "TSZ", "wl")

    for (name, arr, desc) in [
        ("R_tot",    R_tot,    "TOA reflected radiance (multiple scattering)"),
        ("T_tot",    T_tot,    "BOA transmitted radiance (multiple scattering)"),
        ("R_ss",     R_ss,     "TOA reflected radiance (single scattering)"),
        ("T_ss",     T_ss,     "BOA transmitted radiance (single scattering)"),
        ("hemR_tot", hemR_tot, "Hemispheric reflectance (multiple scattering)"),
        ("hemT_tot", hemT_tot, "Hemispheric transmittance (multiple scattering)"),
        ("hemR_ss",  hemR_ss,  "Hemispheric reflectance (single scattering)"),
        ("hemT_ss",  hemT_ss,  "Hemispheric transmittance (single scattering)"),
    ]
        v = defVar(ds, name, Float64, data_dims)
        v.attrib["description"] = desc
        v[:] = arr
    end

    close(ds)
    return path
end

# ---- Main driver ------------------------------------------------------------
#
# Optimization strategy
# ---------------------
# `model_from_parameters` is the dominant cost (~30 min of HITRAN line-by-line
# per call).  We exploit two facts:
#
#   1. GNDALT only truncates the *bottom* of the atmosphere.  Layers 1 to k−1
#      keep identical (T, p), so their absorption cross-sections are unchanged.
#      Only the new bottom layer k needs a fresh σ(ν, T_k, p_k) evaluation.
#      → ONE full model build at GNDALT = 0; seconds per additional GNDALT.
#
#   2. H2OSTR, AOT550, TSZ do not change the (T, p) grid:
#      – H2OSTR scales q → changes vmr_h2o / vcd_dry (fast arithmetic).
#      – AOT550 is a multiplicative factor on τ_aer.
#      – TSZ is currently unused.
#      → Pure scaling within the inner loop.
#
# This reduces the expensive computation from n_total × 30 min
# to ONE × 30 min + seconds × (n_GNDALT − 1).
#

function run_all_scenarios()

    mkpath(DAT_DIR)

    @info "Loading JSON scenario descriptor" JSON_PATH
    cfg    = JSON.parsefile(JSON_PATH)
    lut    = cfg["lut_grid"]
    misc   = cfg["MISC"]

    h2o_axis    = lut_axis(lut["H2OSTR"])[1:STRIDE_H2O:end]
    aot_axis    = lut_axis(lut["AOT550"])[1:STRIDE_AOT:end]
    gndalt_axis = lut_axis(lut["GNDALT"])[1:STRIDE_GNDALT:end]
    tsz_axis    = lut_axis(lut["TSZ"])[1:STRIDE_TSZ:end]

    @info "LUT grid dimensions" H2OSTR=length(h2o_axis) AOT550=length(aot_axis) GNDALT=length(gndalt_axis) TSZ=length(tsz_axis)
    n_total = length(h2o_axis) * length(aot_axis) * length(gndalt_axis) * length(tsz_axis)
    n_inner = length(h2o_axis) * length(aot_axis) * length(tsz_axis)
    @info "Total scenarios to run" n_total n_model_builds=1

    # Load the YAML once to get the baseline parameters.
    params  = parameters_from_yaml(YAML_PATH)
    FT = params.float_type

    # Snapshots of mutable arrays so every scenario starts from the same state.
    base_q = copy(params.q)
    base_p = copy(params.p)
    base_T = copy(params.T)
    base_vmr = Dict(g => v isa AbstractVector ? copy(v) : v
                    for (g, v) in params.absorption_params.vmr)

    @info "Geometry from YAML (should match MISC in JSON)" sza=params.sza vaz=params.vaz

    # Pre-compute wavelength axis for the first (and only default) band.
    ν_axis   = collect(params.spec_bands[1])
    λ_nm     = 1e7 ./ ν_axis
    n_spec   = length(ν_axis)

    # Output containers.  Axes: (h2o, aot, gndalt, tsz, spec).
    out_dims = (length(h2o_axis), length(aot_axis),
                length(gndalt_axis), length(tsz_axis), n_spec)
    R_tot    = zeros(Float64, out_dims)
    T_tot    = zeros(Float64, out_dims)
    R_ss     = zeros(Float64, out_dims)
    T_ss     = zeros(Float64, out_dims)
    hemR_tot = zeros(Float64, out_dims)
    hemT_tot = zeros(Float64, out_dims)
    hemR_ss  = zeros(Float64, out_dims)
    hemT_ss  = zeros(Float64, out_dims)

    scenario_k = 0
    t_start    = time()
    h2o_base   = FT(1)
    aot_base   = aot_axis[1]

    # ==================================================================
    #  ONE-TIME: build full model at GNDALT = 0
    # ==================================================================
    @info "Building full model at GNDALT = 0 (one-time HITRAN computation)"
    apply_scenario!(params, base_q, base_p, base_T, base_vmr;
                    h2o = h2o_base, aot = aot_base,
                    gndalt = FT(0), tsz = tsz_axis[1])
    model = model_from_parameters(params)

    # Cache the full-atmosphere optical depths (at h2o=1, aot=aot_base).
    profile_full  = model.profile
    τ_abs_full    = copy(model.τ_abs[1])
    τ_rayl_full   = copy(model.τ_rayl[1])
    τ_aer_full    = copy(model.τ_aer[1])

    # Pre-build absorption models for every molecule (for bottom-layer recomputation).
    @info "Building per-molecule absorption models"
    mol_models, h2o_q_model = build_absorption_models(params, profile_full)

    # Compute q-based H₂O absorption for the full profile (for scaling).
    @info "Computing q-based H₂O absorption (full atmosphere)"
    τ_abs_qH2O_full = zeros(FT, size(τ_abs_full))
    CoreRT.compute_absorption_profile!(τ_abs_qH2O_full, h2o_q_model,
        params.spec_bands[1], profile_full.vmr_h2o, profile_full)

    # Rayleigh depolarisation (needed for τ_rayl recomputation per GNDALT).
    curr_band_λ = FT.(1e4 ./ params.spec_bands[1])
    νₘ  = FT(0.5) * (params.spec_bands[1][1] + params.spec_bands[1][end])
    λₘ  = FT(1e7) / νₘ
    effT_full = (profile_full.vcd_dry' * profile_full.T) / sum(profile_full.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(νₘ, effT_full)
    _, _   = InelasticScattering.compute_γ_air_Cabannes!(λₘ, n2, o2)
    γ_air_Rayleigh, _ = InelasticScattering.compute_γ_air_Rayleigh!(λₘ, n2, o2)
    depol_rayl = FT(2γ_air_Rayleigh / (1 + γ_air_Rayleigh))

    # Aerosol parameters (for τ_aer recomputation per GNDALT).
    # Capture k_ref by replicating the tiny Mie reference-extinction calculation.
    aer_info = map(params.scattering_params.rt_aerosols) do c_aero
        mie_model = vSmartMOM.Scattering.make_mie_model(
            params.scattering_params.decomp_type,
            vSmartMOM.Scattering.Aerosol(c_aero.aerosol.size_distribution,
                                         c_aero.aerosol.nᵣ, c_aero.aerosol.nᵢ),
            params.scattering_params.λ_ref,
            params.polarization_type,
            vSmartMOM.Scattering.δBGE{FT}(params.l_trunc, params.Δ_angle),
            params.scattering_params.r_max,
            params.scattering_params.nquad_radius)
        k_ref = FT(vSmartMOM.Scattering.compute_ref_aerosol_extinction(mie_model, FT))
        (z₀ = FT(c_aero.z₀), σ₀ = FT(c_aero.σ₀), τ_ref = FT(c_aero.τ_ref),
         k_ref = k_ref,
         k_spec = FT.(model.aerosol_optics[1][1].k))   # spectral extinction (n_spec,)
    end

    elapsed_build = time() - t_start
    @info "Full model + separation done" elapsed_s = round(elapsed_build; digits=1)

    # ==================================================================
    #  Outer loop: GNDALT  (cheap — truncate + recompute bottom layer)
    # ==================================================================
    for (ig, gndalt) in enumerate(gndalt_axis)

        t_gndalt = time()
        @info "=== GNDALT $ig / $(length(gndalt_axis)) ===" gndalt

        p_surf = p_surface_from_gndalt(gndalt)

        trunc = truncate_and_assemble(
                    profile_full, τ_abs_full, τ_abs_qH2O_full,
                    τ_rayl_full, τ_aer_full,
                    mol_models, h2o_q_model, params.spec_bands[1],
                    base_q, base_T,
                    p_surf, aer_info, depol_rayl, curr_band_λ)

        # Patch the model for this GNDALT.
        τ_abs_all    = trunc.τ_abs
        τ_abs_qH2O   = trunc.τ_abs_qH2O
        τ_abs_rest   = τ_abs_all .- τ_abs_qH2O
        τ_rayl_base  = trunc.τ_rayl
        τ_aer_base   = trunc.τ_aer
        q_base_g     = trunc.q_base
        p_half_g     = trunc.profile.p_half
        T_g          = trunc.profile.T
        vcd_dry_base = trunc.profile.vcd_dry
        vmr_h2o_base = trunc.profile.vmr_h2o
        n_layers     = length(T_g)

        # Replace model arrays (sizes may differ between GNDALTs).
        model.τ_abs[1]  = copy(τ_abs_all)
        model.τ_rayl[1] = copy(τ_rayl_base)
        model.τ_aer[1]  = copy(τ_aer_base)
        model.profile    = trunc.profile

        elapsed_trunc = time() - t_gndalt
        @info "GNDALT truncation done" gndalt n_layers elapsed_s=round(elapsed_trunc; digits=1)

        # --- Build RS_type for this GNDALT ----
        n_pol        = model.params.polarization_type.n
        n_spec_model = size(τ_rayl_base, 1)
        F₀_vec       = zeros(FT, n_pol); F₀_vec[1] = 1
        SIF₀_vec     = zeros(FT, n_pol)
        RS_type = InelasticScattering.noRS(
                    fscattRayl  = [FT(1)],
                    ϖ_Cabannes  = [FT(1)],
                    bandSpecLim = [],
                    iBand       = [1],
                    F₀          = zeros(FT, n_pol, n_spec_model),
                    SIF₀        = zeros(FT, n_pol, n_spec_model))
        RS_type.F₀  .= F₀_vec
        RS_type.SIF₀ .= SIF₀_vec

        # ==============================================================
        #  Middle loop: TSZ  (update VZA + quadrature; no absorption cost)
        # ==============================================================
        for (it, tsz) in enumerate(tsz_axis)

            vza_new = FT(180 - tsz)
            model.obs_geom.vza = [vza_new]
            model.quad_points  = CoreRT.rt_set_streams(
                                     params.quadrature_type, params.l_trunc,
                                     model.obs_geom, params.polarization_type,
                                     vSmartMOM.Architectures.array_type(params.architecture))
            @info "  TSZ=$tsz (VZA=$(vza_new)°) — quadrature updated"

            # ==========================================================
            #  Inner loop: H2OSTR × AOT550  (cheap — scaling only)
            # ==========================================================
            for (ih, h2o) in enumerate(h2o_axis),
                (ia, aot) in enumerate(aot_axis)

                scenario_k += 1
                @info "Scenario $scenario_k / $n_total" gndalt tsz h2o aot

                # Recompute vcd_dry and vmr_h2o for this H2OSTR (fast).
                q_new = q_base_g .* (h2o / h2o_base)
                vcd_dry_new, vmr_h2o_new = recompute_vcd(T_g, p_half_g, q_new)

                # Per-layer scaling factors.
                s_vcd = reshape(vcd_dry_new ./ vcd_dry_base, 1, :)
                s_h2o = reshape((vmr_h2o_new ./ vmr_h2o_base) .* vec(s_vcd), 1, :)

                # Scale optical depths.
                model.τ_abs[1]  .= τ_abs_rest .* s_vcd .+ τ_abs_qH2O .* s_h2o
                model.τ_rayl[1] .= τ_rayl_base .* s_vcd
                model.τ_aer[1]  .= τ_aer_base .* (aot / aot_base)

                # Full multiple-scattering run.
                R_SFI_tot, T_SFI_tot, _, _, hem_R_tot, hem_T_tot = rt_run(RS_type, model, 1)

                # Single-scattering run.
                R_SFI_ss,  T_SFI_ss,  _, _, hem_R_ss, hem_T_ss = rt_run_ss(RS_type, model, 1)

                R_tot_I, T_tot_I = extract_I_radiances(R_SFI_tot, T_SFI_tot)
                R_ss_I,  T_ss_I  = extract_I_radiances(R_SFI_ss,  T_SFI_ss)

                R_tot[ih, ia, ig, it, :] .= R_tot_I
                T_tot[ih, ia, ig, it, :] .= T_tot_I
                R_ss[ih, ia, ig, it, :] .= R_ss_I
                T_ss[ih, ia, ig, it, :] .= T_ss_I
                hemR_tot[ih, ia, ig, it, :] .= Float64.(hem_R_tot)
                hemT_tot[ih, ia, ig, it, :] .= Float64.(hem_T_tot)
                hemR_ss[ih, ia, ig, it, :] .= Float64.(hem_R_ss)
                hemT_ss[ih, ia, ig, it, :] .= Float64.(hem_T_ss)

                write_scenario_dat(DAT_DIR, λ_nm, ν_axis,
                                   R_tot_I, T_tot_I, R_ss_I, T_ss_I,
                                   hem_R_tot, hem_T_tot, hem_R_ss, hem_T_ss;
                                   h2o=h2o, aot=aot, gndalt=gndalt, tsz=tsz)

                elapsed = time() - t_start
                @info "Scenario complete" scenario_k elapsed_s=round(elapsed; digits=1)
            end
        end

        elapsed_gndalt = time() - t_gndalt
        @info "=== GNDALT $ig done ===" gndalt n_scenarios=n_inner elapsed_s=round(elapsed_gndalt; digits=1)
    end

    metadata = Dict(
        "json_source"    => JSON_PATH,
        "yaml_source"    => YAML_PATH,
        "timestamp"      => string(Dates.now()),
        "misc_from_json" => misc,
        "n_total"        => n_total,
        "stride"         => (H2O=STRIDE_H2O, AOT=STRIDE_AOT,
                             GNDALT=STRIDE_GNDALT, TSZ=STRIDE_TSZ),
    )

    @info "Saving JLD2" OUTPUT_PATH
    @save OUTPUT_PATH h2o_axis aot_axis gndalt_axis tsz_axis ν_axis λ_nm R_tot T_tot R_ss T_ss hemR_tot hemT_tot hemR_ss hemT_ss metadata

    @info "Writing NetCDF" NC_PATH
    write_results_nc(NC_PATH, h2o_axis, aot_axis, gndalt_axis, tsz_axis,
                     λ_nm, ν_axis, R_tot, T_tot, R_ss, T_ss,
                     hemR_tot, hemT_tot, hemR_ss, hemT_ss, metadata)

    @info "All output written" dat_dir=DAT_DIR nc=NC_PATH jld2=OUTPUT_PATH

    return (; h2o_axis, aot_axis, gndalt_axis, tsz_axis,
              ν_axis, λ_nm, R_tot, T_tot, R_ss, T_ss,
              hemR_tot, hemT_tot, hemR_ss, hemT_ss, metadata)
end

# Execute when run as a script.
results = run_all_scenarios()