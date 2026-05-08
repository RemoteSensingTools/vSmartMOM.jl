# =============================================================================
# emit_modtran_noRS_scenarios.jl — noRS forward LUT driver across MODTRAN
# scenarios (Phase 6 port from sanghavi).
#
# Drives noRS (Rayleigh + aerosol + Lambertian surface α=0) forward simulations
# across the scenarios defined in a JSON LUT-grid descriptor (default
# `$HOME/data/EMIT_MODTRANcomp/emit_20250404.json`). For each scenario the
# script computes multiple-scattering TOA R + BOA T, the total vertical
# atmospheric optical thickness, and hemispheric flux integrals; writes
# per-scenario `.dat` files and an aggregate NetCDF.
#
# 2026-04-24 rewrite:
#   * Single-scattering outputs removed. The MODTRAN-equivalent decomposition
#     now uses the two-albedo flux closure (bhr_dw at α=0 and α=0.2) and no
#     longer needs a single-scatter proxy for the spherical albedo. See
#     `plans/MODTRAN_equivalent_equations.md` for the derivation.
#   * Total column τ (absorption + Rayleigh + aerosol extinction, vertical)
#     is written as a spectrum per scenario.
#   * Default YAML switched to ParamsEMIT_MODTRANcomp_newLUT.yaml (rebuilt
#     HITRAN LUTs under ~/data/HITRAN_LUTs/).
#   * AOT axis is hard-clamped to AOT=0.0001 scenarios (stride 4 on the
#     four-point descriptor axis) — MODTRAN comparisons at this leg use
#     AOT=0.0001 only.
#
# ENV overrides:
#   EMIT_YAML          default test_parameters/ParamsEMIT_MODTRANcomp_newLUT.yaml
#   EMIT_JSON          default ~/data/EMIT_MODTRANcomp/emit_20250404.json
#   EMIT_NORS_DAT_DIR  default ~/data/EMIT_MODTRANcomp/vSmartMOM_out
#   EMIT_NORS_NC       default $EMIT_NORS_DAT_DIR/emit_noRS_vSmartMOM.nc
#   EMIT_NORS_JLD2     default ~/data/EMIT_MODTRANcomp/emit_modtran_noRS_results.jld2
#   EMIT_STRIDE_H2O    default 2
#   EMIT_STRIDE_AOT    default 4     (keep only AOT=0.0001)
#   EMIT_STRIDE_GNDALT default 2
#   EMIT_STRIDE_TSZ    default 1
# =============================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using JSON
using JLD2
using NCDatasets
using Dates
using Printf

const EMIT_YAML          = get(ENV, "EMIT_YAML",
                                 joinpath(pkgdir(vSmartMOM), "test", "test_parameters",
                                          "ParamsEMIT_MODTRANcomp_newLUT.yaml"))
const EMIT_JSON          = get(ENV, "EMIT_JSON",
                                 expanduser("~/data/EMIT_MODTRANcomp/emit_20250404.json"))
const EMIT_NORS_DAT_DIR  = get(ENV, "EMIT_NORS_DAT_DIR",
                                 expanduser("~/data/EMIT_MODTRANcomp/vSmartMOM_out"))
const EMIT_NORS_NC       = get(ENV, "EMIT_NORS_NC",
                                 joinpath(EMIT_NORS_DAT_DIR, "emit_noRS_vSmartMOM.nc"))
const EMIT_NORS_JLD2     = get(ENV, "EMIT_NORS_JLD2",
                                 expanduser("~/data/EMIT_MODTRANcomp/emit_modtran_noRS_results.jld2"))

const STRIDE_H2O    = parse(Int, get(ENV, "EMIT_STRIDE_H2O",    "2"))
const STRIDE_AOT    = parse(Int, get(ENV, "EMIT_STRIDE_AOT",    "4"))
const STRIDE_GNDALT = parse(Int, get(ENV, "EMIT_STRIDE_GNDALT", "2"))
const STRIDE_TSZ    = parse(Int, get(ENV, "EMIT_STRIDE_TSZ",    "1"))

"Build a LUT axis from a {min, max, spacing} dict, inclusive of max."
function lut_axis(d::AbstractDict)
    lo, hi, step = Float32(d["min"]), Float32(d["max"]), Float32(d["spacing"])
    collect(lo:step:hi)
end

p_surface_from_gndalt(gndalt_km::Real) =
    Float32(1013.25f0 * exp(-Float32(gndalt_km) / 8.0f0))

"""
    interp_at_p(profile, p_centers, p_target)

Log-pressure linear interpolation with edge clamping.
"""
function interp_at_p(profile, p_centers, p_target)
    k = searchsortedlast(p_centers, p_target)
    k < 1              && return profile[1]
    k >= length(p_centers) && return profile[end]
    lp1, lp2 = log(p_centers[k]), log(p_centers[k + 1])
    frac = (log(p_target) - lp1) / (lp2 - lp1)
    return profile[k] + frac * (profile[k + 1] - profile[k])
end

"""
    apply_scenario!(params, base_q, base_p, base_T, base_vmr; h2o, aot, gndalt)

Mutate `params` for a single LUT-grid point. VZA/VAZ are set once at the top of
the driver (paired with `tsz_axis`), so this function leaves them untouched —
all TSZ scenarios share a single rt_run via the parallel-VZA scheme.
"""
function apply_scenario!(params, base_q, base_p, base_T, base_vmr;
                         h2o::Real, aot::Real, gndalt::Real)
    p_surf = p_surface_from_gndalt(gndalt)
    params.scattering_params.rt_aerosols[1].τ_ref = aot

    k = searchsortedlast(base_p, p_surf - eps(Float32(p_surf)))
    if k >= length(base_p)
        params.p = copy(base_p); params.p[end] = p_surf
        params.T = copy(base_T); params.q = base_q .* h2o
        for (g, v) in base_vmr
            params.absorption_params.vmr[g] =
                v isa AbstractVector ? copy(v) : v
        end
    else
        new_p = vcat(base_p[1:k], p_surf)
        n_layers = k
        p_centers = [(base_p[i] + base_p[i + 1]) / 2
                     for i in 1:(length(base_p) - 1)]
        function truncate_profile(prof)
            out = similar(prof, n_layers)
            out[1:(k - 1)] .= prof[1:(k - 1)]
            out[k] = interp_at_p(prof, p_centers, p_surf)
            return out
        end
        params.p = new_p
        params.T = truncate_profile(base_T)
        params.q = truncate_profile(base_q) .* h2o
        for (g, v) in base_vmr
            params.absorption_params.vmr[g] =
                v isa AbstractVector ? truncate_profile(v) : v
        end
    end

    nothing
end

"Return (R_TOA_I, T_BOA_I) at VZA index `iv`, Stokes I."
function extract_I_radiances(R_SFI::AbstractArray, T_SFI::AbstractArray, iv::Int = 1)
    return Array(R_SFI)[iv, 1, :], Array(T_SFI)[iv, 1, :]
end

"""
    compute_total_optical_depth(τ_abs, τ_rayl, τ_aer)

Compute total column (vertical) optical depth per wavelength by summing the
per-layer τ's. Each τ array is shaped `(nSpec, nLayers)` (see the allocation
in `src/CoreRT/tools/model_from_parameters.jl`), so we reduce over dim 2.
"""
function compute_total_optical_depth(τ_abs, τ_rayl, τ_aer)
    τ_abs_vec  = vec(sum(τ_abs,  dims = 2))
    τ_rayl_vec = vec(sum(τ_rayl, dims = 2))
    τ_aer_vec  = vec(sum(τ_aer,  dims = 2))
    return τ_abs_vec .+ τ_rayl_vec .+ τ_aer_vec
end

"""
    write_scenario_dat(dir, λ_nm, ν, R_tot, T_tot, τ_total, hem_R, hem_T;
                       h2o, aot, gndalt, tsz)

Write per-scenario `.dat` file (α=0 run): TOA reflectance, BOA transmittance,
vertical total optical thickness, and bi-hemispheric flux ratios.
"""
function write_scenario_dat(dir, λ_nm, ν, R_tot_I, T_tot_I, τ_total, hem_R, hem_T;
                            h2o, aot, gndalt, tsz)
    fname = @sprintf("vSmartMOM_H2O%.4f_AOT%.4f_GNDALT%.3f_TSZ%.1f.dat",
                     h2o, aot, gndalt, tsz)
    path = joinpath(dir, fname)
    open(path, "w") do io
        @printf(io, "# vSmartMOM noRS (alb=0.0)   H2OSTR=%.4f  AOT550=%.4f  GNDALT=%.3f km  TSZ=%.1f\n",
                h2o, aot, gndalt, tsz)
        @printf(io, "# %16s %16s %16s %16s %16s %16s %16s\n",
                "wl_nm", "wn_cm-1", "R_tot", "T_tot", "tau_total", "hem_R", "hem_T")
        for i in eachindex(λ_nm)
            @printf(io, "  %16.6f %16.6f %16.8e %16.8e %16.8e %16.8e %16.8e\n",
                    λ_nm[i], ν[i], R_tot_I[i], T_tot_I[i], τ_total[i], hem_R[i], hem_T[i])
        end
    end
    return path
end

"""
    write_results_nc(path, axes..., data..., metadata)

Write aggregate LUT into NetCDF.
"""
function write_results_nc(path, h2o_axis, aot_axis, gndalt_axis, tsz_axis,
                          λ_nm, ν_axis, R_tot, T_tot, τ_total,
                          hemR_tot, hemT_tot, metadata)
    ds = NCDataset(path, "c")
    ds.attrib["title"]       = "vSmartMOM noRS radiance LUT (albedo=0.0)"
    ds.attrib["source"]      = "vSmartMOM.jl"
    ds.attrib["yaml_source"] = metadata["yaml_source"]
    ds.attrib["json_source"] = metadata["json_source"]
    ds.attrib["timestamp"]   = metadata["timestamp"]
    ds.attrib["n_total"]     = metadata["n_total"]

    defDim(ds, "H2OSTR",              length(h2o_axis))
    defDim(ds, "AOT550",              length(aot_axis))
    defDim(ds, "surface_elevation_km", length(gndalt_axis))
    defDim(ds, "TSZ",                 length(tsz_axis))
    defDim(ds, "wl",                  length(λ_nm))

    v = defVar(ds, "H2OSTR",              Float64, ("H2OSTR",));              v[:] = Float64.(h2o_axis)
    v = defVar(ds, "AOT550",              Float64, ("AOT550",));              v[:] = Float64.(aot_axis)
    v = defVar(ds, "surface_elevation_km", Float64, ("surface_elevation_km",)); v[:] = Float64.(gndalt_axis)
    v = defVar(ds, "TSZ",                 Float64, ("TSZ",));                 v[:] = Float64.(tsz_axis)
    v = defVar(ds, "wl",                  Float64, ("wl",));                  v.attrib["units"] = "nm"; v[:] = Float64.(λ_nm)
    v = defVar(ds, "wn",                  Float64, ("wl",));                  v.attrib["units"] = "cm-1"; v[:] = Float64.(ν_axis)

    data_dims = ("H2OSTR", "AOT550", "surface_elevation_km", "TSZ", "wl")
    for (name, arr, desc) in [
        ("R_tot",    R_tot,    "TOA reflected radiance, Stokes I (F0=1)"),
        ("T_tot",    T_tot,    "BOA transmitted radiance, Stokes I (F0=1)"),
        ("tau_total", τ_total, "Total vertical atmospheric optical depth (absorption + Rayleigh + aerosol extinction)"),
        ("hemR_tot", hemR_tot, "Bi-hemispheric reflectance bhr_uw (flux ratio)"),
        ("hemT_tot", hemT_tot, "Bi-hemispheric transmittance bhr_dw (flux ratio)"),
    ]
        v = defVar(ds, name, Float64, data_dims)
        v.attrib["description"] = desc
        v[:] = arr
    end

    close(ds)
    return path
end

"""
    run_all_scenarios(; yaml=EMIT_YAML, json=EMIT_JSON, ...)

Run the full LUT scan for α=0. Returns `nothing` if inputs are not available.

Port notes (vs sanghavi):
  - Parallel-VZA scheme: every TSZ in `tsz_axis` is loaded into `params.vza`
    once at the top, paired with the YAML's `vaz`. The model is built once per
    GNDALT, and each `rt_run` returns R_SFI/T_SFI of shape
    `(length(tsz_axis), n_stokes, nSpec)`. We slice per-TSZ for output.
  - Inner loop (H2OSTR × AOT550) still uses direct τ mutation for speed.
"""
function run_all_scenarios(; yaml     = EMIT_YAML,
                             json     = EMIT_JSON,
                             dat_dir  = EMIT_NORS_DAT_DIR,
                             nc       = EMIT_NORS_NC,
                             jld      = EMIT_NORS_JLD2,
                             strides  = (; h2o = STRIDE_H2O, aot = STRIDE_AOT,
                                           gndalt = STRIDE_GNDALT, tsz = STRIDE_TSZ))
    isfile(yaml) || (@info "EMIT driver: YAML not found — skipping" yaml; return nothing)
    isfile(json) || (@info "EMIT driver: JSON not found — skipping" json; return nothing)

    mkpath(dat_dir)

    @info "Loading JSON scenario descriptor" json
    cfg  = JSON.parsefile(json)
    lut  = cfg["lut_grid"]
    misc = cfg["MISC"]

    h2o_axis    = lut_axis(lut["H2OSTR"])[1:strides.h2o:end]
    aot_axis    = lut_axis(lut["AOT550"])[1:strides.aot:end]
    gndalt_axis = lut_axis(lut["GNDALT"])[1:strides.gndalt:end]
    tsz_axis    = lut_axis(lut["TSZ"])[1:strides.tsz:end]

    n_total = length(h2o_axis) * length(aot_axis) *
              length(gndalt_axis) * length(tsz_axis)
    @info "LUT grid dimensions" H2OSTR = length(h2o_axis) AOT550 = length(aot_axis) GNDALT = length(gndalt_axis) TSZ = length(tsz_axis) n_total

    params = parameters_from_yaml(yaml)
    FT     = params.float_type

    base_q, base_p, base_T = copy(params.q), copy(params.p), copy(params.T)
    base_vmr = Dict(g => v isa AbstractVector ? copy(v) : v
                    for (g, v) in params.absorption_params.vmr)

    # Parallel-VZA setup: load every TSZ as a paired (vza, vaz) pair so a single
    # rt_run yields all TSZ outputs in one shot.
    base_vaz_scalar = FT(params.vaz[1])
    params.vza = FT.(180 .- tsz_axis)
    params.vaz = fill(base_vaz_scalar, length(tsz_axis))
    @info "Parallel-VZA setup" vza = params.vza vaz = params.vaz

    n_bands = length(params.spec_bands)
    ν_axis  = vcat([collect(params.spec_bands[ib]) for ib in 1:n_bands]...)
    λ_nm    = 1e7 ./ ν_axis
    n_spec  = length(ν_axis)
    @info "Spectral grid" n_bands n_spec n_per_band = [length(collect(b)) for b in params.spec_bands]

    out_dims = (length(h2o_axis), length(aot_axis),
                length(gndalt_axis), length(tsz_axis), n_spec)
    R_tot, T_tot     = zeros(Float64, out_dims), zeros(Float64, out_dims)
    τ_total          = zeros(Float64, out_dims)
    hemR_tot, hemT_tot = zeros(Float64, out_dims), zeros(Float64, out_dims)

    scenario_k = 0
    t_start    = time()

    for (ig, gndalt) in enumerate(gndalt_axis)
        for (ih, h2o) in enumerate(h2o_axis)
            @info "Build model for (GNDALT, H2OSTR) — all TSZs in one rt_run" ig ih gndalt h2o
            # Rebuild the model per (gndalt, h2o) so the H2O scaling is
            # applied ONLY to the H2O VMR (via apply_scenario!'s `params.q
            # = base_q .* h2o`, with all other species at their base VMR).
            # The previous inner-loop trick `τ_abs .*= h2o` scaled every
            # species' contribution to the summed τ_abs (H2O + O2 + CO2 +
            # NO2 + SO2 + CH4 + O3 + N2O + CO), which was wrong.
            apply_scenario!(params, base_q, base_p, base_T, base_vmr;
                            h2o = h2o, aot = aot_axis[1], gndalt = gndalt)
            model = model_from_parameters(params)
            RS    = InelasticScattering.noRS{FT}()

            # Aerosol τ scales linearly with τ_ref, so AOT can still be
            # scaled in the inner loop. τ_abs and τ_rayl are h2o-correct
            # as built and must NOT be rescaled here.
            τ_aer_base = [copy(model.τ_aer[ib]) for ib in 1:n_bands]
            aot_base   = aot_axis[1]

            for (ia, aot) in enumerate(aot_axis)
                t_scen = time()

                # Per-TSZ accumulators for the full concatenated spectrum.
                n_tsz = length(tsz_axis)
                R_concat = [Float64[] for _ in 1:n_tsz]
                T_concat = [Float64[] for _ in 1:n_tsz]
                τ_concat    = Float64[]
                hemR_concat = Float64[]
                hemT_concat = Float64[]

                for ib in 1:n_bands
                    model.τ_aer[ib] .= τ_aer_base[ib] .* (aot / aot_base)
                    # τ_abs[ib] and τ_rayl[ib] are h2o-correct from the
                    # model build above — leave them as-is.

                    # rt_run returns 7-tuple (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw).
                    # R_SFI/T_SFI shape: (length(vza), n_stokes, nSpec_band) — one slice per TSZ.
                    R_SFI_b, T_SFI_b, _, _, _, hem_R_b, hem_T_b = rt_run(RS, model, ib)

                    for it in 1:n_tsz
                        R_band, T_band = extract_I_radiances(R_SFI_b, T_SFI_b, it)
                        append!(R_concat[it], R_band)
                        append!(T_concat[it], T_band)
                    end

                    τ_b = compute_total_optical_depth(model.τ_abs[ib],
                                                       model.τ_rayl[ib],
                                                       model.τ_aer[ib])
                    append!(τ_concat, τ_b)
                    append!(hemR_concat, Float64.(hem_R_b))
                    append!(hemT_concat, Float64.(hem_T_b))
                end

                for (it, tsz) in enumerate(tsz_axis)
                    scenario_k += 1
                    R_tot[ih, ia, ig, it, :]    .= R_concat[it]
                    T_tot[ih, ia, ig, it, :]    .= T_concat[it]
                    τ_total[ih, ia, ig, it, :]  .= τ_concat
                    hemR_tot[ih, ia, ig, it, :] .= hemR_concat
                    hemT_tot[ih, ia, ig, it, :] .= hemT_concat

                    write_scenario_dat(dat_dir, λ_nm, ν_axis,
                                       R_concat[it], T_concat[it], τ_concat,
                                       hemR_concat, hemT_concat;
                                       h2o = h2o, aot = aot, gndalt = gndalt, tsz = tsz)
                end

                dt_scen   = round(time() - t_scen;   digits = 2)
                elapsed_s = round(time() - t_start; digits = 1)
                @info "Inner-loop point complete (all bands × all TSZs)" scenario_k n_total dt_s = dt_scen elapsed_s h2o aot gndalt
            end
        end
    end

    metadata = Dict(
        "json_source"    => json,
        "yaml_source"    => yaml,
        "timestamp"      => string(Dates.now()),
        "misc_from_json" => misc,
        "n_total"        => n_total,
        "stride"         => (; H2O = strides.h2o, AOT = strides.aot,
                               GNDALT = strides.gndalt, TSZ = strides.tsz),
        "surface_albedo" => 0.0,
    )

    @info "Saving JLD2" jld
    @save jld h2o_axis aot_axis gndalt_axis tsz_axis ν_axis λ_nm R_tot T_tot τ_total hemR_tot hemT_tot metadata
    @info "Writing NetCDF" nc
    write_results_nc(nc, h2o_axis, aot_axis, gndalt_axis, tsz_axis,
                     λ_nm, ν_axis, R_tot, T_tot, τ_total,
                     hemR_tot, hemT_tot, metadata)

    return (; h2o_axis, aot_axis, gndalt_axis, tsz_axis,
              ν_axis, λ_nm, R_tot, T_tot, τ_total,
              hemR_tot, hemT_tot, metadata)
end

if abspath(PROGRAM_FILE) == @__FILE__
    res = run_all_scenarios()
    isnothing(res) && @info "emit_modtran_noRS_scenarios: no outputs produced (missing inputs)."
end
