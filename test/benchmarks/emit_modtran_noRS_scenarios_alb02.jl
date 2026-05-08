# =============================================================================
# emit_modtran_noRS_scenarios_alb02.jl — noRS forward LUT driver with α=0.2
#
# Paired companion to emit_modtran_noRS_scenarios.jl. Writes exactly the same
# per-scenario fields (R_tot, T_tot, tau_total, hemR_tot, hemT_tot) so that
# modtran_equivalent_fields.jl can assemble the MODTRAN-equivalent fields via
# the two-albedo Lambertian closure (see plans/MODTRAN_equivalent_equations.md).
#
# 2026-04-24 rewrite:
#   * Default YAML switched to ParamsEMIT_MODTRANcomp_newLUT_alb02.yaml
#     (rebuilt HITRAN LUTs under ~/data/HITRAN_LUTs/).
#   * AOT stride defaults to 4 (AOT=0.0001 only) to match the α=0 driver.
#
# Output files: vSmartMOM_alb02_H2O*_AOT*_GNDALT*_TSZ*.dat
# =============================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using JSON
using JLD2
using NCDatasets
using Dates
using Printf

const EMIT_YAML_ALB02     = get(ENV, "EMIT_YAML_ALB02",
                                 joinpath(pkgdir(vSmartMOM), "test", "test_parameters",
                                          "ParamsEMIT_MODTRANcomp_newLUT_alb02.yaml"))
const EMIT_JSON          = get(ENV, "EMIT_JSON",
                                 expanduser("~/data/EMIT_MODTRANcomp/emit_20250404.json"))
const EMIT_NORS_DAT_DIR  = get(ENV, "EMIT_NORS_DAT_DIR",
                                 expanduser("~/data/EMIT_MODTRANcomp/vSmartMOM_out"))
const EMIT_NORS_NC_ALB02 = get(ENV, "EMIT_NORS_NC_ALB02",
                                 joinpath(EMIT_NORS_DAT_DIR, "emit_noRS_vSmartMOM_alb02.nc"))
const EMIT_NORS_JLD2_ALB02 = get(ENV, "EMIT_NORS_JLD2_ALB02",
                                 expanduser("~/data/EMIT_MODTRANcomp/emit_modtran_noRS_results_alb02.jld2"))

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
    apply_scenario!(params, base_q, base_p, base_T, base_vmr; h2o, aot, gndalt, tsz)

Mutate `params` for a single LUT-grid point.
"""
function apply_scenario!(params, base_q, base_p, base_T, base_vmr;
                         h2o::Real, aot::Real, gndalt::Real, tsz::Real)
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

    params.vza = [180 - tsz]
    nothing
end

"Return (R_TOA_I, T_BOA_I) as spectra for the first VZA, Stokes I."
function extract_I_radiances(R_SFI::AbstractArray, T_SFI::AbstractArray)
    return Array(R_SFI)[1, 1, :], Array(T_SFI)[1, 1, :]
end

"""
    compute_total_optical_depth(τ_abs, τ_rayl, τ_aer)

Compute total column optical depth per wavelength by summing across layers.
τ_abs/τ_rayl/τ_aer are each (nSpec, nLayers) (see `τ_abs` alloc in
`src/CoreRT/tools/model_from_parameters.jl`), so we reduce over dim 2.
"""
function compute_total_optical_depth(τ_abs, τ_rayl, τ_aer)
    τ_abs_vec  = vec(sum(τ_abs,  dims = 2))
    τ_rayl_vec = vec(sum(τ_rayl, dims = 2))
    τ_aer_vec  = vec(sum(τ_aer,  dims = 2))
    return τ_abs_vec .+ τ_rayl_vec .+ τ_aer_vec
end

"""
    write_scenario_dat(dir, λ_nm, ν, R_tot, T_tot, τ_total; h2o, aot, gndalt, tsz)

Write scenario output with total optical depth instead of single-scattering.
"""
function write_scenario_dat(dir, λ_nm, ν, R_tot_I, T_tot_I, τ_total, hem_R, hem_T;
                            h2o, aot, gndalt, tsz)
    fname = @sprintf("vSmartMOM_alb02_H2O%.4f_AOT%.4f_GNDALT%.3f_TSZ%.1f.dat",
                     h2o, aot, gndalt, tsz)
    path = joinpath(dir, fname)
    open(path, "w") do io
        @printf(io, "# vSmartMOM noRS (alb02)   H2OSTR=%.4f  AOT550=%.4f  GNDALT=%.3f km  TSZ=%.1f\n",
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

Write aggregate LUT into NetCDF with total optical depth.
"""
function write_results_nc(path, h2o_axis, aot_axis, gndalt_axis, tsz_axis,
                          λ_nm, ν_axis, R_tot, T_tot, τ_total,
                          hemR_tot, hemT_tot, metadata)
    ds = NCDataset(path, "c")
    ds.attrib["title"]       = "vSmartMOM noRS radiance LUT (albedo=0.2)"
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
        ("R_tot",    R_tot,    "TOA reflected radiance (multiple scattering)"),
        ("T_tot",    T_tot,    "BOA transmitted radiance (multiple scattering)"),
        ("tau_total", τ_total, "Total atmospheric optical depth (Rayleigh + aerosol + absorption)"),
        ("hemR_tot", hemR_tot, "Hemispheric reflectance (multiple scattering)"),
        ("hemT_tot", hemT_tot, "Hemispheric transmittance (multiple scattering)"),
    ]
        v = defVar(ds, name, Float64, data_dims)
        v.attrib["description"] = desc
        v[:] = arr
    end

    close(ds)
    return path
end

"""
    run_all_scenarios(; yaml=EMIT_YAML_ALB02, json=EMIT_JSON, ...)

Run the full LUT scan with albedo 0.2 and optical depth output.
"""
function run_all_scenarios(; yaml     = EMIT_YAML_ALB02,
                             json     = EMIT_JSON,
                             dat_dir  = EMIT_NORS_DAT_DIR,
                             nc       = EMIT_NORS_NC_ALB02,
                             jld      = EMIT_NORS_JLD2_ALB02,
                             strides  = (; h2o = STRIDE_H2O, aot = STRIDE_AOT,
                                           gndalt = STRIDE_GNDALT, tsz = STRIDE_TSZ))
    isfile(yaml) || (@info "EMIT driver (alb02): YAML not found — skipping" yaml; return nothing)
    isfile(json) || (@info "EMIT driver (alb02): JSON not found — skipping" json; return nothing)

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

    ν_axis = collect(params.spec_bands[1])
    λ_nm   = 1e7 ./ ν_axis
    n_spec = length(ν_axis)

    out_dims = (length(h2o_axis), length(aot_axis),
                length(gndalt_axis), length(tsz_axis), n_spec)
    R_tot, T_tot     = zeros(Float64, out_dims), zeros(Float64, out_dims)
    τ_total          = zeros(Float64, out_dims)
    hemR_tot, hemT_tot = zeros(Float64, out_dims), zeros(Float64, out_dims)

    scenario_k = 0
    t_start    = time()

    for (ig, gndalt) in enumerate(gndalt_axis)
        for (it, tsz) in enumerate(tsz_axis)
            for (ih, h2o) in enumerate(h2o_axis)
                @info "Build model for (GNDALT, TSZ, H2OSTR)" ig it ih gndalt tsz h2o
                # Rebuild the model per (gndalt, tsz, h2o) so the H2O
                # scaling is applied ONLY to the H2O VMR (via
                # apply_scenario!'s `params.q = base_q .* h2o`, with all
                # other species at their base VMR). The previous inner-loop
                # trick `τ_abs .*= h2o` scaled every species' contribution
                # to the summed τ_abs.
                apply_scenario!(params, base_q, base_p, base_T, base_vmr;
                                h2o = h2o, aot = aot_axis[1], gndalt = gndalt, tsz = tsz)
                model = model_from_parameters(params)
                RS    = InelasticScattering.noRS{FT}()

                # Aerosol τ still scales linearly with τ_ref, so AOT scaling
                # remains in the inner loop. τ_abs/τ_rayl are h2o-correct
                # from the build above and must NOT be rescaled.
                τ_aer_base = copy(model.τ_aer[1])
                aot_base   = aot_axis[1]

                for (ia, aot) in enumerate(aot_axis)
                    scenario_k += 1
                    t_scen = time()

                    model.τ_aer[1] .= τ_aer_base .* (aot / aot_base)
                    # τ_abs[1] and τ_rayl[1] left as built (h2o-correct).

                    # rt_run returns 7-tuple (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw).
                    R_SFI_tot, T_SFI_tot, _, _, _, hem_R_tot, hem_T_tot = rt_run(RS, model, 1)

                    R_tot_I, T_tot_I = extract_I_radiances(R_SFI_tot, T_SFI_tot)

                    # Compute total optical depth
                    τ_total_spec = compute_total_optical_depth(model.τ_abs[1],
                                                               model.τ_rayl[1],
                                                               model.τ_aer[1])

                    R_tot[ih, ia, ig, it, :] .= R_tot_I
                    T_tot[ih, ia, ig, it, :] .= T_tot_I
                    τ_total[ih, ia, ig, it, :] .= τ_total_spec
                    hemR_tot[ih, ia, ig, it, :] .= Float64.(hem_R_tot)
                    hemT_tot[ih, ia, ig, it, :] .= Float64.(hem_T_tot)

                    write_scenario_dat(dat_dir, λ_nm, ν_axis,
                                       R_tot_I, T_tot_I, τ_total_spec,
                                       hem_R_tot, hem_T_tot;
                                       h2o = h2o, aot = aot, gndalt = gndalt, tsz = tsz)

                    dt_scen   = round(time() - t_scen;   digits = 2)
                    elapsed_s = round(time() - t_start; digits = 1)
                    @info "Scenario complete" scenario_k n_total dt_s = dt_scen elapsed_s h2o aot gndalt tsz
                end
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
        "surface_albedo" => 0.2,
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
    isnothing(res) && @info "emit_modtran_noRS_scenarios_alb02: no outputs produced (missing inputs)."
end
