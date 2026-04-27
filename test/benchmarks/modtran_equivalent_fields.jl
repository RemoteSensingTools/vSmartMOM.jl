# =============================================================================
# modtran_equivalent_fields.jl — Generate MODTRAN-equivalent atmospheric fields
# from paired vSmartMOM runs at surface albedo 0.0 and 0.2.
#
# Uses the 6S/Vermote two-albedo Lambertian decomposition:
#
#     ρ_TOA(ρ_s) = ρ_atm + T↓(μ_s) · ρ_s · T↑(μ_v) / (1 − ρ_s · S)
#
# The full derivation, including the flux-based closure for S, is in
# plans/MODTRAN_equivalent_equations.md.
#
# Output fields (one per (H2O, AOT, GNDALT, TSZ) scenario × wavelength):
#   rhoatm           — path reflectance  (π·L₀/μ_s with ρ_s=0)
#   transm_down_dir  — direct downwelling transmittance exp(−τ/μ_s)
#   transm_down_dif  — diffuse downwelling transmittance T↓(μ_s) − exp(−τ/μ_s)
#   transm_up_dir    — direct upwelling transmittance exp(−τ/μ_v)
#   transm_up_dif    — diffuse upwelling transmittance T↑(μ_v) − exp(−τ/μ_v)
#   sphalb           — atmospheric spherical albedo S
#
# 2026-04-24 rewrite: sphalb is now computed from the two-albedo flux closure
#
#     S = (1 − bhr_dw(0) / bhr_dw(ρ_s2)) / ρ_s2
#
# instead of the single-scatter proxy hemR_ss/μ_s used previously. This gives
# the true multiple-scattering spherical albedo and no longer requires an SS
# run in the α=0 driver.
#
# Inputs (both produced by the emit_modtran_noRS_scenarios*.jl drivers):
#   EMIT_NORS_JLD2        ~/data/EMIT_MODTRANcomp/emit_modtran_noRS_results.jld2
#   EMIT_NORS_JLD2_ALB02  ~/data/EMIT_MODTRANcomp/emit_modtran_noRS_results_alb02.jld2
#
# Outputs:
#   NetCDF LUT aggregate          MODTRAN_EQUIV_NC  (vSmartMOM_out/modtran_equiv.nc)
#   Per-scenario .dat files       vSmartMOM_modtran_equiv_H2O*_AOT*_GNDALT*_TSZ*.dat
# =============================================================================

using vSmartMOM
using JLD2
using NCDatasets
using Printf
using Dates: now

const EMIT_NORS_JLD2        = get(ENV, "EMIT_NORS_JLD2",
    expanduser("~/data/EMIT_MODTRANcomp/emit_modtran_noRS_results.jld2"))
const EMIT_NORS_JLD2_ALB02  = get(ENV, "EMIT_NORS_JLD2_ALB02",
    expanduser("~/data/EMIT_MODTRANcomp/emit_modtran_noRS_results_alb02.jld2"))
const MODTRAN_EQUIV_NC      = get(ENV, "MODTRAN_EQUIV_NC",
    expanduser("~/data/EMIT_MODTRANcomp/vSmartMOM_out/modtran_equiv.nc"))
const MODTRAN_EQUIV_DAT_DIR = get(ENV, "MODTRAN_EQUIV_DAT_DIR",
    dirname(MODTRAN_EQUIV_NC))

const RHO_S2 = 0.2  # reference Lambertian albedo for the second run

clamp_nonneg!(a) = @. a = max(a, 0)

"""
    write_modtran_equiv_dat(dir, λ_nm, ν, fields; h2o, aot, gndalt, tsz)

Write MODTRAN-equivalent output to a `.dat` file with the same column layout
as the reference files in `~/data/EMIT_MODTRANcomp/MODTRAN_out/`.
"""
function write_modtran_equiv_dat(dir, λ_nm, ν, fields; h2o, aot, gndalt, tsz)
    rhoatm, transm_down_dir, transm_down_dif,
        transm_up_dir, transm_up_dif, sphalb = fields

    fname = @sprintf("vSmartMOM_modtran_equiv_H2O%.4f_AOT%.4f_GNDALT%.3f_TSZ%.1f.dat",
                     h2o, aot, gndalt, tsz)
    path = joinpath(dir, fname)
    open(path, "w") do io
        @printf(io, "# vSmartMOM MODTRAN-equivalent  H2OSTR=%.4f  AOT550=%.4f  GNDALT=%.3f km  TSZ=%.1f\n",
                h2o, aot, gndalt, tsz)
        @printf(io, "# %16s %16s %16s %16s %16s %16s %16s %16s\n",
                "wl_nm", "wn_cm-1", "rhoatm",
                "transm_down_dir", "transm_down_dif",
                "transm_up_dir", "transm_up_dif", "sphalb")
        for i in eachindex(λ_nm)
            @printf(io, "  %16.6f %16.6f %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
                    λ_nm[i], ν[i], rhoatm[i],
                    transm_down_dir[i], transm_down_dif[i],
                    transm_up_dir[i], transm_up_dif[i], sphalb[i])
        end
    end
    return path
end

"""
    write_modtran_equiv_nc(path, axes, fields, metadata)

Aggregate NetCDF writer (all scenarios × wavelengths).
"""
function write_modtran_equiv_nc(path, h2o_axis, aot_axis, gndalt_axis, tsz_axis,
                                λ_nm, ν_axis, rhoatm,
                                transm_down_dir, transm_down_dif,
                                transm_up_dir, transm_up_dif, sphalb, metadata)
    ds = NCDataset(path, "c")
    ds.attrib["title"]       = "vSmartMOM MODTRAN-equivalent atmospheric LUT"
    ds.attrib["description"] = "Lambertian two-albedo decomposition of paired vSmartMOM noRS runs (ρ_s=0 and ρ_s=0.2). See plans/MODTRAN_equivalent_equations.md."
    ds.attrib["source"]      = "vSmartMOM.jl"
    ds.attrib["alb0_jld2"]   = metadata["alb0_jld2"]
    ds.attrib["alb02_jld2"]  = metadata["alb02_jld2"]
    ds.attrib["yaml_source"] = metadata["yaml_source"]
    ds.attrib["timestamp"]   = metadata["timestamp"]
    ds.attrib["sza_deg"]     = metadata["sza_deg"]
    ds.attrib["mu_s"]        = metadata["mu_s"]
    ds.attrib["rho_s2"]      = metadata["rho_s2"]

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

    dims = ("H2OSTR", "AOT550", "surface_elevation_km", "TSZ", "wl")
    for (name, arr, desc) in [
        ("rhoatm",          rhoatm,          "Path reflectance at sensor VZA (π·L_TOA/μ_s with ρ_s=0)"),
        ("transm_down_dir", transm_down_dir, "Direct downwelling transmittance exp(-τ/μ_s)"),
        ("transm_down_dif", transm_down_dif, "Diffuse downwelling transmittance"),
        ("transm_up_dir",   transm_up_dir,   "Direct upwelling transmittance exp(-τ/μ_v)"),
        ("transm_up_dif",   transm_up_dif,   "Diffuse upwelling transmittance (Lambertian-to-directional)"),
        ("sphalb",          sphalb,          "Atmospheric spherical albedo S (two-albedo flux closure)"),
    ]
        var = defVar(ds, name, Float64, dims)
        var.attrib["description"] = desc
        var[:] = arr
    end

    close(ds)
    return path
end

"""
    generate_modtran_equivalent_fields()

Load paired α=0 and α=0.2 JLD2 results and produce MODTRAN-equivalent fields
via the two-albedo Lambertian decomposition. Writes a per-scenario `.dat` file
and an aggregate NetCDF.
"""
function generate_modtran_equivalent_fields()
    isfile(EMIT_NORS_JLD2) ||
        (@info "MODTRAN-equivalent: α=0 JLD2 missing — skipping" EMIT_NORS_JLD2; return nothing)
    isfile(EMIT_NORS_JLD2_ALB02) ||
        (@info "MODTRAN-equivalent: α=0.2 JLD2 missing — skipping" EMIT_NORS_JLD2_ALB02; return nothing)

    # ---------- Load α=0 run ----------
    @info "Loading α=0 JLD2" EMIT_NORS_JLD2
    alb0 = JLD2.jldopen(EMIT_NORS_JLD2, "r") do io
        (; h2o_axis    = io["h2o_axis"],
           aot_axis    = io["aot_axis"],
           gndalt_axis = io["gndalt_axis"],
           tsz_axis    = io["tsz_axis"],
           ν_axis      = io["ν_axis"],
           λ_nm        = io["λ_nm"],
           R_tot       = io["R_tot"],
           τ_total     = io["τ_total"],
           hemT_tot    = io["hemT_tot"],
           metadata    = io["metadata"])
    end

    # ---------- Load α=0.2 run ----------
    @info "Loading α=0.2 JLD2" EMIT_NORS_JLD2_ALB02
    alb02 = JLD2.jldopen(EMIT_NORS_JLD2_ALB02, "r") do io
        (; R_tot    = io["R_tot"],
           τ_total  = io["τ_total"],
           hemT_tot = io["hemT_tot"],
           metadata = io["metadata"])
    end

    @assert size(alb0.R_tot)   == size(alb02.R_tot)    "α=0 / α=0.2 R_tot shapes differ"
    @assert size(alb0.R_tot)   == size(alb02.τ_total)  "R_tot / τ_total shapes differ"
    @assert size(alb0.hemT_tot) == size(alb02.hemT_tot) "hemT_tot shapes differ"

    # ---------- Load YAML for SZA ----------
    yaml = alb0.metadata["yaml_source"]
    isfile(yaml) || (@info "YAML not found — skipping" yaml; return nothing)
    params = parameters_from_yaml(yaml)
    sza    = Float64(params.sza)
    μ_s    = cos(deg2rad(sza))
    @info "Geometry" sza=sza μ_s=μ_s

    # TSZ → VZA: driver sets params.vza = [180 - tsz]; μ_v = cos(VZA).
    μ_v_axis = [cos(deg2rad(180.0 - Float64(tsz))) for tsz in alb0.tsz_axis]

    # ---------- Allocate output arrays ----------
    dims = size(alb0.R_tot)  # (nH2O, nAOT, nGNDALT, nTSZ, nλ)
    rhoatm          = zeros(Float64, dims)
    transm_down_dir = zeros(Float64, dims)
    transm_down_dif = zeros(Float64, dims)
    transm_up_dir   = zeros(Float64, dims)
    transm_up_dif   = zeros(Float64, dims)
    sphalb          = zeros(Float64, dims)

    # ---------- Per-scenario decomposition ----------
    nH2O, nAOT, nGND, nTSZ, nλ = dims
    @info "Computing MODTRAN-equivalent fields" nH2O nAOT nGND nTSZ nλ
    for it in 1:nTSZ
        μ_v = μ_v_axis[it]
        for ig in 1:nGND, ia in 1:nAOT, ih in 1:nH2O
            τ     = @view alb0.τ_total[ih, ia, ig, it, :]
            L0    = @view alb0.R_tot[ih, ia, ig, it, :]   # radiance, F₀=1
            L2    = @view alb02.R_tot[ih, ia, ig, it, :]
            hT0   = @view alb0.hemT_tot[ih, ia, ig, it, :]   # bhr_dw(α=0)
            hT2   = @view alb02.hemT_tot[ih, ia, ig, it, :]  # bhr_dw(α=0.2)

            rho     = @view rhoatm[ih, ia, ig, it, :]
            tdw_dir = @view transm_down_dir[ih, ia, ig, it, :]
            tdw_dif = @view transm_down_dif[ih, ia, ig, it, :]
            tup_dir = @view transm_up_dir[ih, ia, ig, it, :]
            tup_dif = @view transm_up_dif[ih, ia, ig, it, :]
            S       = @view sphalb[ih, ia, ig, it, :]

            # Path reflectance (radiance → reflectance).
            @. rho = π * L0 / μ_s

            # Spherical albedo from two-albedo flux closure:
            #   bhr_dw(ρ_s) = bhr_dw(0) / (1 − ρ_s · S)
            # → S = (1 − bhr_dw(0)/bhr_dw(ρ_s)) / ρ_s
            # In deep line cores both fluxes collapse to ~0 and the ratio
            # is 0/0 → NaN. Fall back to S=0 when bhr_dw(ρ_s) is below a
            # tiny threshold — there's no light to bounce back anyway.
            @. S = ifelse(hT2 > eps(Float64), (1 - hT0 / hT2) / RHO_S2, 0.0)

            # Direct-beam transmittances (Beer-Lambert).
            @. tdw_dir = exp(-τ / μ_s)
            @. tup_dir = exp(-τ / μ_v)

            # Total downward transmittance (diffuse + direct) from the α=0 BOA
            # flux: bhr_dw(0) = μ_s · T↓(μ_s).
            T_down_tot = hT0 ./ μ_s
            @. tdw_dif = T_down_tot - tdw_dir

            # Two-albedo closure for total upward transmittance:
            #   π·(L2 − L0)/μ_s = ρ_s · T↓(μ_s) · T↑(μ_v) / (1 − ρ_s·S)
            # → T↑(μ_v) = π·(L2 − L0)·(1 − ρ_s·S) / (ρ_s · bhr_dw(0))
            # Same guard: when bhr_dw(0) ≈ 0, no direct-beam photons reach
            # BOA so T↑ is unresolved; set T↑ = exp(-τ/μ_v) so diffuse = 0.
            α_ρ      = @. π * (L2 - L0) / μ_s
            T_up_tot = @. ifelse(hT0 > eps(Float64),
                                  α_ρ * (1 - RHO_S2 * S) / (RHO_S2 * T_down_tot),
                                  tup_dir)
            @. tup_dif = T_up_tot - tup_dir
        end
    end

    clamp_nonneg!(transm_down_dif)
    clamp_nonneg!(transm_up_dif)
    clamp_nonneg!(sphalb)

    # ---------- Persist outputs ----------
    mkpath(dirname(MODTRAN_EQUIV_NC))
    metadata = Dict(
        "alb0_jld2"   => EMIT_NORS_JLD2,
        "alb02_jld2"  => EMIT_NORS_JLD2_ALB02,
        "yaml_source" => yaml,
        "timestamp"   => string(now()),
        "sza_deg"     => sza,
        "mu_s"        => μ_s,
        "rho_s2"      => RHO_S2,
    )

    @info "Writing aggregate NetCDF" MODTRAN_EQUIV_NC
    write_modtran_equiv_nc(MODTRAN_EQUIV_NC,
        alb0.h2o_axis, alb0.aot_axis, alb0.gndalt_axis, alb0.tsz_axis,
        alb0.λ_nm, alb0.ν_axis,
        rhoatm, transm_down_dir, transm_down_dif,
        transm_up_dir, transm_up_dif, sphalb, metadata)

    @info "Writing per-scenario .dat files to" MODTRAN_EQUIV_DAT_DIR
    mkpath(MODTRAN_EQUIV_DAT_DIR)
    n = 0
    for (it, tsz) in enumerate(alb0.tsz_axis),
        (ig, gndalt) in enumerate(alb0.gndalt_axis),
        (ia, aot) in enumerate(alb0.aot_axis),
        (ih, h2o) in enumerate(alb0.h2o_axis)

        fields = (rhoatm[ih, ia, ig, it, :],
                  transm_down_dir[ih, ia, ig, it, :],
                  transm_down_dif[ih, ia, ig, it, :],
                  transm_up_dir[ih, ia, ig, it, :],
                  transm_up_dif[ih, ia, ig, it, :],
                  sphalb[ih, ia, ig, it, :])
        write_modtran_equiv_dat(MODTRAN_EQUIV_DAT_DIR, alb0.λ_nm, alb0.ν_axis,
                                fields;
                                h2o = h2o, aot = aot, gndalt = gndalt, tsz = tsz)
        n += 1
        n % 50 == 0 && @info "Wrote scenario $n / $(prod(dims[1:4]))"
    end
    @info "Complete" scenarios = n

    return (; rhoatm, transm_down_dir, transm_down_dif,
              transm_up_dir, transm_up_dif, sphalb, metadata)
end

if abspath(PROGRAM_FILE) == @__FILE__
    res = generate_modtran_equivalent_fields()
    isnothing(res) && @info "modtran_equivalent_fields: no outputs produced (missing inputs)."
end
