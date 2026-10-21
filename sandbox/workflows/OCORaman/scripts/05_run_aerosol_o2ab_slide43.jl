#!/usr/bin/env julia

# Recreate the aerosol/Rayleigh-error ASCII products behind the old
# evalAer_O2ABbands_Raman.jl slide figure.
#
# Output columns:
#   rrs:  wn, Icab, Qcab, Ucab, IRRS, QRRS, URRS, F0
#   nors: wn, IRayl, QRayl, URayl, F0
#
# Default case set matches the uploaded slide screenshot:
#   SZA=0, VZA=0, VAZ=0, psurf=1000 hPa
#   rho = 0, 0.1, 0.5, 1.0
#   AOD = 0, 0.1, 0.5
#   aerosol scene for nonzero AOD: z0=2, r0=0.1

using CUDA
using Dates
using DelimitedFiles
using Distributions
using Interpolations
using Printf
using Statistics

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel

const OUTDIR = get(ENV, "OUTDIR", "/home/sanghavi/data/RamanSIFgrid/aerosol_test_new")
const OLD_ROOT = get(ENV, "OLD_SMARTMOM_ROOT",
    "/home/sanghavi/Raman_misc/workflows/worktrees/uni_vSmartMOM_sanghavi_2025-03-18")
const OLD_PARAM_DIR = joinpath(OLD_ROOT, "test", "test_parameters")
const SOLAR_OUT = get(ENV, "SOLAR_OUT", joinpath(OLD_ROOT, "src", "SolarModel", "solar.out"))
const O2_ABSCO = get(ENV, "O2_ABSCO", "/net/fluo/data1/ABSCO_CS_Database/v5.2_final/o2_v52_v2.jld2")
const OVERLAP_NU = parse(Float64, get(ENV, "OVERLAP_NU", "200"))
const RUN_TAGS = split(get(ENV, "BANDS", "ABO2,BBO2"), ","; keepempty=false)
const SMOKE = get(ENV, "SMOKE", "0") == "1"

mkpath(OUTDIR)
mkpath(joinpath(OUTDIR, "config"))
mkpath(joinpath(OUTDIR, "logs"))

if CUDA.functional()
    device!(parse(Int, get(ENV, "CUDA_DEVICE", "0")))
    CUDA.allowscalar(false)
else
    error("CUDA is not functional; refusing to start expensive aerosol run without GPU.")
end

function patched_yaml(src_name::AbstractString)
    src = joinpath(OLD_PARAM_DIR, src_name)
    dst = joinpath(OUTDIR, "config", replace(src_name, ".yaml" => "_current_schema.yaml"))
    text = read(src, String)
    text = replace(text, "GaussQuadHemisphere()" => "GaussLegQuad()")
    text = replace(text,
        "/net/fluo/data3/data/Databases/ABSCO_CS_Database/v5.2_final/o2_v52_extra.jld2" => O2_ABSCO)
    # Current vSmartMOM treats nonnegative `depol` as an explicit override for
    # both Rayleigh and Cabannes Greek coefficients. The old sanghavi Raman
    # YAMLs carried `depol: 0.028` as a mostly dummy/default value while the
    # Raman code computed the actual Rayleigh/Cabannes depolarizations
    # internally. Preserve that sanghavi behavior by using the current
    # molecular auto mode.
    text = replace(text, r"(?m)^(\s*)depol:\s*[-+0-9.eE]+.*$" => s"\1depol:              -1")
    write(dst, text)
    return dst
end

const BAND_CONFIGS = Dict(
    "ABO2" => patched_yaml("testAer_O2Aband_Raman.yaml"),
    "BBO2" => patched_yaml("testAer_O2Bband_Raman.yaml"),
)

function old_float_label(x::Real)
    if isinteger(x)
        return string(Int(x))
    end
    s = @sprintf("%.1f", Float64(x))
    return replace(s, "." => "p")
end

function old_albedo_label(x::Real)
    s = @sprintf("%.1f", Float64(x))
    return replace(s, "." => "p")
end

function with_lambertian_albedo(model::CoreRT.RTModel, alpha)
    FT = CoreRT.float_type(model)
    surfaces = [CoreRT.LambertianSurfaceScalar{FT}(FT(alpha)) for _ in 1:length(CoreRT.get_surfaces(model))]
    return CoreRT.RTModel(
        model.architecture,
        model.solver,
        model.numerics,
        model.geometry,
        model.quad_points,
        model.atmosphere,
        model.optics,
        surfaces,
        model.sources,
    )
end

function solar_interpolator(FT)
    if isfile(SOLAR_OUT)
        solar = solar_transmission_from_file(SOLAR_OUT)
        @info "Using solar transmission" SOLAR_OUT size=size(solar)
        return LinearInterpolation(FT.(solar[4:end, 1]), FT.(solar[4:end, 2]),
                                   extrapolation_bc = Line())
    end
    @warn "Solar transmission file missing; using unit solar transmission" SOLAR_OUT
    return x -> one(FT)
end

function empty_rrs(FT, n_pol, n_spec)
    zmat = zeros(FT, n_pol, n_spec)
    return InelasticScattering.RRS{FT}(
        n2 = nothing,
        o2 = nothing,
        greek_raman = InelasticScattering.GreekCoefs(
            FT[1], FT[1], FT[1], FT[1], FT[1], FT[1]),
        fscattRayl = FT[1],
        ϖ_Cabannes = FT[1],
        ϖ_λ₁λ₀ = zeros(FT, 1),
        i_λ₁λ₀ = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1),
        Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
        i_ref = 0,
        n_Raman = 0,
        F₀ = copy(zmat),
        SIF₀ = copy(zmat),
    )
end

function make_rrs(FT, n_pol, ν, ν_center, effT, F0, cabannes_by_band)
    n_model_bands = length(cabannes_by_band)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν_center, effT)
    rs = InelasticScattering.RRS(
        n2 = n2,
        o2 = o2,
        greek_raman = InelasticScattering.GreekCoefs(
            FT[1], FT[1], FT[1], FT[1], FT[1], FT[1]),
        fscattRayl = ones(FT, n_model_bands),
        ϖ_Cabannes = ones(FT, n_model_bands),
        ϖ_λ₁λ₀ = zeros(FT, 1),
        i_λ₁λ₀ = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1),
        Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
        i_ref = argmin(abs.(ν .- mean(ν))),
        n_Raman = 0,
        iBand = [1],
        F₀ = copy(F0),
        SIF₀ = zeros(FT, n_pol, length(ν)),
    )
    CoreRT.getRamanSSProp!(rs, 1e7 / mean(ν), ν)
    rs.ϖ_Cabannes .= FT.(cabannes_by_band)
    return rs
end

function make_nors(FT, n_pol, F0, n_model_bands)
    return InelasticScattering.noRS{FT}(
        fscattRayl = FT[1],
        ϖ_Cabannes = ones(FT, n_model_bands),
        bandSpecLim = UnitRange{Int}[],
        iBand = [1],
        F₀ = copy(F0),
        SIF₀ = zeros(FT, n_pol, size(F0, 2)),
    )
end

function checked_writedlm(path, data)
    mkpath(dirname(path))
    if any(!isfinite, data)
        error("Refusing to write non-finite data to $path")
    end
    writedlm(path, data)
end

function case_paths(scene, albedo, band_tag)
    sza_str = old_float_label(scene.sza)
    vza_str = "0p0"
    vaz_str = "0p0"
    alb_str = old_albedo_label(albedo)
    psurf_str = "1000"

    if scene.aod == 0
        stem = "aer_AOD0p0_sza$(sza_str)_vza$(vza_str)_vaz$(vaz_str)_alb$(alb_str)_psurf$(psurf_str)hpa"
    else
        aod_str = old_albedo_label(scene.aod)
        z_str = old_float_label(scene.z0)
        r_str = old_albedo_label(scene.r0)
        stem = "aer_AOD$(aod_str)_z$(z_str)_r$(r_str)_sza$(sza_str)_vza$(vza_str)_vaz$(vaz_str)_alb$(alb_str)_psurf$(psurf_str)hpa"
    end

    rrs_path = joinpath(OUTDIR, "$(stem)_rrs_$(band_tag).dat")
    nors_path = joinpath(OUTDIR, "$(stem)_nors_$(band_tag).dat")
    return rrs_path, nors_path
end

function valid_dat(path, nrows, ncols)
    isfile(path) || return false
    try
        data = readdlm(path)
        return size(data) == (nrows, ncols) && all(isfinite, data)
    catch
        return false
    end
end

function configure_scene!(params, scene)
    FT = params.float_type
    params.sza = FT(scene.sza)
    params.vza = FT[0]
    params.vaz = FT[0]

    aer = params.scattering_params.rt_aerosols[1]
    aer.τ_ref = FT(scene.aod)
    if scene.aod > 0
        aer.profile = LogNormal(log(FT(scene.z0)), FT(0.4))
        aer.aerosol.size_distribution = LogNormal(log(FT(scene.r0)), FT(1.12))
    end

    return params
end

function run_band_case!(base_model, band_tag, scene, albedo, solar_T)
    FT = CoreRT.float_type(base_model)
    model = with_lambertian_albedo(base_model, albedo)

    n_bands = length(CoreRT.get_spec_bands(model))
    Δν = CoreRT.get_spec_bands(model)[1][2] - CoreRT.get_spec_bands(model)[1][1]
    n_overlap = Int(floor(FT(OVERLAP_NU) / Δν))
    tot_nspec = sum(length.(CoreRT.get_spec_bands(model))) - 2 * n_bands * n_overlap
    tot_nspec > 0 || error("Non-positive output spectrum length for $band_tag")
    rrs_path, nors_path = case_paths(scene, albedo, band_tag)
    if get(ENV, "FORCE", "0") != "1" &&
       valid_dat(rrs_path, tot_nspec, 8) &&
       valid_dat(nors_path, tot_nspec, 5)
        @info "Skipping existing valid outputs" rrs_path nors_path rows=tot_nspec
        return nothing
    end

    n_cam = length(model.obs_geom.vza)
    n_pol = CoreRT.polarization_type(model).n
    R = zeros(FT, n_cam, n_pol, tot_nspec)
    RnoRS = zeros(FT, n_cam, n_pol, tot_nspec)
    ieR = zeros(FT, n_cam, n_pol, tot_nspec)
    tot_ν = zeros(FT, tot_nspec)
    tot_F0 = zeros(FT, tot_nspec)

    ν_center = FT(0.5) * (CoreRT.get_spec_bands(model)[1][1] + CoreRT.get_spec_bands(model)[end][end])
    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)

    spec_end = 0
    for i_band in 1:n_bands
        ν = CoreRT.get_spec_bands(model)[i_band]
        spec_start = spec_end + 1
        spec_end += length(ν) - 2 * n_overlap
        keep = (n_overlap + 1):(length(ν) - n_overlap)

        P = FT.(planck_spectrum_wn(FT(5777), collect(ν)) .* FT(2.1629e-05 * π))
        F_vec = FT.(solar_T.(ν)) .* P
        F0 = zeros(FT, n_pol, length(ν))
        F0[1, :] .= F_vec
        sources = CoreRT.SolarBeam(F₀ = F0) + CoreRT.SurfaceSIF(SIF₀ = zeros(FT, n_pol, length(ν)))

        cabannes_by_band = collect(FT, model.ϖ_Cabannes)
        rs_rrs = make_rrs(FT, n_pol, ν, ν_center, effT, F0, cabannes_by_band)
        rs_nors = make_nors(FT, n_pol, F0, length(cabannes_by_band))

        @info "RT" band_tag i_band n_bands aod=scene.aod albedo=albedo nSpec=length(ν)
        rrs_result = CoreRT.rt_run_test(rs_rrs, model, i_band; sources)
        nors_result = CoreRT.rt_run_test(rs_nors, model, i_band; sources)
        length(rrs_result) == 7 || error("Unexpected RRS result tuple length: $(length(rrs_result))")
        length(nors_result) == 7 || error("Unexpected noRS result tuple length: $(length(nors_result))")
        R1 = rrs_result[1]
        ieR1 = rrs_result[3]
        R0 = nors_result[1]

        tot_ν[spec_start:spec_end] .= ν[keep]
        tot_F0[spec_start:spec_end] .= F_vec[keep]
        R[:, :, spec_start:spec_end] .= R1[:, :, keep]
        RnoRS[:, :, spec_start:spec_end] .= R0[:, :, keep]
        ieR[:, :, spec_start:spec_end] .= ieR1[:, :, keep]
    end

    for vctr in 1:n_cam
        rrs = [tot_ν R[vctr, 1, :] R[vctr, 2, :] R[vctr, 3, :] ieR[vctr, 1, :] ieR[vctr, 2, :] ieR[vctr, 3, :] tot_F0]
        nors = [tot_ν RnoRS[vctr, 1, :] RnoRS[vctr, 2, :] RnoRS[vctr, 3, :] tot_F0]
        checked_writedlm(rrs_path, rrs)
        checked_writedlm(nors_path, nors)
        @info "Wrote" rrs_path nors_path rows=size(rrs, 1)
    end

    return nothing
end

function run!()
    albedos = SMOKE ? [0.0] : [0.0, 0.1, 0.5, 1.0]
    scenes = [
        (; aod = 0.5, z0 = 2.0, r0 = 0.1, sza = 0.0),
        (; aod = 0.1, z0 = 2.0, r0 = 0.1, sza = 0.0),
        (; aod = 0.0, z0 = 2.0, r0 = 0.1, sza = 0.0),
    ]
    if SMOKE
        scenes = scenes[1:1]
    end

    manifest = joinpath(OUTDIR, "run_manifest_$(Dates.format(now(), "yyyymmdd_HHMMSS")).txt")
    open(manifest, "w") do io
        println(io, "created = $(now())")
        println(io, "script = $(abspath(PROGRAM_FILE))")
        println(io, "outdir = $OUTDIR")
        println(io, "old_root = $OLD_ROOT")
        println(io, "solar_out = $SOLAR_OUT")
        println(io, "o2_absco = $O2_ABSCO")
        println(io, "bands = $(join(RUN_TAGS, ","))")
        println(io, "smoke = $SMOKE")
        println(io, "overlap_nu = $OVERLAP_NU")
        println(io, "albedos = $(join(albedos, ","))")
        println(io, "scenes = $(scenes)")
    end
    @info "Manifest" manifest

    band_params = Dict{String, Any}()
    for band_tag in RUN_TAGS
        yaml_path = BAND_CONFIGS[band_tag]
        params = parameters_from_yaml(yaml_path)
        FT = params.float_type
        @assert FT === Float64
        @assert params.architecture isa vSmartMOM.Architectures.GPU
        band_params[band_tag] = params
    end

    solar_T = solar_interpolator(Float64)
    for scene in scenes
        for band_tag in RUN_TAGS
            params = band_params[band_tag]
            scene_params = deepcopy(params)
            configure_scene!(scene_params, scene)
            @info "Building model" band_tag aod=scene.aod z0=scene.z0 r0=scene.r0 sza=scene.sza
            base_model = model_from_parameters(scene_params)
            for albedo in albedos
                run_band_case!(base_model, band_tag, scene, albedo, solar_T)
                GC.gc()
                CUDA.reclaim()
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run!()
end
