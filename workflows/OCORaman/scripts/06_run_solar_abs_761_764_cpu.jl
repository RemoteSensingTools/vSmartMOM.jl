using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel
using DelimitedFiles
using Dates
using Interpolations
using Printf
using Statistics

const VERSION_TAG = get(ENV, "VERSION_TAG", "current")
const CASE_TAG = "solar_abs"
const OUTDIR = get(ENV, "OUTDIR",
    joinpath(homedir(), "Raman_misc", "workflows", "diagnostics",
             "current_vs_june2024_$(CASE_TAG)_761_764_cpu_" *
             Dates.format(now(), "yyyymmdd_HHMMSS")))
const YAML_PATH = joinpath(OUTDIR, "$(VERSION_TAG)_$(CASE_TAG)_cpu.yaml")
const OLD_SOLAR_OUT = "/home/sanghavi/Raman_misc/workflows/worktrees/uni_vSmartMOM_sanghavi_2024-06-21/src/SolarModel/solar.out"
const SOLAR_OUT = get(ENV, "SOLAR_OUT", OLD_SOLAR_OUT)
const O2_ABSCO = get(ENV, "O2_ABSCO", "/net/fluo/data1/ABSCO_CS_Database/v5.2_final/o2_v52_v2.jld2")
const SOLVE_WN_MIN = 1e7 / (764 + 10)
const SOLVE_WN_MAX = 1e7 / (761 - 10)
const OUTPUT_WN_MIN = 1e7 / 764
const OUTPUT_WN_MAX = 1e7 / 761

const T_PROFILE = [
    231.62, 244.33, 251.34, 258.09, 264.25, 269.15,
    272.59, 274.07, 273.30, 269.65, 264.27, 258.11,
    251.52, 245.22, 239.20, 234.05, 229.71, 225.70,
    222.70, 220.62, 219.32, 217.93, 216.98, 217.10,
    218.35, 223.33, 234.19, 249.34, 264.12, 277.20,
    280.77, 282.60, 284.40, 285.80,
]

const P_PROFILE = [
      0.14,   0.22,   0.30,   0.39,   0.53,   0.71,
      0.96,   1.28,   1.70,   2.27,   3.03,   4.03,
      5.44,   7.26,   9.67,  12.90,  17.23,  23.30,
     31.00,  42.07,  56.09,  74.78,  99.69, 131.00,
    176.85, 236.64, 314.58, 418.87, 557.76, 735.00,
    800.12, 849.00, 912.00, 980.00, 1000.0,
]

function parse_albedos(text)
    vals = Float64[]
    for part in split(text, ","; keepempty=false)
        push!(vals, parse(Float64, strip(part)))
    end
    return vals
end

const ALBEDOS = parse_albedos(get(ENV, "ALBEDOS", "0.0,0.3,1.0"))

fmt_float(x) = replace(@sprintf("%.1f", x), "." => "p")

function yaml_array(xs)
    return "[" * join(string.(xs), ", ") * "]"
end

function write_yaml()
    mkpath(OUTDIR)
    if VERSION_TAG == "june2024"
        open(YAML_PATH, "w") do io
            println(io, """
radiative_transfer:
  spec_bands:
    - $(SOLVE_WN_MIN):0.1:$(SOLVE_WN_MAX)
  surface:
    - LambertianSurfaceScalar{Float64}(0.0)
  quadrature_type: GaussQuadHemisphere()
  polarization_type: Stokes_IQU()
  max_m: 3
  Δ_angle: 2.0
  l_trunc: 5
  depol: 0.028
  float_type: Float64
  architecture: Architectures.CPU()

geometry:
  sza: 19.0
  vza: [0.0]
  vaz: [0.0]
  obs_alt: [0]

atmospheric_profile:
  T: $(yaml_array(T_PROFILE))
  p: $(yaml_array(P_PROFILE))
  profile_reduction: 1

absorption:
  molecules:
    - [O2]
  LUTfiles:
    - ["$(O2_ABSCO)"]
  vmr:
    O2: 0.21
  broadening: Voigt()
  CEF: HumlicekWeidemann32SDErrorFunction()
  wing_cutoff: 150

scattering:
  aerosols:
    - τ_ref:          0.0
      μ:              0.01
      σ:              1.12
      nᵣ:             1.5
      nᵢ:             0.00000001
      z₀:             2.0
      σ₀:             0.4
  r_max:          50.0
  nquad_radius:   2500
  λ_ref:          0.770
  decomp_type:    NAI2()
""")
        end
    else
        open(YAML_PATH, "w") do io
            println(io, """
radiative_transfer:
  spec_bands:
    - $(SOLVE_WN_MIN):0.1:$(SOLVE_WN_MAX)
  surface:
    - LambertianSurfaceScalar{Float64}(0.0)
  polarization_type: Stokes_IQU()
  nstreams: 3
  max_m: 3
  truncation: auto
  depol: -1
  float_type: Float64
  architecture: Architectures.CPU()

geometry:
  sza: 19.0
  vza: [0.0]
  vaz: [0.0]
  obs_alt: [0]

atmospheric_profile:
  T: $(yaml_array(T_PROFILE))
  p: $(yaml_array(P_PROFILE))
  profile_reduction: 1

absorption:
  fixed_molecules:
    - [O2]
  variable_molecules:
    - []
  LUTfiles:
    - ["$(O2_ABSCO)"]
  vmr:
    O2: 0.21
  broadening: Voigt()
  CEF: HumlicekWeidemann32SDErrorFunction()
  wing_cutoff: 150

scattering:
  aerosols:
    - τ_ref:          0.0
      μ:              0.01
      σ:              1.12
      nᵣ:             1.5
      nᵢ:             0.00000001
      z₀:             2.0
      σ₀:             0.4
  r_max:          50.0
  nquad_radius:   2500
  λ_ref:          0.770
  decomp_type:    NAI2()
""")
        end
    end
end

function with_lambertian_albedo(model, albedo)
    FT = VERSION_TAG == "june2024" ? model.params.float_type : CoreRT.float_type(model)
    if VERSION_TAG == "june2024"
        model.params.brdf[1] = CoreRT.LambertianSurfaceScalar{FT}(FT(albedo))
        return model
    end
    surfaces = [CoreRT.LambertianSurfaceScalar{FT}(FT(albedo)) for _ in 1:length(CoreRT.get_surfaces(model))]
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

function spec_band(model)
    return VERSION_TAG == "june2024" ? model.params.spec_bands[1] : CoreRT.get_spec_bands(model)[1]
end

function n_pol(model)
    return VERSION_TAG == "june2024" ? model.params.polarization_type.n : CoreRT.polarization_type(model).n
end

function solar_interpolator(FT)
    isfile(SOLAR_OUT) || error("SOLAR_OUT not found: $SOLAR_OUT")
    solar = solar_transmission_from_file(SOLAR_OUT)
    @info "Using solar transmission" SOLAR_OUT size=size(solar)
    return LinearInterpolation(FT.(solar[4:end, 1]), FT.(solar[4:end, 2]), extrapolation_bc = Line())
end

function solar_F0(FT, nstokes, ν, solar_T)
    P = FT.(planck_spectrum_wn(FT(5777), collect(ν)) .* FT(2.1629e-05 * π))
    Fvec = FT.(solar_T.(ν)) .* P
    F0 = zeros(FT, nstokes, length(ν))
    F0[1, :] .= Fvec
    return F0
end

function make_rrs(FT, nstokes, ν, model, F0)
    ν̃ = mean(ν)
    i_ref = argmin(abs.(ν .- ν̃))
    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)
    cab = VERSION_TAG == "june2024" ? FT[1] : FT.(copy(model.ϖ_Cabannes))
    rs = InelasticScattering.RRS(
        n2 = n2,
        o2 = o2,
        greek_raman = InelasticScattering.GreekCoefs(
            FT[1], FT[1], FT[1], FT[1], FT[1], FT[1]),
        fscattRayl = FT[1],
        ϖ_Cabannes = cab,
        ϖ_λ₁λ₀ = zeros(FT, 1),
        i_λ₁λ₀ = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1),
        Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
        i_ref = i_ref,
        n_Raman = 0,
        iBand = [1],
        F₀ = copy(F0),
        SIF₀ = zeros(FT, nstokes, length(ν)),
    )
    CoreRT.getRamanSSProp!(rs, 1e7 / ν̃, ν)
    VERSION_TAG == "current" && (rs.ϖ_Cabannes .= FT.(model.ϖ_Cabannes))
    return rs, effT
end

function make_nors(FT, nstokes, ν, model, F0)
    if VERSION_TAG == "june2024"
        return InelasticScattering.noRS{FT}(
            fscattRayl = FT[1],
            ϖ_Cabannes = FT[1],
            bandSpecLim = [],
            iBand = [1],
            F₀ = copy(F0),
            SIF₀ = zeros(FT, nstokes, length(ν)),
        )
    end
    return InelasticScattering.noRS{FT}(
        fscattRayl = FT[1],
        ϖ_Cabannes = ones(FT, length(model.ϖ_Cabannes)),
        bandSpecLim = UnitRange{Int}[],
        iBand = [1],
        F₀ = copy(F0),
        SIF₀ = zeros(FT, nstokes, length(ν)),
    )
end

function run_rrs(rs, model, F0)
    if VERSION_TAG == "june2024"
        return CoreRT.rt_run_test(rs, model, 1)
    end
    src = CoreRT.SolarBeam(F₀ = F0) +
          CoreRT.SurfaceSIF(SIF₀ = zeros(eltype(F0), size(F0, 1), size(F0, 2)))
    return CoreRT.rt_run_test(rs, model, [1]; sources = src)
end

function metadata_lines(params, model, rs, ν, effT)
    keep = findall(x -> OUTPUT_WN_MIN <= x <= OUTPUT_WN_MAX, ν)
    i_ref = argmin(abs.(ν .- mean(ν)))
    project_dir = dirname(Base.active_project())
    git_commit = readchomp(`git -C $project_dir rev-parse HEAD`)
    base = [
        "version=$VERSION_TAG",
        "case=$CASE_TAG",
        "project=$(Base.active_project())",
        "git_commit=$git_commit",
        "yaml=$YAML_PATH",
        "solar_out=$SOLAR_OUT",
        "o2_absco=$O2_ABSCO",
        "absorption=O2_ABSCO_only",
        "spectral_grid_solve_cm-1=$(first(ν)):$(ν[2]-ν[1]):$(last(ν))",
        "spectral_grid_output_cm-1=$(ν[first(keep)]):$(ν[2]-ν[1]):$(ν[last(keep)])",
        "n_solve=$(length(ν))",
        "n_output=$(length(keep))",
        "n_layers=$(length(model.profile.p_full))",
        "rrs_weight_sum_before_rt=$(sum(rs.ϖ_λ₁λ₀))",
        "n_raman=$(rs.n_Raman)",
        "tau_rayl_center_column=$(sum(model.τ_rayl[1][i_ref, :]))",
        "tau_rayl_output_min=$(minimum(sum(model.τ_rayl[1][keep, :], dims=2)))",
        "tau_rayl_output_max=$(maximum(sum(model.τ_rayl[1][keep, :], dims=2)))",
        "tau_abs_center_column=$(sum(model.τ_abs[1][i_ref, :]))",
        "tau_abs_output_min=$(minimum(sum(model.τ_abs[1][keep, :], dims=2)))",
        "tau_abs_output_max=$(maximum(sum(model.τ_abs[1][keep, :], dims=2)))",
        "effective_temperature=$effT",
    ]
    if VERSION_TAG == "june2024"
        append!(base, [
            "max_m_count=$(params.max_m)",
            "m_loop_orders=0:$(params.max_m - 1)",
            "depol_config=$(params.depol)",
            "rho_ray_yaml=$(params.depol)",
            "rho_ray_model_greek=$(params.depol)",
            "rho_cab_model_greek=$(params.depol)",
            "rho_cab_internal=not_separate_in_this_commit",
            "varpi_cab_rs=$(rs.ϖ_Cabannes[1])",
            "note=June-2024 commit has one Rayleigh Greek matrix from YAML depol; Cabannes does not have a separate depol path here.",
        ])
    else
        λ_center = 1e7 / mean(ν)
        n2, o2 = InelasticScattering.getRamanAtmoConstants(mean(ν), effT)
        γ_cab, _ = InelasticScattering.compute_γ_air_Cabannes!(λ_center, n2, o2)
        γ_ray, _ = InelasticScattering.compute_γ_air_Rayleigh!(λ_center, n2, o2)
        append!(base, [
            "n_fourier_moments=$(CoreRT.n_fourier_moments_bands(model)[1])",
            "m_max_order=$(CoreRT.m_max_bands(model)[1])",
            "depol_config=$(params.depol)",
            "rho_cab_internal=$(2γ_cab / (1 + γ_cab))",
            "rho_ray_internal=$(2γ_ray / (1 + γ_ray))",
            "varpi_cab_model=$(model.ϖ_Cabannes[1])",
            "varpi_cab_rs=$(rs.ϖ_Cabannes[1])",
        ])
    end
    return base
end

function run_case(base_model, params, solar_T, albedo)
    model = with_lambertian_albedo(base_model, albedo)
    FT = VERSION_TAG == "june2024" ? params.float_type : CoreRT.float_type(model)
    ν = spec_band(model)
    keep = findall(x -> OUTPUT_WN_MIN <= x <= OUTPUT_WN_MAX, ν)
    F0 = solar_F0(FT, n_pol(model), ν, solar_T)

    rs, effT = make_rrs(FT, n_pol(model), ν, model, F0)
    rrs_result = run_rrs(rs, model, F0)
    nors = make_nors(FT, n_pol(model), ν, model, F0)
    nors_result = run_rrs(nors, model, F0)

    R = Array(rrs_result[1])
    ieR = Array(rrs_result[3])
    R0 = Array(nors_result[1])

    alb = fmt_float(albedo)
    rrs_path = joinpath(OUTDIR, "$(VERSION_TAG)_$(CASE_TAG)_alb$(alb)_rrs.dat")
    nors_path = joinpath(OUTDIR, "$(VERSION_TAG)_$(CASE_TAG)_alb$(alb)_nors.dat")
    writedlm(rrs_path, hcat(ν[keep], R[1, 1, keep], R[1, 2, keep], R[1, 3, keep],
                            ieR[1, 1, keep], ieR[1, 2, keep], ieR[1, 3, keep], F0[1, keep]))
    writedlm(nors_path, hcat(ν[keep], R0[1, 1, keep], R0[1, 2, keep], R0[1, 3, keep], F0[1, keep]))

    open(joinpath(OUTDIR, "$(VERSION_TAG)_$(CASE_TAG)_alb$(alb)_metadata.txt"), "w") do io
        for line in metadata_lines(params, model, rs, ν, effT)
            println(io, line)
        end
        println(io, "albedo=$albedo")
        println(io, "rrs_weight_sum_after_rt=$(sum(rs.ϖ_λ₁λ₀))")
        println(io, "F0_output_min=$(minimum(F0[1, keep]))")
        println(io, "F0_output_max=$(maximum(F0[1, keep]))")
    end
    @info "wrote $(VERSION_TAG) $(CASE_TAG) case" albedo rrs_path nors_path
end

write_yaml()
params = parameters_from_yaml(YAML_PATH)
solar_T = solar_interpolator(params.float_type)
base_model = model_from_parameters(params)
for albedo in ALBEDOS
    run_case(base_model, params, solar_T, albedo)
end
println(OUTDIR)
