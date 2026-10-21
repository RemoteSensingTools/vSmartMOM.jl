using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using DelimitedFiles
using Dates
using Printf
using Statistics

const OUTDIR = get(ENV, "OUTDIR",
    joinpath(homedir(), "Raman_misc", "workflows", "diagnostics",
             "current_vs_june2024_flat_noabs_761_764_cpu_" *
             Dates.format(now(), "yyyymmdd_HHMMSS")))
const YAML_PATH = joinpath(OUTDIR, "june2024_flat_noabs_cpu.yaml")
const VERSION_TAG = "june2024"
const ALBEDOS = [0.0, 0.3, 1.0]
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

fmt_float(x) = replace(@sprintf("%.1f", x), "." => "p")

function yaml_array(xs)
    return "[" * join(string.(xs), ", ") * "]"
end

function write_yaml()
    mkpath(OUTDIR)
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

geometry:
  sza: 19.0
  vza: [0.0]
  vaz: [0.0]
  obs_alt: [0]

atmospheric_profile:
  T: $(yaml_array(T_PROFILE))
  p: $(yaml_array(P_PROFILE))
  profile_reduction: 1
""")
    end
end

function make_rrs(FT, n_pol, ν, model)
    ν̃ = mean(ν)
    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)
    F0 = zeros(FT, n_pol, length(ν))
    F0[1, :] .= one(FT)
    SIF0 = zeros(FT, n_pol, length(ν))
    rs = InelasticScattering.RRS(
        n2 = n2,
        o2 = o2,
        greek_raman = InelasticScattering.GreekCoefs(
            FT[1], FT[1], FT[1], FT[1], FT[1], FT[1]),
        fscattRayl = FT[1],
        ϖ_Cabannes = FT[1],
        ϖ_λ₁λ₀ = zeros(FT, 1),
        i_λ₁λ₀ = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1),
        Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
        i_ref = argmin(abs.(ν .- ν̃)),
        n_Raman = 0,
        iBand = [1],
        F₀ = F0,
        SIF₀ = SIF0,
    )
    CoreRT.getRamanSSProp!(rs, 1e7 / ν̃, ν)
    return rs, effT
end

function metadata_lines(params, model, rs, ν, effT)
    keep = findall(x -> OUTPUT_WN_MIN <= x <= OUTPUT_WN_MAX, ν)
    i_ref = argmin(abs.(ν .- mean(ν)))
    return [
        "version=$VERSION_TAG",
        "project=$(Base.active_project())",
        "git_commit=$(readchomp(`git rev-parse HEAD`))",
        "yaml=$YAML_PATH",
        "spectral_grid_solve_cm-1=$(first(ν)):$(ν[2]-ν[1]):$(last(ν))",
        "spectral_grid_output_cm-1=$(ν[first(keep)]):$(ν[2]-ν[1]):$(ν[last(keep)])",
        "n_solve=$(length(ν))",
        "n_output=$(length(keep))",
        "n_layers=$(length(model.profile.p_full))",
        "max_m_count=$(params.max_m)",
        "m_loop_orders=0:$(params.max_m - 1)",
        "rho_ray_yaml=$(params.depol)",
        "rho_ray_model_greek=$(params.depol)",
        "rho_cab_model_greek=$(params.depol)",
        "rho_cab_internal=not_separate_in_this_commit",
        "varpi_cab_rs=$(rs.ϖ_Cabannes[1])",
        "rrs_weight_sum_before_rt=$(sum(rs.ϖ_λ₁λ₀))",
        "n_raman=$(rs.n_Raman)",
        "tau_rayl_center_column=$(sum(model.τ_rayl[1][i_ref, :]))",
        "tau_rayl_output_min=$(minimum(sum(model.τ_rayl[1][keep, :], dims=2)))",
        "tau_rayl_output_max=$(maximum(sum(model.τ_rayl[1][keep, :], dims=2)))",
        "effective_temperature=$effT",
    ]
end

function run_case(albedo)
    params = parameters_from_yaml(YAML_PATH)
    FT = params.float_type
    params.brdf[1] = CoreRT.LambertianSurfaceScalar{FT}(FT(albedo))

    model = model_from_parameters(params)

    ν = model.params.spec_bands[1]
    keep = findall(x -> OUTPUT_WN_MIN <= x <= OUTPUT_WN_MAX, ν)
    n_pol = model.params.polarization_type.n
    F0 = zeros(FT, n_pol, length(ν))
    F0[1, :] .= one(FT)
    SIF0 = zeros(FT, n_pol, length(ν))

    rs, effT = make_rrs(FT, n_pol, ν, model)
    rrs_result = CoreRT.rt_run_test(rs, model, 1)

    nors = InelasticScattering.noRS{FT}(
        fscattRayl = FT[1],
        ϖ_Cabannes = FT[1],
        bandSpecLim = [],
        iBand = [1],
        F₀ = F0,
        SIF₀ = SIF0,
    )
    nors_result = CoreRT.rt_run_test(nors, model, 1)

    R = Array(rrs_result[1])
    ieR = Array(rrs_result[3])
    R0 = Array(nors_result[1])
    F = ones(FT, length(keep))

    alb = fmt_float(albedo)
    rrs_path = joinpath(OUTDIR, "$(VERSION_TAG)_flat_noabs_alb$(alb)_rrs.dat")
    nors_path = joinpath(OUTDIR, "$(VERSION_TAG)_flat_noabs_alb$(alb)_nors.dat")
    writedlm(rrs_path, hcat(ν[keep], R[1, 1, keep], R[1, 2, keep], R[1, 3, keep],
                            ieR[1, 1, keep], ieR[1, 2, keep], ieR[1, 3, keep], F))
    writedlm(nors_path, hcat(ν[keep], R0[1, 1, keep], R0[1, 2, keep], R0[1, 3, keep], F))

    open(joinpath(OUTDIR, "$(VERSION_TAG)_flat_noabs_alb$(alb)_metadata.txt"), "w") do io
        for line in metadata_lines(params, model, rs, ν, effT)
            println(io, line)
        end
        println(io, "albedo=$albedo")
        println(io, "rrs_weight_sum_after_rt=$(sum(rs.ϖ_λ₁λ₀))")
        println(io, "note=June-2024 commit has one Rayleigh Greek matrix from YAML depol; Cabannes does not have a separate depol path here.")
    end
    @info "wrote june2024 flat/noabs case" albedo rrs_path nors_path
end

write_yaml()
for albedo in ALBEDOS
    run_case(albedo)
end
println(OUTDIR)
