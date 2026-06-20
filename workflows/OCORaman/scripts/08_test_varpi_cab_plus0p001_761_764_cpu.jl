using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel
using DelimitedFiles
using Dates
using Interpolations
using LaTeXStrings
using Plots
using Printf
using Statistics

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

const VERSION_TAG = "current"
const VARPI_CAB_OFFSET = parse(Float64, get(ENV, "VARPI_CAB_OFFSET", "0.001"))
const VARPI_CAB_SET_TEXT = strip(get(ENV, "VARPI_CAB_SET", ""))
const USE_VARPI_CAB_SET = !isempty(VARPI_CAB_SET_TEXT)
const VARPI_CAB_SET_VALUE = USE_VARPI_CAB_SET ? parse(Float64, VARPI_CAB_SET_TEXT) : NaN
function offset_tag(x::Real)
    s = @sprintf("%.6f", x)
    s = replace(s, r"0+$" => "")
    s = replace(s, r"\.$" => "")
    return replace(s, "-" => "m", "." => "p", "+" => "p")
end
offset_label(x::Real) = @sprintf("%.6g", x)
const VARPI_CAB_TAG = (
    USE_VARPI_CAB_SET ?
    "varpiCab_set$(offset_tag(VARPI_CAB_SET_VALUE))" :
    "varpiCab_plus$(offset_tag(VARPI_CAB_OFFSET))"
)
const VARPI_CAB_LABEL = (
    USE_VARPI_CAB_SET ?
    "ϖ_Cabannes = $(offset_label(VARPI_CAB_SET_VALUE))" :
    "ϖ_Cabannes + $(offset_label(VARPI_CAB_OFFSET))"
)
const VARPI_CAB_MODE = USE_VARPI_CAB_SET ? "set" : "offset"
const VARPI_CAB_SET_METADATA = USE_VARPI_CAB_SET ? string(VARPI_CAB_SET_VALUE) : "none"
const CASE_TAG = "solar_abs_$(VARPI_CAB_TAG)"
const BASELINE_DIR = get(ENV, "BASELINE_DIR",
    joinpath(homedir(), "Raman_misc", "workflows", "diagnostics",
             "current_vs_june2024_solar_abs_761_764_cpu_20260614"))
const OUTDIR = get(ENV, "OUTDIR",
    joinpath(homedir(), "Raman_misc", "workflows", "diagnostics",
             "$(CASE_TAG)_761_764_cpu_" *
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
  obs_alt: 0

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

function with_lambertian_albedo(model, albedo)
    FT = CoreRT.float_type(model)
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

spec_band(model) = CoreRT.get_spec_bands(model)[1]

n_pol(model) = CoreRT.polarization_type(model).n

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

function make_rrs(FT, nstokes, ν, model, F0; varpi_cab_offset = 0.0, varpi_cab_set = nothing)
    ν̃ = mean(ν)
    i_ref = argmin(abs.(ν .- ν̃))
    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)
    cab = FT.(copy(model.ϖ_Cabannes))
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
    # Diagnostic-only perturbation: Raman redistribution weights are computed
    # first; only the Cabannes elastic fraction used by the RRS elastic path is
    # changed afterward.
    if varpi_cab_set === nothing
        rs.ϖ_Cabannes .= FT.(model.ϖ_Cabannes) .+ FT(varpi_cab_offset)
    else
        rs.ϖ_Cabannes .= FT(varpi_cab_set)
    end
    return rs, effT
end

function make_nors(FT, nstokes, ν, model, F0)
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
    src = CoreRT.SolarBeam(F₀ = F0) +
          CoreRT.SurfaceSIF(SIF₀ = zeros(eltype(F0), size(F0, 1), size(F0, 2)))
    return CoreRT.rt_run_test(rs, model, [1]; sources = src)
end

function metadata_lines(params, model, rs_baseline, rs_perturbed, ν, effT)
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
        "rrs_weight_sum_baseline=$(sum(rs_baseline.ϖ_λ₁λ₀))",
        "rrs_weight_sum_perturbed=$(sum(rs_perturbed.ϖ_λ₁λ₀))",
        "n_raman=$(rs_perturbed.n_Raman)",
        "tau_rayl_center_column=$(sum(model.τ_rayl[1][i_ref, :]))",
        "tau_rayl_output_min=$(minimum(sum(model.τ_rayl[1][keep, :], dims=2)))",
        "tau_rayl_output_max=$(maximum(sum(model.τ_rayl[1][keep, :], dims=2)))",
        "tau_abs_center_column=$(sum(model.τ_abs[1][i_ref, :]))",
        "tau_abs_output_min=$(minimum(sum(model.τ_abs[1][keep, :], dims=2)))",
        "tau_abs_output_max=$(maximum(sum(model.τ_abs[1][keep, :], dims=2)))",
        "effective_temperature=$effT",
    ]
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
        "varpi_cab_baseline_rs=$(rs_baseline.ϖ_Cabannes[1])",
        "varpi_cab_mode=$VARPI_CAB_MODE",
        "varpi_cab_offset=$(VARPI_CAB_OFFSET)",
        "varpi_cab_set=$(VARPI_CAB_SET_METADATA)",
        "varpi_cab_perturbed_rs=$(rs_perturbed.ϖ_Cabannes[1])",
        "diagnostic_note=Only RS_type.ϖ_Cabannes is changed after getRamanSSProp!; Raman weights are left at their computed baseline values.",
    ])
    return base
end

function run_case(base_model, params, solar_T, albedo)
    model = with_lambertian_albedo(base_model, albedo)
    FT = CoreRT.float_type(model)
    ν = spec_band(model)
    keep = findall(x -> OUTPUT_WN_MIN <= x <= OUTPUT_WN_MAX, ν)
    F0 = solar_F0(FT, n_pol(model), ν, solar_T)

    rs_baseline, effT = make_rrs(FT, n_pol(model), ν, model, F0; varpi_cab_offset = 0.0)
    rs_perturbed, _ = make_rrs(FT, n_pol(model), ν, model, F0;
        varpi_cab_offset = VARPI_CAB_OFFSET,
        varpi_cab_set = USE_VARPI_CAB_SET ? VARPI_CAB_SET_VALUE : nothing)
    perturbed_result = run_rrs(rs_perturbed, model, F0)

    R_pert = Array(perturbed_result[1])
    ieR_pert = Array(perturbed_result[3])

    alb = fmt_float(albedo)
    baseline_src = joinpath(BASELINE_DIR, "current_solar_abs_alb$(alb)_rrs.dat")
    nors_src = joinpath(BASELINE_DIR, "current_solar_abs_alb$(alb)_nors.dat")
    isfile(baseline_src) || error("Baseline RRS file not found: $baseline_src")
    isfile(nors_src) || error("Baseline noRS file not found: $nors_src")

    base_path = joinpath(OUTDIR, "baseline_$(CASE_TAG)_alb$(alb)_rrs.dat")
    pert_path = joinpath(OUTDIR, "perturbed_$(CASE_TAG)_alb$(alb)_rrs.dat")
    nors_path = joinpath(OUTDIR, "$(CASE_TAG)_alb$(alb)_nors.dat")
    cp(baseline_src, base_path; force = true)
    cp(nors_src, nors_path; force = true)
    writedlm(pert_path, hcat(ν[keep], R_pert[1, 1, keep], R_pert[1, 2, keep], R_pert[1, 3, keep],
                             ieR_pert[1, 1, keep], ieR_pert[1, 2, keep], ieR_pert[1, 3, keep], F0[1, keep]))

    open(joinpath(OUTDIR, "$(CASE_TAG)_alb$(alb)_metadata.txt"), "w") do io
        for line in metadata_lines(params, model, rs_baseline, rs_perturbed, ν, effT)
            println(io, line)
        end
        println(io, "albedo=$albedo")
        println(io, "baseline_dir=$BASELINE_DIR")
        println(io, "baseline_rrs_source=$baseline_src")
        println(io, "baseline_nors_source=$nors_src")
        println(io, "rrs_weight_sum_after_rt_baseline=$(sum(rs_baseline.ϖ_λ₁λ₀))")
        println(io, "rrs_weight_sum_after_rt_perturbed=$(sum(rs_perturbed.ϖ_λ₁λ₀))")
        println(io, "F0_output_min=$(minimum(F0[1, keep]))")
        println(io, "F0_output_max=$(maximum(F0[1, keep]))")
    end
    @info "wrote diagnostic $(CASE_TAG) case" albedo base_path pert_path nors_path
end

function read_rgb(mode::AbstractString, albedo::Float64)
    alb = fmt_float(albedo)
    rrs_path = mode == "baseline" ?
        joinpath(OUTDIR, "baseline_$(CASE_TAG)_alb$(alb)_rrs.dat") :
        joinpath(OUTDIR, "perturbed_$(CASE_TAG)_alb$(alb)_rrs.dat")
    nors_path = joinpath(OUTDIR, "$(CASE_TAG)_alb$(alb)_nors.dat")
    rrs = Float64.(readdlm(rrs_path))
    nors = Float64.(readdlm(nors_path))
    maximum(abs.(rrs[:, 1] .- nors[:, 1])) <= 1e-8 ||
        error("Wavenumber mismatch for albedo=$albedo mode=$mode")
    wn = rrs[:, 1]
    wl = 1e7 ./ wn
    order = sortperm(wl)
    F0 = rrs[:, 8]
    red = (rrs[:, 2] .- nors[:, 2]) ./ F0
    green = rrs[:, 5] ./ F0
    blue = (rrs[:, 2] .+ rrs[:, 5] .- nors[:, 2]) ./ F0
    return (; wavelength = wl[order], red = red[order], green = green[order], blue = blue[order])
end

function add_rgb!(plt, idx::Int, data; title = "", ylabel = "", xlabel = "", legend = false)
    plot!(plt[idx], data.wavelength, data.red; color=:red3, linewidth=1.25,
          label=L"(I_\mathrm{Cab}-I_\mathrm{Rayl})/F_0")
    plot!(plt[idx], data.wavelength, data.green; color=:seagreen, linewidth=1.25,
          label=L"I_\mathrm{RRS}/F_0")
    plot!(plt[idx], data.wavelength, data.blue; color=:royalblue3, linewidth=1.25,
          label=L"(I_\mathrm{Cab}+I_\mathrm{RRS}-I_\mathrm{Rayl})/F_0")
    hline!(plt[idx], [0.0]; color=:gray45, linewidth=0.7, alpha=0.55, label="")
    plot!(plt[idx]; title, ylabel, xlabel, xlims=(761.0, 764.0), grid=true,
          framestyle=:box, legend=(legend ? :bottomright : false),
          tickfontsize=8, guidefontsize=9, titlefontsize=9, legendfontsize=8)
end

function plot_results()
    plt = plot(
        layout=(length(ALBEDOS), 3),
        size=(1800, 950),
        margin=5Plots.mm,
        plot_title="w/Fraunhofer, w/O₂ absorption: diagnostic $(VARPI_CAB_LABEL)",
        plot_titlefontsize=15,
    )
    for (row, albedo) in enumerate(ALBEDOS)
        base = read_rgb("baseline", albedo)
        pert = read_rgb("perturbed", albedo)
        diff = (;
            wavelength = base.wavelength,
            red = pert.red .- base.red,
            green = pert.green .- base.green,
            blue = pert.blue .- base.blue,
        )
        ylabel = @sprintf("ρ=%g\nradiance/F₀", albedo)
        xlabel = row == length(ALBEDOS) ? "wavelength [nm]" : ""
        add_rgb!(plt, (row - 1) * 3 + 1, base;
                 title = row == 1 ? "baseline" : "", ylabel, xlabel,
                 legend = row == length(ALBEDOS))
        add_rgb!(plt, (row - 1) * 3 + 2, pert;
                 title = row == 1 ? VARPI_CAB_LABEL : "", xlabel)
        add_rgb!(plt, (row - 1) * 3 + 3, diff;
                 title = row == 1 ? "perturbed - baseline" : "", xlabel)
    end
    png = joinpath(OUTDIR, "$(CASE_TAG)_rgb.png")
    pdf = joinpath(OUTDIR, "$(CASE_TAG)_rgb.pdf")
    savefig(plt, png)
    savefig(plt, pdf)

    summary = joinpath(OUTDIR, "$(CASE_TAG)_rgb_summary.txt")
    open(summary, "w") do io
        println(io, "case=$CASE_TAG")
        println(io, "varpi_cab_mode=$VARPI_CAB_MODE")
        println(io, "varpi_cab_offset=$VARPI_CAB_OFFSET")
        println(io, "varpi_cab_set=$(VARPI_CAB_SET_METADATA)")
        println(io, "geometry=sza19_vza0_vaz0_psurf1000")
        println(io, "source=w/Fraunhofer")
        println(io, "absorption=O2_ABSCO")
        println(io, "output_window_nm=761:764")
        for albedo in ALBEDOS
            base = read_rgb("baseline", albedo)
            pert = read_rgb("perturbed", albedo)
            for (name, curve) in (
                ("red", pert.red .- base.red),
                ("green", pert.green .- base.green),
                ("blue", pert.blue .- base.blue),
            )
                println(io, @sprintf(
                    "albedo=%g curve=%s diff_min=%.12e diff_max=%.12e diff_mean=%.12e diff_max_abs=%.12e",
                    albedo, name, minimum(curve), maximum(curve), mean(curve), maximum(abs.(curve))))
            end
        end
        println(io, "plot_png=$png")
        println(io, "plot_pdf=$pdf")
    end
    println(png)
    println(pdf)
    println(summary)
end

write_yaml()
params = parameters_from_yaml(YAML_PATH)
solar_T = solar_interpolator(params.float_type)
base_model = model_from_parameters(params)
for albedo in ALBEDOS
    run_case(base_model, params, solar_T, albedo)
end
plot_results()
println(OUTDIR)
