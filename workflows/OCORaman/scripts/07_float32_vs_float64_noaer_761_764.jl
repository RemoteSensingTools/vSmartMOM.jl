#!/usr/bin/env julia

# Precision diagnostic for current vSmartMOM Raman/Rayleigh calculations.
#
# Scene:
#   - black Lambertian surface
#   - no aerosol block
#   - solar Fraunhofer structure + O2 absorption, with q-driven H2O enabled
#   - solve with Raman shoulders, but write/plot only 761--764 nm
#   - VZA=0, vSmartMOM vaz=0, SZA in {0, 30, 60, 70}
#
# Output data columns:
#   wn, Icab, Qcab, IRRS, QRRS, IRayl, QRayl, F0,
#   Itot=(Icab+IRRS), Qtot=(Qcab+QRRS),
#   DeltaI=(Icab+IRRS-IRayl), DeltaQ=(Qcab+QRRS-QRayl)

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using CUDA
using Dates
using DelimitedFiles
using Interpolations
using LaTeXStrings
using Printf
using Statistics
using Plots

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using vSmartMOM.SolarModel

const OUTDIR = get(ENV, "OUTDIR",
    joinpath(homedir(), "Raman_misc", "workflows", "diagnostics",
             "float32_vs_float64_noaer_761_764_" *
             Dates.format(now(), "yyyymmdd_HHMMSS")))
const SOLAR_OUT = get(ENV, "SOLAR_OUT",
    "/home/sanghavi/Raman_misc/workflows/worktrees/uni_vSmartMOM_sanghavi_2025-03-18/src/SolarModel/solar.out")
const O2_ABSCO = get(ENV, "O2_ABSCO",
    "/net/fluo/data1/ABSCO_CS_Database/v5.2_final/o2_v52_v2.jld2")
const ARCH = uppercase(get(ENV, "ARCH", "GPU"))
const GPU_DEVICE = parse(Int, get(ENV, "CUDA_DEVICE", "0"))
const FORCE = get(ENV, "FORCE", "0") == "1"

const SZAS = [0.0, 30.0, 60.0, 70.0]
const OUTPUT_WL_MIN = 761.0
const OUTPUT_WL_MAX = 764.0
const SOLVE_WN_MIN = 1e7 / (764.0 + 10.0)
const SOLVE_WN_MAX = 1e7 / (761.0 - 10.0)
const OUTPUT_WN_MIN = 1e7 / OUTPUT_WL_MAX
const OUTPUT_WN_MAX = 1e7 / OUTPUT_WL_MIN
const WN_STEP = 0.1

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

const Q_PROFILE = [
    2.1e-6, 2.85e-6, 3.45e-6, 4.5e-6, 5.7e-6,
    6.45e-6, 7.35e-6, 9.3e-6, 11.4e-6, 13.95e-6,
    15.3e-6, 16.5e-6, 17.7e-6, 20.1e-6, 20.7e-6,
    23.1e-6, 25.2e-6, 26.4e-6, 27.3e-6, 19.65e-6,
    14.1e-6, 8.7e-6, 5.85e-6, 4.95e-6, 4.95e-6,
    16.95e-6, 195.0e-6, 577.5e-6, 1366.5e-6, 3673.5e-6,
    4350.0e-6, 5100.0e-6, 5994.0e-6, 6789.0e-6, 7029.0e-6,
]

const FTYPES = (Float64, Float32)
const LINE_COLORS = Dict(
    :total => "#0072B2",
    :rayl => "#D55E00",
    :delta => "#009E73",
)

fmt_float(x) = replace(@sprintf("%.1f", Float64(x)), "." => "p")
ft_name(::Type{Float64}) = "Float64"
ft_name(::Type{Float32}) = "Float32"
yaml_array(xs) = "[" * join(string.(xs), ", ") * "]"

function configure_architecture!()
    if ARCH == "GPU"
        CUDA.functional() || error("CUDA is not functional, but ARCH=GPU")
        CUDA.device!(GPU_DEVICE)
        CUDA.allowscalar(false)
    elseif ARCH != "CPU"
        error("ARCH must be GPU or CPU; got $ARCH")
    end
end

function write_yaml(::Type{FT}, sza::Float64) where {FT}
    mkpath(joinpath(OUTDIR, "config"))
    ft = ft_name(FT)
    arch_expr = ARCH == "GPU" ? "Architectures.GPU()" : "Architectures.CPU()"
    yaml_path = joinpath(OUTDIR, "config", "noaer_761_764_sza$(fmt_float(sza))_$(ft).yaml")
    open(yaml_path, "w") do io
        println(io, """
radiative_transfer:
  spec_bands:
    - $(SOLVE_WN_MIN):$(WN_STEP):$(SOLVE_WN_MAX)
  surface:
    - LambertianSurfaceScalar{$(ft)}(0.0)
  polarization_type: Stokes_IQ()
  nstreams: 3
  truncation: auto
  depol: -1
  float_type: $(ft)
  architecture: $(arch_expr)

geometry:
  sza: $(sza)
  vza: [0.0]
  vaz: [0.0]
  obs_alt: [0]

atmospheric_profile:
  T: $(yaml_array(T_PROFILE))
  p: $(yaml_array(P_PROFILE))
  q: $(yaml_array(Q_PROFILE))
  profile_reduction: 5

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
""")
    end
    return yaml_path
end

function solar_interpolator(FT)
    isfile(SOLAR_OUT) || error("SOLAR_OUT not found: $SOLAR_OUT")
    solar = solar_transmission_from_file(SOLAR_OUT)
    @info "Using solar transmission" SOLAR_OUT size=size(solar)
    return LinearInterpolation(FT.(solar[4:end, 1]), FT.(solar[4:end, 2]),
                               extrapolation_bc = Line())
end

function solar_F0(FT, nstokes, ν, solar_T)
    P = FT.(planck_spectrum_wn(FT(5777), collect(ν)) .* FT(2.1629e-05 * π))
    Fvec = FT.(solar_T.(ν)) .* P
    F0 = zeros(FT, nstokes, length(ν))
    F0[1, :] .= Fvec
    return F0
end

function make_rrs(FT, nstokes, ν, model, F0)
    ν_center = FT(0.5) * (first(ν) + last(ν))
    i_ref = argmin(abs.(ν .- ν_center))
    effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
    n2, o2 = InelasticScattering.getRamanAtmoConstants(ν_center, effT)
    rs = InelasticScattering.RRS(
        n2 = n2,
        o2 = o2,
        greek_raman = InelasticScattering.GreekCoefs(
            FT[1], FT[1], FT[1], FT[1], FT[1], FT[1]),
        fscattRayl = FT[1],
        ϖ_Cabannes = FT.(copy(model.ϖ_Cabannes)),
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
    CoreRT.getRamanSSProp!(rs, 1e7 / mean(ν), ν)
    rs.ϖ_Cabannes .= FT.(model.ϖ_Cabannes)
    return rs, effT
end

function make_nors(FT, nstokes, model, F0)
    return InelasticScattering.noRS{FT}(
        fscattRayl = FT[1],
        ϖ_Cabannes = ones(FT, length(model.ϖ_Cabannes)),
        bandSpecLim = UnitRange{Int}[],
        iBand = [1],
        F₀ = copy(F0),
        SIF₀ = zeros(FT, nstokes, size(F0, 2)),
    )
end

function out_path(::Type{FT}, sza::Float64) where {FT}
    return joinpath(OUTDIR, "noaer_761_764_sza$(fmt_float(sza))_$(ft_name(FT)).dat")
end

function metadata_path(::Type{FT}, sza::Float64) where {FT}
    return joinpath(OUTDIR, "noaer_761_764_sza$(fmt_float(sza))_$(ft_name(FT))_metadata.txt")
end

function valid_output(path)
    isfile(path) || return false
    try
        data = readdlm(path)
        return size(data, 2) == 12 && all(isfinite, data)
    catch
        return false
    end
end

function run_case!(::Type{FT}, sza::Float64) where {FT}
    path = out_path(FT, sza)
    if !FORCE && valid_output(path)
        @info "Skipping existing valid output" path
        return
    end

    yaml_path = write_yaml(FT, sza)
    params = parameters_from_yaml(yaml_path)
    @info "Building model" float_type=ft_name(FT) sza yaml_path
    model = model_from_parameters(params)
    nstokes = CoreRT.polarization_type(model).n
    ν = CoreRT.get_spec_bands(model)[1]
    keep = findall(x -> OUTPUT_WN_MIN <= x <= OUTPUT_WN_MAX, ν)
    isempty(keep) && error("No output points inside 761--764 nm")

    solar_T = solar_interpolator(FT)
    F0 = solar_F0(FT, nstokes, ν, solar_T)
    sources = CoreRT.SolarBeam(F₀ = F0) +
              CoreRT.SurfaceSIF(SIF₀ = zeros(FT, nstokes, length(ν)))

    rs_rrs, effT = make_rrs(FT, nstokes, ν, model, F0)
    rs_nors = make_nors(FT, nstokes, model, F0)

    @info "Running RRS" float_type=ft_name(FT) sza nSpec=length(ν) nKeep=length(keep)
    rrs_result = CoreRT.rt_run_test(rs_rrs, model, [1]; sources)
    @info "Running noRS" float_type=ft_name(FT) sza nSpec=length(ν)
    nors_result = CoreRT.rt_run_test(rs_nors, model, [1]; sources)

    R = Array(rrs_result[1])
    ieR = Array(rrs_result[3])
    R0 = Array(nors_result[1])

    Icab = Float64.(R[1, 1, keep])
    Qcab = Float64.(R[1, 2, keep])
    IRRS = Float64.(ieR[1, 1, keep])
    QRRS = Float64.(ieR[1, 2, keep])
    IRayl = Float64.(R0[1, 1, keep])
    QRayl = Float64.(R0[1, 2, keep])
    F0out = Float64.(F0[1, keep])
    Itot = Icab .+ IRRS
    Qtot = Qcab .+ QRRS
    DeltaI = Itot .- IRayl
    DeltaQ = Qtot .- QRayl

    data = hcat(Float64.(ν[keep]), Icab, Qcab, IRRS, QRRS,
                IRayl, QRayl, F0out, Itot, Qtot, DeltaI, DeltaQ)
    any(!isfinite, data) && error("Non-finite output for $(ft_name(FT)) SZA=$sza")
    writedlm(path, data)

    open(metadata_path(FT, sza), "w") do io
        println(io, "created=$(now())")
        println(io, "script=$(abspath(PROGRAM_FILE))")
        println(io, "project=$(Base.active_project())")
        println(io, "git_commit=$(readchomp(`git rev-parse HEAD`))")
        println(io, "float_type=$(ft_name(FT))")
        println(io, "architecture=$ARCH")
        println(io, "cuda_device=$(ARCH == "GPU" ? string(GPU_DEVICE) : "none")")
        println(io, "sza=$sza")
        println(io, "vza=0")
        println(io, "vaz=0")
        println(io, "surface_albedo=0")
        println(io, "aerosol=none")
        println(io, "solar_out=$SOLAR_OUT")
        println(io, "o2_absco=$O2_ABSCO")
        println(io, "q_driven_h2o=true")
        println(io, "yaml=$yaml_path")
        println(io, "solve_wn=$(first(ν)):$(ν[2]-ν[1]):$(last(ν))")
        println(io, "output_wn=$(ν[first(keep)]):$(ν[2]-ν[1]):$(ν[last(keep)])")
        println(io, "n_solve=$(length(ν))")
        println(io, "n_output=$(length(keep))")
        println(io, "n_layers=$(length(model.profile.p_full))")
        println(io, "n_fourier_moments=$(CoreRT.n_fourier_moments_bands(model)[1])")
        println(io, "m_max_order=$(CoreRT.m_max_bands(model)[1])")
        println(io, "effective_temperature=$effT")
        println(io, "tau_rayl_center_column=$(sum(model.τ_rayl[1][argmin(abs.(ν .- mean(ν))), :]))")
        println(io, "tau_abs_center_column=$(sum(model.τ_abs[1][argmin(abs.(ν .- mean(ν))), :]))")
        println(io, "varpi_cab_model=$(model.ϖ_Cabannes[1])")
        println(io, "varpi_cab_rs=$(rs_rrs.ϖ_Cabannes[1])")
        println(io, "rrs_weight_sum_after_rt=$(sum(rs_rrs.ϖ_λ₁λ₀))")
        println(io, "data=$path")
    end
    @info "Wrote" path rows=size(data, 1)
    if ARCH == "GPU"
        GC.gc()
        CUDA.reclaim()
    end
end

function load_case(::Type{FT}, sza::Float64) where {FT}
    data = Float64.(readdlm(out_path(FT, sza)))
    order = sortperm(1e7 ./ data[:, 1])
    return data[order, :]
end

function metric_line(name, a32, a64)
    d = a32 .- a64
    scale = maximum(abs.(a64))
    rel = maximum(abs.(d)) / max(scale, eps(Float64))
    return maximum(abs.(d)), sqrt(mean(d .^ 2)), scale, rel
end

function write_metrics()
    metric_file = joinpath(OUTDIR, "float32_vs_float64_metrics.txt")
    quantities = [
        ("Icab", 2), ("Qcab", 3), ("IRRS", 4), ("QRRS", 5),
        ("IRayl", 6), ("QRayl", 7), ("Itot", 9), ("Qtot", 10),
        ("DeltaI", 11), ("DeltaQ", 12),
    ]
    open(metric_file, "w") do io
        println(io, "Float32 - Float64 precision diagnostics")
        println(io, "outdir=$OUTDIR")
        println(io, "scene=black_surface_no_aerosol_solar_abs_761_764")
        println(io, "columns: sza quantity max_abs_diff rms_diff max_abs_float64 max_abs_rel_diff")
        for sza in SZAS
            f64 = load_case(Float64, sza)
            f32 = load_case(Float32, sza)
            size(f64) == size(f32) || error("Shape mismatch for SZA=$sza: $(size(f64)) vs $(size(f32))")
            max_dwn = maximum(abs.(f64[:, 1] .- f32[:, 1]))
            max_dλ = maximum(abs.(1e7 ./ f64[:, 1] .- 1e7 ./ f32[:, 1]))
            println(io, @sprintf("# sza %.1f grid max_abs_dwn_cm-1 %.8e max_abs_dlambda_nm %.8e",
                sza, max_dwn, max_dλ))
            for (name, col) in quantities
                maxabs, rms, scale, rel = metric_line(name, f32[:, col], f64[:, col])
                println(io, @sprintf("%5.1f %-8s %.8e %.8e %.8e %.8e",
                    sza, name, maxabs, rms, scale, rel))
            end
        end
    end
    return metric_file
end

function padded_lims(values; include_zero=true, frac=0.08, zero_pad=1e-10)
    vals = collect(Float64, values)
    isempty(vals) && return nothing
    lo = minimum(vals)
    hi = maximum(vals)
    if include_zero
        lo = min(lo, 0.0)
        hi = max(hi, 0.0)
    end
    span = hi - lo
    if !isfinite(span) || span == 0.0
        pad = max(maximum(abs.(vals)), zero_pad)
        return (-pad, pad)
    end
    pad = max(span * frac, zero_pad)
    return (lo - pad, hi + pad)
end

function setup_panel!(plt, idx; title="", ylabel="", xlabel="", ylims=nothing)
    plot!(
        plt[idx];
        title,
        ylabel,
        xlabel,
        xlims=(OUTPUT_WL_MIN, OUTPUT_WL_MAX),
        ylims,
        framestyle=:box,
        grid=true,
        gridalpha=0.16,
        tickfontsize=8,
        guidefontsize=9,
        titlefontsize=10,
        legendfontsize=7,
        foreground_color_axis=:black,
        foreground_color_border=:black,
    )
end

function add_precision_overlay!(plt, idx, λ, y64a, y32a, y64b, y32b; labela, labelb)
    plot!(plt[idx], λ, y64a; color=LINE_COLORS[:total], linewidth=1.3, label="$(labela) F64")
    plot!(plt[idx], λ, y32a; color=LINE_COLORS[:total], linewidth=1.0, linestyle=:dash, label="$(labela) F32")
    plot!(plt[idx], λ, y64b; color=LINE_COLORS[:rayl], linewidth=1.3, label="$(labelb) F64")
    plot!(plt[idx], λ, y32b; color=LINE_COLORS[:rayl], linewidth=1.0, linestyle=:dash, label="$(labelb) F32")
end

function add_precision_diff!(plt, idx, λ, y64a, y32a, y64b, y32b; labela, labelb)
    plot!(plt[idx], λ, y32a .- y64a; color=LINE_COLORS[:total], linewidth=1.2, label=labela)
    plot!(plt[idx], λ, y32b .- y64b; color=LINE_COLORS[:rayl], linewidth=1.2, label=labelb)
    plot!(plt[idx], [OUTPUT_WL_MIN, OUTPUT_WL_MAX], [0.0, 0.0]; color=:gray45, linewidth=0.7, label="")
end

function plot_stokes_outputs()
    plt = plot(
        layout=(length(SZAS), 4),
        size=(2100, 1250),
        dpi=250,
        margin=3Plots.mm,
        plot_title="Float32 vs Float64, No Aerosol, Black Surface",
        plot_titlefontsize=16,
        background_color=:white,
        foreground_color=:black,
    )
    for (i, sza) in enumerate(SZAS)
        f64 = load_case(Float64, sza)
        f32 = load_case(Float32, sza)
        λ = 1e7 ./ f64[:, 1]
        row = (i - 1) * 4
        xlabel = i == length(SZAS) ? L"\lambda\;[\mathrm{nm}]" : ""
        setup_panel!(plt, row + 1; title="SZA=$(Int(sza))°  I", ylabel=i == 1 ? "radiance" : "", xlabel,
            ylims=padded_lims(vcat(f64[:, 9], f32[:, 9], f64[:, 6], f32[:, 6])))
        setup_panel!(plt, row + 2; title="F32 - F64  I", xlabel,
            ylims=padded_lims(vcat(f32[:, 9] .- f64[:, 9], f32[:, 6] .- f64[:, 6])))
        setup_panel!(plt, row + 3; title="SZA=$(Int(sza))°  Q", ylabel=i == 1 ? "radiance" : "", xlabel,
            ylims=padded_lims(vcat(f64[:, 10], f32[:, 10], f64[:, 7], f32[:, 7])))
        setup_panel!(plt, row + 4; title="F32 - F64  Q", xlabel,
            ylims=padded_lims(vcat(f32[:, 10] .- f64[:, 10], f32[:, 7] .- f64[:, 7])))

        add_precision_overlay!(plt, row + 1, λ, f64[:, 9], f32[:, 9], f64[:, 6], f32[:, 6];
            labela="Cab+RRS", labelb="Rayl")
        add_precision_diff!(plt, row + 2, λ, f64[:, 9], f32[:, 9], f64[:, 6], f32[:, 6];
            labela="Cab+RRS", labelb="Rayl")
        add_precision_overlay!(plt, row + 3, λ, f64[:, 10], f32[:, 10], f64[:, 7], f32[:, 7];
            labela="Cab+RRS", labelb="Rayl")
        add_precision_diff!(plt, row + 4, λ, f64[:, 10], f32[:, 10], f64[:, 7], f32[:, 7];
            labela="Cab+RRS", labelb="Rayl")

        i == 1 || foreach(j -> plot!(plt[row + j]; legend=false), 1:4)
    end
    pdf = joinpath(OUTDIR, "float32_vs_float64_stokes_outputs_761_764.pdf")
    png = joinpath(OUTDIR, "float32_vs_float64_stokes_outputs_761_764.png")
    savefig(plt, pdf)
    savefig(plt, png)
    return pdf, png
end

function plot_rayleigh_error()
    plt = plot(
        layout=(length(SZAS), 4),
        size=(2100, 1250),
        dpi=250,
        margin=3Plots.mm,
        plot_title="Float32 vs Float64 Rayleigh-Approximation Error",
        plot_titlefontsize=16,
        background_color=:white,
        foreground_color=:black,
    )
    for (i, sza) in enumerate(SZAS)
        f64 = load_case(Float64, sza)
        f32 = load_case(Float32, sza)
        λ = 1e7 ./ f64[:, 1]
        row = (i - 1) * 4
        xlabel = i == length(SZAS) ? L"\lambda\;[\mathrm{nm}]" : ""
        setup_panel!(plt, row + 1; title="SZA=$(Int(sza))°  ΔI", ylabel=i == 1 ? "radiance" : "", xlabel,
            ylims=padded_lims(vcat(f64[:, 11], f32[:, 11])))
        setup_panel!(plt, row + 2; title="F32 - F64  ΔI", xlabel,
            ylims=padded_lims(f32[:, 11] .- f64[:, 11]))
        setup_panel!(plt, row + 3; title="SZA=$(Int(sza))°  ΔQ", ylabel=i == 1 ? "radiance" : "", xlabel,
            ylims=padded_lims(vcat(f64[:, 12], f32[:, 12])))
        setup_panel!(plt, row + 4; title="F32 - F64  ΔQ", xlabel,
            ylims=padded_lims(f32[:, 12] .- f64[:, 12]))

        plot!(plt[row + 1], λ, f64[:, 11]; color=LINE_COLORS[:delta], linewidth=1.35, label="F64")
        plot!(plt[row + 1], λ, f32[:, 11]; color=LINE_COLORS[:rayl], linewidth=1.0, linestyle=:dash, label="F32")
        plot!(plt[row + 2], λ, f32[:, 11] .- f64[:, 11]; color=:black, linewidth=1.2, label="")
        plot!(plt[row + 2], [OUTPUT_WL_MIN, OUTPUT_WL_MAX], [0.0, 0.0]; color=:gray45, linewidth=0.7, label="")
        plot!(plt[row + 3], λ, f64[:, 12]; color=LINE_COLORS[:delta], linewidth=1.35, label="F64")
        plot!(plt[row + 3], λ, f32[:, 12]; color=LINE_COLORS[:rayl], linewidth=1.0, linestyle=:dash, label="F32")
        plot!(plt[row + 4], λ, f32[:, 12] .- f64[:, 12]; color=:black, linewidth=1.2, label="")
        plot!(plt[row + 4], [OUTPUT_WL_MIN, OUTPUT_WL_MAX], [0.0, 0.0]; color=:gray45, linewidth=0.7, label="")

        i == 1 || foreach(j -> plot!(plt[row + j]; legend=false), 1:4)
    end
    pdf = joinpath(OUTDIR, "float32_vs_float64_rayleigh_error_761_764.pdf")
    png = joinpath(OUTDIR, "float32_vs_float64_rayleigh_error_761_764.png")
    savefig(plt, pdf)
    savefig(plt, png)
    return pdf, png
end

function write_manifest(paths)
    manifest = joinpath(OUTDIR, "manifest.txt")
    open(manifest, "w") do io
        println(io, "created=$(now())")
        println(io, "script=$(abspath(PROGRAM_FILE))")
        println(io, "outdir=$OUTDIR")
        println(io, "arch=$ARCH")
        println(io, "gpu_device=$(ARCH == "GPU" ? string(GPU_DEVICE) : "none")")
        println(io, "solar_out=$SOLAR_OUT")
        println(io, "o2_absco=$O2_ABSCO")
        println(io, "sza=$(join(SZAS, ","))")
        println(io, "vza=0")
        println(io, "vaz=0")
        println(io, "surface_albedo=0")
        println(io, "aerosol=none")
        println(io, "solve_wn=$(SOLVE_WN_MIN):$(WN_STEP):$(SOLVE_WN_MAX)")
        println(io, "output_nm=$(OUTPUT_WL_MIN):$(OUTPUT_WL_MAX)")
        for path in paths
            println(io, "output=$path")
        end
    end
    return manifest
end

function main()
    mkpath(OUTDIR)
    configure_architecture!()
    for FT in FTYPES
        for sza in SZAS
            run_case!(FT, sza)
        end
    end
    metrics = write_metrics()
    stokes_pdf, stokes_png = plot_stokes_outputs()
    err_pdf, err_png = plot_rayleigh_error()
    manifest = write_manifest([metrics, stokes_pdf, stokes_png, err_pdf, err_png])
    println(OUTDIR)
    println(metrics)
    println(stokes_pdf)
    println(stokes_png)
    println(err_pdf)
    println(err_png)
    println(manifest)
end

main()
