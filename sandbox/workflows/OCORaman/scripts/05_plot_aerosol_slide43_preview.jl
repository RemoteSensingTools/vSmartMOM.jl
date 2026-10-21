#!/usr/bin/env julia

# Slide-43-style figure for the aerosol Raman/Rayleigh ASCII outputs.
#
# Expected input columns:
#   rrs:  wn, Icab, Qcab, Ucab, IRRS, QRRS, URRS, F0
#   nors: wn, IRayl, QRayl, URayl, F0
#
# The script is intentionally tolerant of missing cases. Re-run it while the
# expensive GPU job progresses and the figure will fill in case by case.

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using DelimitedFiles
using LaTeXStrings
using Printf
using Plots
using Statistics

const OUTDIR = get(ENV, "OUTDIR", "/home/sanghavi/data/RamanSIFgrid/aerosol_test_new")
const BAND_TAG = get(ENV, "BAND_TAG", "ABO2")
const PLOTDIR = get(ENV, "PLOTDIR", joinpath(OUTDIR, "plots"))

const ALBEDOS = [1.0, 0.5, 0.1, 0.0]
const AODS = [0.0, 0.1, 0.5]
# ABO2 output includes Raman shoulders around the 761--764 nm analysis window.
const LAMBDA_LIMS = (755.0, 775.0)

const AOD_COLORS = Dict(
    0.0 => "#0072B2", # Okabe-Ito blue
    0.1 => "#009E73", # Okabe-Ito bluish green
    0.5 => "#D55E00", # Okabe-Ito vermillion
)

const ALBEDO_ALPHAS = Dict(
    1.0 => 0.98,
    0.5 => 0.74,
    0.1 => 0.50,
    0.0 => 0.34,
)

const ROW_LABELS = Dict(
    1.0 => L"\rho = 1",
    0.5 => L"\rho = 0.5",
    0.1 => L"\rho = 0.1",
    0.0 => L"\rho = 0",
)

const LAMBDA_LABEL = L"\lambda\;[\mathrm{nm}]"
const RADIANCE_LABEL = L"\mathrm{Radiance}\;[\mathrm{mW}\,\mathrm{m}^{-2}\,\mathrm{sr}^{-1}\,\mathrm{nm}^{-1}]"
const DELTA_LABEL = L"\Delta I\;[\mathrm{mW}\,\mathrm{m}^{-2}\,\mathrm{sr}^{-1}\,\mathrm{nm}^{-1}]"

label_float(x) = replace(@sprintf("%.1f", Float64(x)), "." => "p")
label_value(x) = iszero(x - round(x)) ? string(Int(round(x))) : @sprintf("%.1f", Float64(x))
row_label_latex(albedo) = latexstring("\\rho = $(label_value(albedo))")
delta_title(albedo) = latexstring("\\Delta I = I_{\\mathrm{Cab}} + I_{\\mathrm{RRS}} - I_{\\mathrm{Rayl}}\\quad(\\rho = $(label_value(albedo)))")
delta_elastic_title(albedo) = latexstring("\\Delta I_{\\mathrm{el}} = I_{\\mathrm{Cab}} - I_{\\mathrm{Rayl}}\\quad(\\rho = $(label_value(albedo)))")

function case_paths(aod::Float64, albedo::Float64)
    alb_str = label_float(albedo)
    if iszero(aod)
        stem = "aer_AOD0p0_sza0_vza0p0_vaz0p0_alb$(alb_str)_psurf1000hpa"
    else
        stem = "aer_AOD$(label_float(aod))_z2_r0p1_sza0_vza0p0_vaz0p0_alb$(alb_str)_psurf1000hpa"
    end
    return (
        joinpath(OUTDIR, "$(stem)_rrs_$(BAND_TAG).dat"),
        joinpath(OUTDIR, "$(stem)_nors_$(BAND_TAG).dat"),
    )
end

function load_case(aod::Float64, albedo::Float64)
    rrs_path, nors_path = case_paths(aod, albedo)
    (isfile(rrs_path) && isfile(nors_path)) || return nothing

    rrs = Float64.(readdlm(rrs_path))
    nors = Float64.(readdlm(nors_path))
    size(rrs, 2) >= 8 || error("Expected at least 8 columns in $rrs_path")
    size(nors, 2) >= 5 || error("Expected at least 5 columns in $nors_path")
    size(rrs, 1) == size(nors, 1) || error("Row mismatch between $rrs_path and $nors_path")

    wn = rrs[:, 1]
    if maximum(abs.(wn .- nors[:, 1])) > 1e-7
        error("Wavenumber mismatch between $rrs_path and $nors_path")
    end

    lambda = 1e7 ./ wn
    jac = 1e7 ./ lambda.^2
    order = sortperm(lambda)
    icab = rrs[:, 2]
    irrs = rrs[:, 5]
    irayl = nors[:, 2]

    return (;
        lambda=lambda[order],
        # vSmartMOM ASCII output is spectral radiance per cm^-1. These panels
        # are drawn against wavelength, so convert to per nm with dν/dλ.
        delta=((icab .+ irrs .- irayl) .* jac)[order],
        delta_elastic=((icab .- irayl) .* jac)[order],
        irrs=(irrs .* jac)[order],
        rrs_path,
        nors_path,
    )
end

function quantity_values(data, quantity::Symbol)
    if quantity === :delta
        return data.delta
    elseif quantity === :delta_elastic
        return data.delta_elastic
    elseif quantity === :irrs
        return data.irrs
    end
    error("Unknown quantity: $quantity")
end

function padded_lims(values; include_zero=true, frac=0.08)
    vals = collect(skipmissing(values))
    isempty(vals) && return nothing
    lo = minimum(vals)
    hi = maximum(vals)
    if include_zero
        lo = min(lo, 0.0)
        hi = max(hi, 0.0)
    end
    span = hi - lo
    if !(isfinite(span)) || span == 0
        pad = max(abs(hi), 1.0) * 0.05
    else
        pad = span * frac
    end
    return (lo - pad, hi + pad)
end

function ylims_for(albedo::Float64, quantity::Symbol)
    y = Float64[]
    for aod in AODS
        data = load_case(aod, albedo)
        data === nothing && continue
        append!(y, quantity_values(data, quantity))
    end
    return padded_lims(y)
end

function ylims_for_irrs()
    y = Float64[]
    for albedo in ALBEDOS, aod in AODS
        data = load_case(aod, albedo)
        data === nothing && continue
        append!(y, data.irrs)
    end
    return padded_lims(y)
end

function setup_panel!(plt, idx::Int; title="", ylabel="", xlabel="", ylims=nothing)
    left_pad = idx in (2, 4, 6, 8) ? 8Plots.mm : 3Plots.mm
    plot!(
        plt[idx],
        collect(LAMBDA_LIMS),
        [0.0, 0.0];
        color=:gray45,
        linewidth=0.7,
        alpha=0.65,
        label="",
    )
    plot!(
        plt[idx];
        title,
        ylabel,
        xlabel,
        ylims,
        xlims=LAMBDA_LIMS,
        framestyle=:box,
        grid=true,
        gridalpha=0.16,
        minorgrid=false,
        tickfontsize=11,
        guidefontsize=13,
        titlefontsize=12,
        legendfontsize=10,
        left_margin=left_pad,
        right_margin=2Plots.mm,
        top_margin=1Plots.mm,
        bottom_margin=7Plots.mm,
        foreground_color_axis=:black,
        foreground_color_border=:black,
    )
end

function setup_shared_delta_label!(plt, idx::Int)
    plot!(
        plt[idx],
        [0.0, 1.0],
        [0.0, 1.0];
        color=:white,
        alpha=0.0,
        label="",
    )
    plot!(
        plt[idx];
        framestyle=:none,
        grid=false,
        xticks=false,
        yticks=false,
        xlims=(0, 1),
        ylims=(0, 1),
        guidefontsize=13,
        left_margin=0Plots.mm,
        right_margin=0Plots.mm,
        top_margin=0Plots.mm,
        bottom_margin=0Plots.mm,
    )
    annotate!(
        plt[idx],
        0.5,
        0.5,
        text("ΔI [mW m⁻² sr⁻¹ nm⁻¹]", 13, :black, rotation=90),
    )
end

function mark_pending!(plt, idx::Int)
    annotate!(
        plt[idx],
        mean(LAMBDA_LIMS),
        0.0,
        text("pending", 12, :gray40, :center),
    )
end

function add_left_panel!(plt, idx::Int, albedo::Float64, quantity::Symbol; title="", ylabel="", xlabel="")
    setup_panel!(plt, idx; title, ylabel, xlabel, ylims=ylims_for(albedo, quantity))

    any_loaded = false
    for aod in AODS
        data = load_case(aod, albedo)
        data === nothing && continue
        y = quantity_values(data, quantity)
        plot!(
            plt[idx],
            data.lambda,
            y;
            color=AOD_COLORS[aod],
            linewidth=1.65,
            alpha=0.96,
            label=L"\tau_\mathrm{aer} = %$(label_value(aod))",
        )
        any_loaded = true
    end
    any_loaded || mark_pending!(plt, idx)
    plot!(plt[idx]; legend=false)
end

function add_irrs_panel!(plt, idx::Int)
    plot!(
        plt[idx];
        title=L"I_\mathrm{RRS}",
        xlabel=LAMBDA_LABEL,
        ylabel=RADIANCE_LABEL,
        xlims=LAMBDA_LIMS,
        ylims=ylims_for_irrs(),
        framestyle=:box,
        grid=true,
        gridalpha=0.16,
        tickfontsize=11,
        guidefontsize=13,
        titlefontsize=12,
        legend=:top,
        legendfontsize=10,
        background_color_legend=:white,
        foreground_color_legend=:gray35,
        left_margin=9Plots.mm,
        right_margin=3Plots.mm,
        top_margin=1Plots.mm,
        bottom_margin=7Plots.mm,
        foreground_color_axis=:black,
        foreground_color_border=:black,
    )

    for aod in AODS
        plot!(
            plt[idx],
            [NaN],
            [NaN];
            color=AOD_COLORS[aod],
            linestyle=:solid,
            linewidth=2.2,
            label=L"\tau_\mathrm{aer} = %$(label_value(aod))",
        )
    end
    plot!(plt[idx], [NaN], [NaN]; color=:white, linewidth=0.0, label=" ")
    for albedo in ALBEDOS
        plot!(
            plt[idx],
            [NaN],
            [NaN];
            color=:gray20,
            linestyle=:solid,
            linewidth=2.0,
            alpha=ALBEDO_ALPHAS[albedo],
            label=ROW_LABELS[albedo],
        )
    end

    any_loaded = false
    for albedo in ALBEDOS
        for aod in AODS
            data = load_case(aod, albedo)
            data === nothing && continue
            plot!(
                plt[idx],
                data.lambda,
                data.irrs;
                color=AOD_COLORS[aod],
                linestyle=:solid,
                linewidth=1.45,
                alpha=ALBEDO_ALPHAS[albedo],
                label="",
            )
            any_loaded = true
        end
    end

    any_loaded || mark_pending!(plt, idx)
end

function write_case_summary(path::AbstractString)
    open(path, "w") do io
        println(io, "Slide-43-style aerosol preview")
        println(io, "outdir = $OUTDIR")
        println(io, "band_tag = $BAND_TAG")
        println(io, "conversion = plotted quantities are file quantities * (1e7 / wavelength_nm^2)")
        println(io, "units = file quantities are mW m^-2 sr^-1 cm; plotted quantities are mW m^-2 sr^-1 nm^-1")
        println(io)
        for albedo in ALBEDOS
            for aod in AODS
                rrs_path, nors_path = case_paths(aod, albedo)
                status = (isfile(rrs_path) && isfile(nors_path)) ? "available" : "pending"
                println(io, @sprintf("albedo=%s AOD=%s %s", label_value(albedo), label_value(aod), status))
            end
        end
    end
end

function main()
    mkpath(PLOTDIR)

    layout = @layout [a{0.04w} grid(4, 2) c{0.32w}]
    plt = plot(
        layout=layout,
        size=(2400, 1330),
        dpi=300,
        margin=1Plots.mm,
        left_margin=2Plots.mm,
        right_margin=3Plots.mm,
        bottom_margin=8Plots.mm,
        top_margin=3Plots.mm,
        plot_title="Rayleigh-Approximation Error in the O₂ A Band",
        plot_titlefontsize=19,
        fontfamily="DejaVu Sans",
        background_color=:white,
        foreground_color=:black,
    )

    setup_shared_delta_label!(plt, 1)

    for (row, albedo) in enumerate(ALBEDOS)
        left = (row - 1) * 2 + 2
        right = left + 1
        row_label = row_label_latex(albedo)
        xlabel = row == length(ALBEDOS) ? LAMBDA_LABEL : ""

        add_left_panel!(
            plt,
            left,
            albedo,
            :delta;
            title=(row == 1 ? delta_title(albedo) : row_label),
            xlabel,
        )
        add_left_panel!(
            plt,
            right,
            albedo,
            :delta_elastic;
            title=(row == 1 ? delta_elastic_title(albedo) : row_label),
            xlabel,
        )
    end

    add_irrs_panel!(plt, 10)

    pdf = joinpath(PLOTDIR, "aerosol_slide43_preview_$(BAND_TAG).pdf")
    png = joinpath(PLOTDIR, "aerosol_slide43_preview_$(BAND_TAG).png")
    pub_pdf = joinpath(PLOTDIR, "aerosol_slide43_publication_$(BAND_TAG).pdf")
    pub_png = joinpath(PLOTDIR, "aerosol_slide43_publication_$(BAND_TAG).png")
    summary = joinpath(PLOTDIR, "aerosol_slide43_preview_$(BAND_TAG)_cases.txt")
    savefig(plt, pdf)
    savefig(plt, png)
    savefig(plt, pub_pdf)
    savefig(plt, pub_png)
    write_case_summary(summary)

    println(pdf)
    println(png)
    println(pub_pdf)
    println(pub_png)
    println(summary)
end

main()
