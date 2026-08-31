#!/usr/bin/env julia

# Difference plot for the completed aerosol Raman/Rayleigh ASCII outputs.
#
# For each albedo, this plots:
#   (AOD=0.5 quantity) - (AOD=0.1 quantity)
#
# Expected input columns:
#   rrs:  wn, Icab, Qcab, Ucab, IRRS, QRRS, URRS, F0
#   nors: wn, IRayl, QRayl, URayl, F0

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using DelimitedFiles
using Printf
using Plots

const OUTDIR = get(ENV, "OUTDIR", "/home/sanghavi/data/RamanSIFgrid/aerosol_test_new")
const BAND_TAG = get(ENV, "BAND_TAG", "ABO2")
const PLOTDIR = get(ENV, "PLOTDIR", joinpath(OUTDIR, "plots"))

const ALBEDOS = [1.0, 0.5, 0.1, 0.0]
const AOD_HIGH = 0.5
const AOD_LOW = 0.1
const LAMBDA_LIMS = (755.0, 775.0)

const ALBEDO_COLORS = Dict(
    1.0 => :red3,
    0.5 => :darkorange3,
    0.1 => :seagreen4,
    0.0 => :royalblue3,
)

const ALBEDO_STYLES = Dict(
    1.0 => :solid,
    0.5 => :dash,
    0.1 => :dot,
    0.0 => :dashdot,
)

label_float(x) = replace(@sprintf("%.1f", Float64(x)), "." => "p")
label_value(x) = iszero(x - round(x)) ? string(Int(round(x))) : @sprintf("%.1f", Float64(x))

function case_paths(aod::Float64, albedo::Float64)
    alb_str = label_float(albedo)
    stem = "aer_AOD$(label_float(aod))_z2_r0p1_sza0_vza0p0_vaz0p0_alb$(alb_str)_psurf1000hpa"
    return (
        joinpath(OUTDIR, "$(stem)_rrs_$(BAND_TAG).dat"),
        joinpath(OUTDIR, "$(stem)_nors_$(BAND_TAG).dat"),
    )
end

function load_case(aod::Float64, albedo::Float64)
    rrs_path, nors_path = case_paths(aod, albedo)
    (isfile(rrs_path) && isfile(nors_path)) ||
        error("Missing case files for AOD=$(aod), albedo=$(albedo)")

    rrs = Float64.(readdlm(rrs_path))
    nors = Float64.(readdlm(nors_path))
    size(rrs, 2) >= 8 || error("Expected at least 8 columns in $rrs_path")
    size(nors, 2) >= 5 || error("Expected at least 5 columns in $nors_path")
    size(rrs, 1) == size(nors, 1) || error("Row mismatch between $rrs_path and $nors_path")

    wn = rrs[:, 1]
    maximum(abs.(wn .- nors[:, 1])) <= 1e-7 ||
        error("Wavenumber mismatch between $rrs_path and $nors_path")

    lambda = 1e7 ./ wn
    jac = 1e7 ./ lambda.^2
    order = sortperm(lambda)
    icab = rrs[:, 2]
    irrs = rrs[:, 5]
    irayl = nors[:, 2]

    return (;
        lambda=lambda[order],
        # vSmartMOM ASCII output is per cm^-1. Convert to per nm before
        # plotting against wavelength.
        delta=((icab .+ irrs .- irayl) .* jac)[order],
        delta_elastic=((icab .- irayl) .* jac)[order],
        irrs=(irrs .* jac)[order],
    )
end

function diff_case(albedo::Float64)
    high = load_case(AOD_HIGH, albedo)
    low = load_case(AOD_LOW, albedo)
    maximum(abs.(high.lambda .- low.lambda)) <= 1e-8 ||
        error("Wavelength grid mismatch for albedo=$(albedo)")

    return (;
        lambda=high.lambda,
        delta=high.delta .- low.delta,
        delta_elastic=high.delta_elastic .- low.delta_elastic,
        irrs=high.irrs .- low.irrs,
    )
end

function setup_panel!(plt, idx::Int; title="", ylabel="", xlabel="")
    hline!(plt[idx], [0.0]; color=:gray45, linewidth=0.7, alpha=0.65, label="")
    plot!(
        plt[idx];
        title,
        ylabel,
        xlabel,
        xlims=LAMBDA_LIMS,
        framestyle=:box,
        grid=true,
        tickfontsize=8,
        guidefontsize=11,
        titlefontsize=10,
        legend=false,
    )
end

function add_left_panel!(plt, idx::Int, albedo::Float64, quantity::Symbol; title="", ylabel="", xlabel="")
    setup_panel!(plt, idx; title, ylabel, xlabel)
    data = diff_case(albedo)
    y = quantity === :delta ? data.delta : data.delta_elastic
    plot!(
        plt[idx],
        data.lambda,
        y;
        color=ALBEDO_COLORS[albedo],
        linewidth=1.35,
        label="",
    )
end

function add_irrs_panel!(plt, idx::Int)
    setup_panel!(
        plt,
        idx;
        title="ΔI_RRS = I_RRS(AOD 0.5) - I_RRS(AOD 0.1)",
        ylabel="",
        xlabel="wavelength [nm]",
    )
    for albedo in ALBEDOS
        data = diff_case(albedo)
        plot!(
            plt[idx],
            data.lambda,
            data.irrs;
            color=ALBEDO_COLORS[albedo],
            linestyle=ALBEDO_STYLES[albedo],
            linewidth=1.35,
            label="rho=$(label_value(albedo))",
        )
    end
    plot!(plt[idx]; legend=:topright, legendfontsize=8)
end

function write_summary(path::AbstractString)
    open(path, "w") do io
        println(io, "Aerosol AOD-difference preview")
        println(io, "outdir = $OUTDIR")
        println(io, "band_tag = $BAND_TAG")
        println(io, "difference = AOD $(label_value(AOD_HIGH)) - AOD $(label_value(AOD_LOW))")
        println(io, "conversion = plotted quantities are file quantities * (1e7 / wavelength_nm^2)")
        println(io, "units = file quantities are mW m^-2 sr^-1 cm; plotted quantities are mW m^-2 sr^-1 nm^-1")
        println(io)
        for albedo in ALBEDOS
            data = diff_case(albedo)
            for (name, curve) in (
                ("DeltaI", data.delta),
                ("DeltaIe", data.delta_elastic),
                ("IRRS", data.irrs),
            )
                println(
                    io,
                    @sprintf(
                        "albedo=%s quantity=%s min=%.12e max=%.12e mean=%.12e",
                        label_value(albedo),
                        name,
                        minimum(curve),
                        maximum(curve),
                        sum(curve) / length(curve),
                    ),
                )
            end
        end
    end
end

function main()
    mkpath(PLOTDIR)
    layout = @layout [grid(4, 2) c{0.34w}]
    plt = plot(
        layout=layout,
        size=(1650, 980),
        margin=8Plots.mm,
        left_margin=12Plots.mm,
        right_margin=9Plots.mm,
        plot_title="Aerosol sensitivity: AOD 0.5 - AOD 0.1",
        plot_titlefontsize=15,
    )

    for (row, albedo) in enumerate(ALBEDOS)
        left = (row - 1) * 2 + 1
        right = left + 1
        row_label = "rho=$(label_value(albedo))"
        ylabel = ""
        xlabel = row == length(ALBEDOS) ? "wavelength [nm]" : ""

        add_left_panel!(
            plt,
            left,
            albedo,
            :delta;
            title=(row == 1 ? "Δ(I_Cab + I_RRS - I_Rayl), $row_label" : row_label),
            ylabel,
            xlabel,
        )
        add_left_panel!(
            plt,
            right,
            albedo,
            :delta_elastic;
            title=(row == 1 ? "Δ(I_Cab - I_Rayl), $row_label" : row_label),
            xlabel,
        )
    end
    add_irrs_panel!(plt, 9)

    pdf = joinpath(PLOTDIR, "aerosol_slide43_AOD0p5_minus_AOD0p1_$(BAND_TAG).pdf")
    png = joinpath(PLOTDIR, "aerosol_slide43_AOD0p5_minus_AOD0p1_$(BAND_TAG).png")
    summary = joinpath(PLOTDIR, "aerosol_slide43_AOD0p5_minus_AOD0p1_$(BAND_TAG)_summary.txt")
    savefig(plt, pdf)
    savefig(plt, png)
    write_summary(summary)

    println(pdf)
    println(png)
    println(summary)
end

main()
