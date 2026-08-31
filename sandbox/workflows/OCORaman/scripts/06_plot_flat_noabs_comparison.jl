#!/usr/bin/env julia

using DelimitedFiles
using Printf
using Plots

const ALBEDOS = [0.0, 0.3, 1.0]
const CASE_TAG = get(ENV, "CASE_TAG", "flat_noabs")
const PLOT_TITLE = get(ENV, "PLOT_TITLE", "Flat Sun, No Absorption, 761-764 nm Output Window")

label_float(x) = replace(@sprintf("%.1f", x), "." => "p")

function read_pair(outdir::AbstractString, version::AbstractString, albedo::Float64)
    tag = label_float(albedo)
    rrs = Float64.(readdlm(joinpath(outdir, "$(version)_$(CASE_TAG)_alb$(tag)_rrs.dat")))
    nors = Float64.(readdlm(joinpath(outdir, "$(version)_$(CASE_TAG)_alb$(tag)_nors.dat")))
    wn = rrs[:, 1]
    if maximum(abs.(wn .- nors[:, 1])) > 1e-8
        error("wavenumber mismatch for $version albedo $albedo")
    end
    f0 = rrs[:, 8]
    red = (rrs[:, 2] .- nors[:, 2]) ./ f0
    green = rrs[:, 5] ./ f0
    blue = (rrs[:, 2] .+ rrs[:, 5] .- nors[:, 2]) ./ f0
    return (; wn, red, green, blue)
end

function read_meta(outdir::AbstractString, version::AbstractString, albedo::Float64=0.0)
    path = joinpath(outdir, "$(version)_$(CASE_TAG)_alb$(label_float(albedo))_metadata.txt")
    meta = Dict{String,String}()
    isfile(path) || return meta
    for line in readlines(path)
        occursin("=", line) || continue
        key, value = split(line, "=", limit=2)
        meta[strip(key)] = strip(value)
    end
    return meta
end

function fmt_meta_value(meta::Dict{String,String}, keys::String...)
    for key in keys
        haskey(meta, key) || continue
        value = meta[key]
        parsed = tryparse(Float64, value)
        return parsed === nothing ? value : @sprintf("%.8f", parsed)
    end
    return "n/a"
end

function column_title(outdir::AbstractString, version::AbstractString)
    meta = read_meta(outdir, version)
    if version == "current"
        return "current\nrho_Cab=" *
               fmt_meta_value(meta, "rho_cab_internal") *
               ", rho_Rayl=" *
               fmt_meta_value(meta, "rho_ray_internal")
    end
    return "June 2024\nrho_Rayl=" *
           fmt_meta_value(meta, "rho_ray_model_greek", "rho_ray_yaml", "rho_ray_internal") *
           ", rho_Cab=" *
           fmt_meta_value(meta, "rho_cab_model_greek", "rho_cab_internal")
end

function add_rgb!(p, subplot_index::Int, data; title="", ylabel="", xlabel="", legend=false)
    plot!(
        p[subplot_index],
        data.wn,
        data.red;
        color=:red3,
        linewidth=1.1,
        label="Icab - IRayl",
    )
    plot!(
        p[subplot_index],
        data.wn,
        data.green;
        color=:seagreen,
        linewidth=1.1,
        label="IRRS",
    )
    plot!(
        p[subplot_index],
        data.wn,
        data.blue;
        color=:royalblue3,
        linewidth=1.1,
        label="Delta I",
    )
    hline!(p[subplot_index], [0.0]; color=:gray45, linewidth=0.6, alpha=0.55, label="")
    plot!(
        p[subplot_index];
        title,
        ylabel,
        xlabel,
        grid=true,
        legend=(legend ? :bottomright : false),
        framestyle=:box,
        tickfontsize=8,
        guidefontsize=9,
        titlefontsize=9,
        legendfontsize=8,
    )
end

function main(args)
    length(args) == 1 || error("usage: julia 06_plot_flat_noabs_comparison.jl OUTDIR")
    outdir = abspath(expanduser(args[1]))

    p = plot(
        layout=(length(ALBEDOS), 3),
        size=(1650, 950),
        margin=5Plots.mm,
        plot_title=PLOT_TITLE,
        plot_titlefontsize=15,
    )

    for (row, albedo) in enumerate(ALBEDOS)
        current = read_pair(outdir, "current", albedo)
        june = read_pair(outdir, "june2024", albedo)
        if maximum(abs.(current.wn .- june.wn)) > 1e-8
            error("current/june2024 grid mismatch for albedo $albedo")
        end
        diff = (;
            wn=current.wn,
            red=current.red .- june.red,
            green=current.green .- june.green,
            blue=current.blue .- june.blue,
        )

        ylabel = @sprintf("rho=%g\nradiance / F0", albedo)
        xlabel = row == length(ALBEDOS) ? "wavenumber [cm^-1]" : ""
        add_rgb!(p, (row - 1) * 3 + 1, current; title=(row == 1 ? column_title(outdir, "current") : ""), ylabel, xlabel, legend=(row == length(ALBEDOS)))
        add_rgb!(p, (row - 1) * 3 + 2, june; title=(row == 1 ? column_title(outdir, "june2024") : ""), xlabel)
        add_rgb!(p, (row - 1) * 3 + 3, diff; title=(row == 1 ? "current - June 2024" : ""), xlabel)
    end

    pdf = joinpath(outdir, "$(CASE_TAG)_761_764_current_vs_june2024_rgb.pdf")
    png = joinpath(outdir, "$(CASE_TAG)_761_764_current_vs_june2024_rgb.png")
    savefig(p, pdf)
    savefig(p, png)

    summary = joinpath(outdir, "$(CASE_TAG)_761_764_current_minus_june2024_summary.txt")
    open(summary, "w") do io
        for albedo in ALBEDOS
            current = read_pair(outdir, "current", albedo)
            june = read_pair(outdir, "june2024", albedo)
            for (name, curve) in (("red", :red), ("green", :green), ("blue", :blue))
                diff = getproperty(current, curve) .- getproperty(june, curve)
                println(
                    io,
                    @sprintf(
                        "albedo=%g curve=%s max_abs_diff=%.12e mean_abs_diff=%.12e min_diff=%.12e max_diff=%.12e",
                        albedo,
                        name,
                        maximum(abs.(diff)),
                        sum(abs.(diff)) / length(diff),
                        minimum(diff),
                        maximum(diff),
                    ),
                )
            end
        end
    end

    println(pdf)
    println(png)
    println(summary)
end

main(ARGS)
