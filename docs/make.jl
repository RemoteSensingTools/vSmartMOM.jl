using Documenter
using DocumenterVitepress
using Literate
using vSmartMOM
using vSmartMOM.Absorption
using vSmartMOM.Scattering
using vSmartMOM.CoreRT
using CairoMakie
using Distributions
using FastGaussQuadrature

const PLOTLY_CDN_URL = "https://cdn.plot.ly/plotly-2.35.2.min.js"
const MIE_PREVIEW_CACHE = Ref{Any}(nothing)

function _json_escape(s::AbstractString)
    return replace(s,
        "\\" => "\\\\",
        "\"" => "\\\"",
        "\n" => "\\n",
        "\r" => "\\r",
        "\t" => "\\t",
    )
end

_json_value(x::Nothing) = "null"
_json_value(x::Bool) = x ? "true" : "false"
_json_value(x::Integer) = string(x)
_json_value(x::AbstractFloat) = isfinite(x) ? string(x) : "null"
_json_value(x::AbstractString) = "\"" * _json_escape(x) * "\""
_json_value(x::Symbol) = _json_value(String(x))
_json_value(x::Tuple) = "[" * join((_json_value(v) for v in x), ",") * "]"
_json_value(x::AbstractVector) = "[" * join((_json_value(v) for v in x), ",") * "]"
_json_value(x::NamedTuple) = "{" * join((_json_value(String(k)) * ":" * _json_value(v) for (k, v) in pairs(x)), ",") * "}"
_json_value(x::AbstractDict) = "{" * join((_json_value(String(k)) * ":" * _json_value(v) for (k, v) in pairs(x)), ",") * "}"

function _asset_destinations(subdir)
    destinations = [
        joinpath(@__DIR__, "build", ".documenter", "assets", subdir),
        joinpath(@__DIR__, "build", ".documenter", "public", "assets", subdir),
    ]

    bases_file = joinpath(@__DIR__, "build", "bases.txt")
    if isfile(bases_file)
        append!(destinations,
                [joinpath(@__DIR__, "build", string(i), "assets", subdir)
                 for i in eachindex(readlines(bases_file))])
    end

    return unique(destinations)
end

function _write_plotly_html(path, data, layout; title)
    config = (;
        responsive = true,
        displaylogo = false,
        modeBarButtonsToRemove = ["lasso2d", "select2d"],
    )

    html = """
    <!doctype html>
    <html lang="en">
    <head>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1">
      <title>$title</title>
      <script src="$PLOTLY_CDN_URL"></script>
      <style>
        html, body, #plot {
          width: 100%;
          height: 100%;
          margin: 0;
          background: #ffffff;
        }
      </style>
    </head>
    <body>
      <div id="plot"></div>
      <script>
        const data = $(_json_value(data));
        const layout = $(_json_value(layout));
        const config = $(_json_value(config));
        Plotly.newPlot("plot", data, layout, config);
      </script>
    </body>
    </html>
    """

    mkpath(dirname(path))
    write(path, html)
    return nothing
end

function _write_plotly_asset(filename, data, layout; title)
    for dest in _asset_destinations("plots")
        _write_plotly_html(joinpath(dest, filename), data, layout; title)
    end
    return nothing
end

function _bottom_legend()
    return (; orientation = "h", x = 0.0, y = -0.24,
              xanchor = "left", yanchor = "top")
end

function _standard_font()
    return (; family = "Inter, system-ui, sans-serif")
end

function _quietly(f)
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            return f()
        end
    end
end

function _quickstart_plot()
    params = read_parameters(joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml"))
    params.architecture = vSmartMOM.Architectures.CPU()
    R, T = _quietly() do
        model = model_from_parameters(params)
        rt_run(model)
    end

    data = [
        (;
            type = "bar",
            x = ["TOA reflectance R", "BOA transmittance T"],
            y = Float64[R[1, 1, 1], T[1, 1, 1]],
            marker = (; color = ["#2563eb", "#16a34a"]),
            hovertemplate = "%{x}<br>Stokes I=%{y:.6f}<extra></extra>",
        ),
    ]
    layout = (;
        title = (; text = "Quickstart RT Output", x = 0.02),
        margin = (; l = 56, r = 18, t = 72, b = 52),
        xaxis = (; title = ""),
        yaxis = (; title = "Stokes I", rangemode = "tozero"),
        paper_bgcolor = "#ffffff",
        plot_bgcolor = "#ffffff",
        font = _standard_font(),
    )
    return data, layout
end

function _core_rt_vza_plot()
    params = read_parameters(joinpath(pkgdir(vSmartMOM),
                                      "test", "test_parameters", "PureRayleighParameters.yaml"))
    params.architecture = vSmartMOM.Architectures.CPU()
    params.max_m = 2
    params.l_trunc = 20
    R, _ = _quietly() do
        model = model_from_parameters(params)
        rt_run(model)
    end

    signed_vza = ifelse.(Float64.(params.vaz) .≈ 180.0,
                         .-Float64.(params.vza),
                         Float64.(params.vza))
    order = sortperm(signed_vza)
    data = [
        (;
            type = "scatter",
            mode = "lines+markers",
            name = "Stokes I",
            x = signed_vza[order],
            y = Float64.(R[order, 1, 1]),
            line = (; color = "#2563eb", width = 3),
            marker = (; size = 7),
            hovertemplate = "VZA=%{x:.0f} deg<br>R_I=%{y:.6f}<extra></extra>",
        ),
        (;
            type = "scatter",
            mode = "lines+markers",
            name = "Stokes Q",
            x = signed_vza[order],
            y = Float64.(R[order, 2, 1]),
            line = (; color = "#dc2626", width = 3),
            marker = (; size = 7),
            hovertemplate = "VZA=%{x:.0f} deg<br>R_Q=%{y:.6f}<extra></extra>",
        ),
    ]
    layout = (;
        title = (; text = "Pure-Rayleigh Reflectance by Viewing Angle", x = 0.02),
        margin = (; l = 60, r = 20, t = 78, b = 92),
        xaxis = (; title = "Signed viewing zenith angle (deg)"),
        yaxis = (; title = "TOA reflectance", zeroline = true),
        legend = _bottom_legend(),
        paper_bgcolor = "#ffffff",
        plot_bgcolor = "#ffffff",
        font = _standard_font(),
    )
    return data, layout
end

function _mie_preview()
    if MIE_PREVIEW_CACHE[] === nothing
        size_distribution = LogNormal(log(0.3), log(2.0))
        aerosol = Aerosol(size_distribution, 1.3, 0.001)
        model = make_mie_model(
            NAI2(),
            aerosol,
            0.55,
            Stokes_IQUV(),
            δBGE(30, 2.0),
            30.0,
            500,
        )
        optics = _quietly() do
            compute_aerosol_optical_properties(model)
        end
        MIE_PREVIEW_CACHE[] = optics
    end
    return MIE_PREVIEW_CACHE[]
end

function _scattering_greek_plot()
    optics = _mie_preview()
    greek = optics.greek_coefs
    coefs = [
        ("alpha", greek.α, "#2563eb"),
        ("beta",  greek.β, "#16a34a"),
        ("gamma", greek.γ, "#7c3aed"),
        ("delta", greek.δ, "#dc2626"),
        ("epsilon", greek.ϵ, "#f59e0b"),
        ("zeta", greek.ζ, "#475569"),
    ]
    n = min(80, length(greek.α))
    ell = collect(0:(n - 1))
    data = [
        (;
            type = "scatter",
            mode = "lines",
            name = name,
            x = ell,
            y = Float64.(coef[1:n]),
            line = (; color = color, width = 2),
            hovertemplate = "l=%{x}<br>$name=%{y:.3e}<extra></extra>",
        )
        for (name, coef, color) in coefs
    ]
    layout = (;
        title = (; text = "Mie Greek Coefficients", x = 0.02),
        margin = (; l = 64, r = 20, t = 78, b = 92),
        xaxis = (; title = "Legendre order l"),
        yaxis = (; title = "Coefficient value", zeroline = true),
        legend = _bottom_legend(),
        paper_bgcolor = "#ffffff",
        plot_bgcolor = "#ffffff",
        font = _standard_font(),
    )
    return data, layout
end

function _scattering_phase_plot()
    optics = _mie_preview()

    μ_quad, _ = gausslegendre(420)
    scattering_matrix = reconstruct_phase(optics.greek_coefs, μ_quad)
    Θ = rad2deg.(acos.(μ_quad))
    order = sortperm(Θ)
    ratio = abs.(scattering_matrix.f₁₂ ./ scattering_matrix.f₁₁)

    data = [
        (;
            type = "scatter",
            mode = "lines",
            name = "f11",
            x = Float64.(Θ[order]),
            y = Float64.(scattering_matrix.f₁₁[order]),
            xaxis = "x",
            yaxis = "y",
            line = (; color = "#2563eb", width = 3),
            hovertemplate = "Θ=%{x:.1f} deg<br>f11=%{y:.3e}<extra></extra>",
        ),
        (;
            type = "scatter",
            mode = "lines",
            name = "|f12/f11|",
            x = Float64.(Θ[order]),
            y = Float64.(ratio[order]),
            xaxis = "x2",
            yaxis = "y2",
            line = (; color = "#7c3aed", width = 3),
            hovertemplate = "Θ=%{x:.1f} deg<br>|f12/f11|=%{y:.3f}<extra></extra>",
        ),
    ]
    layout = (;
        title = (; text = "Mie Phase Function Preview", x = 0.02),
        margin = (; l = 66, r = 18, t = 58, b = 54),
        xaxis = (; domain = [0.0, 1.0], anchor = "y", showticklabels = false),
        yaxis = (; title = "f11", type = "log", domain = [0.56, 1.0]),
        xaxis2 = (; title = "Scattering angle (deg)", domain = [0.0, 1.0], anchor = "y2"),
        yaxis2 = (; title = "|f12/f11|", domain = [0.0, 0.42], range = [0.0, 1.0]),
        paper_bgcolor = "#ffffff",
        plot_bgcolor = "#ffffff",
        font = _standard_font(),
        showlegend = false,
    )
    return data, layout
end

function _absorption_lineshape_plot()
    Δν = collect(-0.22:0.001:0.22)
    σ_d = 0.028
    γ_l = 0.035
    gaussian(x) = exp(-(x / σ_d)^2) / (σ_d * sqrt(pi))
    lorentz(x) = γ_l / (pi * (x^2 + γ_l^2))
    doppler = gaussian.(Δν)
    lorentzian = lorentz.(Δν)
    step_ν = Δν[2] - Δν[1]
    voigt = [sum(gaussian.(Δν) .* lorentz.(x .- Δν)) * step_ν for x in Δν]

    data = [
        (;
            type = "scatter",
            mode = "lines",
            name = "Doppler",
            x = Δν,
            y = doppler ./ maximum(doppler),
            line = (; color = "#2563eb", width = 3),
        ),
        (;
            type = "scatter",
            mode = "lines",
            name = "Lorentz",
            x = Δν,
            y = lorentzian ./ maximum(lorentzian),
            line = (; color = "#dc2626", width = 3),
        ),
        (;
            type = "scatter",
            mode = "lines",
            name = "Voigt convolution",
            x = Δν,
            y = voigt ./ maximum(voigt),
            line = (; color = "#16a34a", width = 3),
        ),
    ]
    layout = (;
        title = (; text = "Normalized Line-Shape Families", x = 0.02),
        margin = (; l = 60, r = 20, t = 78, b = 92),
        xaxis = (; title = "Offset from line center (cm^-1)"),
        yaxis = (; title = "Normalized line shape", rangemode = "tozero"),
        legend = _bottom_legend(),
        paper_bgcolor = "#ffffff",
        plot_bgcolor = "#ffffff",
        font = _standard_font(),
    )
    return data, layout
end

function _absorption_temperature_plot()
    ν = collect(6300.0:0.05:6380.0)
    centers = [6308.5, 6321.0, 6335.5, 6352.0, 6366.0]
    strengths_220 = [0.80, 1.00, 0.72, 0.34, 0.18]
    strengths_290 = [0.34, 0.55, 0.76, 0.92, 0.82]
    width_220 = 0.34
    width_290 = 0.44
    band_220 = zeros(length(ν))
    band_290 = zeros(length(ν))
    for (center, s220, s290) in zip(centers, strengths_220, strengths_290)
        band_220 .+= s220 .* exp.(-0.5 .* ((ν .- center) ./ width_220).^2)
        band_290 .+= s290 .* exp.(-0.5 .* ((ν .- center) ./ width_290).^2)
    end
    band_220 ./= maximum(band_220)
    band_290 ./= maximum(band_290)
    diff = band_220 .- band_290
    data = [
        (;
            type = "scatter",
            mode = "lines",
            name = "220 K",
            x = ν,
            y = band_220,
            line = (; color = "#2563eb", width = 3),
            yaxis = "y",
        ),
        (;
            type = "scatter",
            mode = "lines",
            name = "290 K",
            x = ν,
            y = band_290,
            line = (; color = "#dc2626", width = 3),
            yaxis = "y",
        ),
        (;
            type = "scatter",
            mode = "lines",
            name = "220 K - 290 K",
            x = ν,
            y = diff,
            line = (; color = "#16a34a", width = 3),
            yaxis = "y2",
        ),
    ]
    layout = (;
        title = (; text = "Temperature Redistribution Across a CO2 Band", x = 0.02),
        margin = (; l = 64, r = 20, t = 78, b = 96),
        xaxis = (; title = "Wavenumber (cm^-1)"),
        yaxis = (; title = "Normalized cross section", domain = [0.48, 1.0]),
        yaxis2 = (; title = "Difference", domain = [0.0, 0.34], zeroline = true),
        legend = _bottom_legend(),
        paper_bgcolor = "#ffffff",
        plot_bgcolor = "#ffffff",
        font = _standard_font(),
    )
    return data, layout
end

function _jacobian_aod_plot()
    params = read_parameters(joinpath(pkgdir(vSmartMOM),
                                      "test", "test_parameters", "JacobianTestFast.yaml"))
    params.architecture = vSmartMOM.Architectures.CPU()
    model, lin_model = _quietly() do
        model_from_parameters_lin(params)
    end
    NAer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, _, dR, _ = _quietly() do
        rt_run_lin(model, lin_model, NAer, NGas, NSurf)
    end
    layout_params = CoreRT.ParameterLayout(aerosol_params = 7, n_aerosols = NAer,
                                           n_gases = NGas, n_surface = NSurf)
    aod_idx = first(CoreRT.aerosol_range(layout_params, 1))
    x = collect(1:size(R, 3))
    traces = Any[
        (;
            type = "scatter",
            mode = "lines+markers",
            name = "R_I, VZA 0 deg",
            x = x,
            y = Float64.(R[1, 1, :]),
            line = (; color = "#2563eb", width = 3),
            marker = (; size = 7),
            yaxis = "y",
        ),
        (;
            type = "scatter",
            mode = "lines+markers",
            name = "dR_I/dtau_ref, VZA 0 deg",
            x = x,
            y = Float64.(dR[1, 1, :, aod_idx]),
            line = (; color = "#16a34a", width = 3),
            marker = (; size = 7),
            yaxis = "y2",
        ),
    ]
    if size(R, 1) > 1
        push!(traces,
            (;
                type = "scatter",
                mode = "lines+markers",
                name = "R_I, VZA 30 deg",
                x = x,
                y = Float64.(R[2, 1, :]),
                line = (; color = "#7c3aed", width = 2, dash = "dash"),
                marker = (; size = 7),
                yaxis = "y",
            ),
            (;
                type = "scatter",
                mode = "lines+markers",
                name = "dR_I/dtau_ref, VZA 30 deg",
                x = x,
                y = Float64.(dR[2, 1, :, aod_idx]),
                line = (; color = "#dc2626", width = 2, dash = "dash"),
                marker = (; size = 7),
                yaxis = "y2",
            ),
        )
    end
    layout = (;
        title = (; text = "Aerosol Optical Depth Jacobian", x = 0.02),
        margin = (; l = 72, r = 20, t = 78, b = 104),
        xaxis = (; title = "Spectral grid index"),
        yaxis = (; title = "TOA reflectance", domain = [0.54, 1.0]),
        yaxis2 = (; title = "dR_I / dtau_ref", domain = [0.0, 0.34], zeroline = true),
        legend = _bottom_legend(),
        paper_bgcolor = "#ffffff",
        plot_bgcolor = "#ffffff",
        font = _standard_font(),
    )
    return traces, layout
end

function _surface_brdf_plot()
    vza = collect(-80.0:1.0:80.0)
    lambertian = fill(0.06, length(vza))
    rpv = @. 0.04 + 0.18 / (1 + ((vza + 35) / 16)^2)
    rossli = @. 0.09 + 0.035 * cosd(vza)^2 + 0.025 * exp(-((vza + 45) / 20)^2)
    coxmunk = @. 0.012 + 0.34 * exp(-((vza - 35) / 7)^2)
    data = [
        (; type = "scatter", mode = "lines", name = "Lambertian", x = vza, y = lambertian, line = (; color = "#475569", width = 3)),
        (; type = "scatter", mode = "lines", name = "RPV hotspot", x = vza, y = rpv, line = (; color = "#2563eb", width = 3)),
        (; type = "scatter", mode = "lines", name = "Ross-Li", x = vza, y = rossli, line = (; color = "#16a34a", width = 3)),
        (; type = "scatter", mode = "lines", name = "Cox-Munk glint", x = vza, y = coxmunk, line = (; color = "#f59e0b", width = 3)),
    ]
    layout = (;
        title = (; text = "Surface BRDF Angular Slices", x = 0.02),
        margin = (; l = 64, r = 20, t = 78, b = 96),
        xaxis = (; title = "Signed view zenith angle (deg)"),
        yaxis = (; title = "Relative BRDF / reflectance"),
        legend = _bottom_legend(),
        paper_bgcolor = "#ffffff",
        plot_bgcolor = "#ffffff",
        font = _standard_font(),
    )
    return data, layout
end

function _canopy_spectral_plot()
    λ = collect(450.0:5.0:900.0)
    leaf_reflectance = @. 0.055 + 0.40 / (1 + exp(-(λ - 710) / 18))
    leaf_transmittance = @. 0.025 + 0.35 / (1 + exp(-(λ - 720) / 20))
    bare_soil = @. 0.10 + 0.00012 * (λ - 450)
    canopy_toa = @. 0.035 + 0.31 / (1 + exp(-(λ - 705) / 22)) + 0.015 * exp(-((λ - 550) / 55)^2)
    data = [
        (; type = "scatter", mode = "lines", name = "Leaf reflectance", x = λ, y = leaf_reflectance, line = (; color = "#16a34a", width = 3)),
        (; type = "scatter", mode = "lines", name = "Leaf transmittance", x = λ, y = leaf_transmittance, line = (; color = "#65a30d", width = 3, dash = "dash")),
        (; type = "scatter", mode = "lines", name = "Canopy TOA response", x = λ, y = canopy_toa, line = (; color = "#2563eb", width = 3)),
        (; type = "scatter", mode = "lines", name = "Bare soil", x = λ, y = bare_soil, line = (; color = "#92400e", width = 3)),
    ]
    layout = (;
        title = (; text = "Canopy Red-Edge Response", x = 0.02),
        margin = (; l = 64, r = 20, t = 78, b = 96),
        xaxis = (; title = "Wavelength (nm)"),
        yaxis = (; title = "Reflectance / transmittance", range = [0, 0.55]),
        legend = _bottom_legend(),
        paper_bgcolor = "#ffffff",
        plot_bgcolor = "#ffffff",
        font = _standard_font(),
    )
    return data, layout
end

function generate_tutorial_plot_assets()
    plot_specs = [
        ("quickstart_rt_response.html", _quickstart_plot, "Quickstart RT Output"),
        ("core_rt_vza_response.html", _core_rt_vza_plot, "CoreRT Viewing-Angle Response"),
        ("scattering_greek_coefficients.html", _scattering_greek_plot, "Scattering Greek Coefficients"),
        ("scattering_phase_preview.html", _scattering_phase_plot, "Scattering Phase Preview"),
        ("absorption_lineshape_families.html", _absorption_lineshape_plot, "Absorption Line-Shape Families"),
        ("absorption_temperature_band.html", _absorption_temperature_plot, "Absorption Temperature Band"),
        ("jacobian_aod_response.html", _jacobian_aod_plot, "AOD Jacobian Response"),
        ("surface_brdf_slices.html", _surface_brdf_plot, "Surface BRDF Slices"),
        ("canopy_spectral_response.html", _canopy_spectral_plot, "Canopy Spectral Response"),
    ]

    for (filename, make_plot, title) in plot_specs
        data, layout = make_plot()
        _write_plotly_asset(filename, data, layout; title)
    end
    return nothing
end

function copy_landing_page_icons()
    icon_src = joinpath(@__DIR__, "src", "assets", "icons")
    isdir(icon_src) || return nothing

    destinations = [joinpath(@__DIR__, "build", ".documenter", "public", "assets", "icons")]
    bases_file = joinpath(@__DIR__, "build", "bases.txt")
    if isfile(bases_file)
        bases = readlines(bases_file)
        append!(destinations,
                [joinpath(@__DIR__, "build", string(i), "assets", "icons")
                 for i in eachindex(bases)])
    end

    for dest in destinations
        mkpath(dest)
        for icon in readdir(icon_src)
            (endswith(icon, ".svg") || endswith(icon, ".png")) || continue
            cp(joinpath(icon_src, icon), joinpath(dest, icon); force = true)
        end
    end
    return nothing
end

function build()

    tutorials = ["Tutorial_QuickStart.jl", "Tutorial_Absorption.jl", "Tutorial_Scattering.jl", "Tutorial_MieDeepDive.jl", "Tutorial_IO.jl", "Tutorial_CoreRT.jl", "Tutorial_Surfaces.jl", "Tutorial_Canopy.jl", "Tutorial_Jacobians.jl", "Tutorial_GPU.jl", "Tutorial_HybridAD.jl"]
    tutorials_paths = [joinpath(@__DIR__, "src", "pages", "tutorials", tutorial) for tutorial in tutorials]

    for tutorial in tutorials_paths
        Literate.markdown(tutorial, joinpath(@__DIR__, "src", "pages", "tutorials");
                          codefence = "```julia" => "```")
    end

    tutorials_md = [joinpath("pages", "tutorials", tutorial[1:end-3]) * ".md" for tutorial in tutorials]

    pages = Any[
        "Home"                  => "index.md",
        "User Guide"            => Any[
                                    "Quick Start (5 min)" => "pages/quickstart.md",
                                    "Conventions ⚠ (read first)" => "pages/conventions.md",
                                    "Configure a Scene" => "pages/IO/Overview.md",
                                    "Compute Jacobians" => "pages/jacobians.md",
                                    "Run on GPU" => "pages/gpu.md",
                                   ],
        "RT basics"             => Any[
                                    "1 · Problem & MOM Thesis" => "pages/concepts/01_overview.md",
                                    "2 · Vector RTE & Discretization" => "pages/concepts/02_rt_theory.md",
                                    "3 · Layer Optical Properties" => "pages/concepts/03_layer_optics.md",
                                    "3a · Gas Absorption" => "pages/concepts/03a_absorption.md",
                                    "3b · Mie & Rayleigh" => "pages/concepts/03b_scattering.md",
                                    "3c · Mixing & δ-M Truncation" => "pages/concepts/03c_mixing.md",
                                    "4 · The MOM Solver" => "pages/concepts/04_mom_solver.md",
                                    "5 · Surfaces" => "pages/concepts/05_surfaces.md",
                                    "6 · Linearization" => "pages/concepts/06_linearization.md",
                                    "7 · Architecture-Agnostic Code" => "pages/concepts/07_architecture.md",
                                    "8 · Inelastic Extension" => "pages/concepts/08_inelastic.md",
                                    "Aerosols (microphysics)" => "pages/Aerosols/Overview.md",
                                    "Solar Model (irradiance)" => "pages/SolarModel/Overview.md",
                                   ],
        "Validation"            => Any[
                                    "Benchmarks" => "pages/benchmarks.md",
                                   ],
        "Developer Guides"      => Any[
                                    "Add a Surface BRDF" => "pages/extending/surfaces.md",
                                    "Add a Raman Mode" => "pages/extending/raman.md",
                                    "GEOS-Chem Integration" => "pages/geoschem_integration.md",
                                   ],
        "Tutorials"             => tutorials_md,
        "Library"               => Any[
                                    "Overview" => "pages/api_reference.md",
                                    "Top-Level API" => "pages/api/top_level.md",
                                    "CoreRT" => "pages/api/core_rt.md",
                                    "IO" => "pages/api/io.md",
                                    "Absorption" => "pages/api/absorption.md",
                                    "Scattering" => "pages/api/scattering.md",
                                    "Surfaces" => "pages/api/surfaces.md",
                                    "Inelastic Scattering" => "pages/api/inelastic.md",
                                    "Aerosols" => "pages/api/aerosols.md",
                                    "SolarModel" => "pages/api/solar_model.md",
                                    "Experimental Helpers" => "pages/api/experimental.md",
                                    "Developer Coverage" => "pages/internal_api_coverage.md",
                                    "Function Index" => "pages/api/function_index.md",
                                   ],
        "Resources"             => Any[
                                    "Input Schema" => "pages/IO/Schema.md",
                                    "IO Examples" => "pages/IO/Examples.md",
                                    "HITRAN Data Management" => "pages/Absorption/HITRAN_Data.md",
                                    "References" => "pages/vSmartMOM/References.md",
                                    "Release Notes" => "pages/release_notes.md",
                                   ],
    ]

    ref_name = get(ENV, "GITHUB_REF_NAME", "")
    deploy_devbranch = if ref_name == "unified-vsmartmom"
        "unified-vsmartmom"
    elseif ref_name == "sanghavi-unified"
        "sanghavi-unified"
    else
        "main"
    end
    deploy_devurl = if ref_name == "unified-vsmartmom"
        "unified-vsmartmom"
    elseif ref_name == "sanghavi-unified"
        "sanghavi-unified"
    else
        "dev"
    end

    format = MarkdownVitepress(
        repo = "github.com/RemoteSensingTools/vSmartMOM.jl",
        devbranch = deploy_devbranch,
        devurl = deploy_devurl,
        deploy_decision = get(ENV, "CI", "false") == "true" ?
            nothing : Documenter.DeployDecision(all_ok = false),
        # NOTE: do NOT pass `deploy_url` here — DocumenterVitepress 0.3's
        # url-parser miscomputes the Vitepress `base` when it's set to a full
        # `https://host/path` URL (it strips only `"https:"` and leaves the
        # hostname in the path, producing
        # `base: '/remotesensingtools.github.io/vSmartMOM.jl/...'`). When the
        # arg is omitted, the fallback path `"/" * splitpath(repo)[end]`
        # correctly yields `/vSmartMOM.jl/`.
        description = "Polarized atmospheric radiative transfer and remote sensing tools in Julia.",
        assets = [
            "assets/favicon.ico",
            "assets/logo.png",
            "assets/icons/logo.png",
            "assets/icons/scattering.png",
            "assets/icons/absorption.png",
            "assets/icons/radiative_transfer.png",
        ],
        sidebar_drawer = true,
    )

    # This way it shows warnings of functions that have not been documented
    makedocs(
        sitename = "vSmartMOM",
        format = format,
        clean = true,
        modules = [vSmartMOM],
        pages = pages,
        checkdocs = :exports,
    )

    generate_tutorial_plot_assets()
    copy_landing_page_icons()

end

build()

# Deploy only in CI contexts; local docs builds should not attempt git deployment.
if get(ENV, "CI", "false") == "true"
    # Keep main as the canonical dev docs, but allow branch-specific dev deployment
    # for unified-vsmartmom so pushes to that branch publish docs as well.
    ref_name = get(ENV, "GITHUB_REF_NAME", "")
    deploy_devbranch = if ref_name == "unified-vsmartmom"
        "unified-vsmartmom"
    elseif ref_name == "sanghavi-unified"
        "sanghavi-unified"
    else
        "main"
    end
    deploy_devurl = if ref_name == "unified-vsmartmom"
        "unified-vsmartmom"
    elseif ref_name == "sanghavi-unified"
        "sanghavi-unified"
    else
        "dev"
    end

    DocumenterVitepress.deploydocs(
        root = @__DIR__,
        repo = "github.com/RemoteSensingTools/vSmartMOM.jl.git",
        target = "build",
        push_preview = true,
        devbranch = deploy_devbranch,
        devurl = deploy_devurl,
    )
end
