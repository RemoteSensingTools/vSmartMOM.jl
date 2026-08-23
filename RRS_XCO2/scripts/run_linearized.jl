#!/usr/bin/env julia

include(joinpath(@__DIR__, "common.jl"))
using .RRSXCO2Common
using vSmartMOM
using vSmartMOM.CoreRT
using JLD2
using DelimitedFiles

function parameter_labels(layout, band)
    result = ["surface_pressure"]
    aerosols = ("sulfate", "organic_carbon", "stratospheric_sulfate")
    fields = ("tau_ref", "n_real", "n_imag", "median_radius",
              "geometric_width", "profile_location", "profile_width")
    for aerosol in aerosols, field in fields
        push!(result, "$(aerosol)_$(field)")
    end
    n_layers = layout.n_gases ÷ 2
    for species in ("h2o", "co2"), iz in 1:n_layers
        push!(result, "$(species)_layer_$(lpad(iz, 2, '0'))")
    end
    append!(result, ("surface_P0", "surface_P1", "surface_P2"))
    layout.n_sif == 2 && append!(result, ("SIF760", "mSIF"))
    length(result) == CoreRT.n_total(layout) || error("Jacobian label/layout mismatch")
    return result
end

params = load_parameters()
NAer = length(params.scattering_params.rt_aerosols)
NAer == 3 || error("This experiment requires three aerosol species")
# The analytic-linearized SFI path uses the legacy embedded-solar
# representation: μ₀ is retained as a zero-weight angular node. This makes
# the polarized operator 21×21 for six diffuse streams and IQU polarization.
model, lin_model = model_from_parameters(LinMode(), params; external_solar=false)
NSurf = 3 # P0, P1, P2 for LambertianSurfaceLegendre
mkpath(joinpath(ROOT, "output"))

for (iband, name) in enumerate(BAND_NAMES)
    NGas = size(lin_model.τ̇_abs[iband], 1)
    @info "Running linearized band" iband name NAer NGas NSurf
    result = rt_run(model, lin_model, NAer, NGas, NSurf; i_band=iband,
                    sources=sources_for_band(params, iband))
    ν = collect(params.spec_bands[iband])
    λ_nm = wavelengths_nm(ν)
    toa = Array(result.toa)
    jacobian = Array(result.toa_jacobian)
    layout = result.layout
    labels = parameter_labels(layout, name)
    jldsave(joinpath(ROOT, "output", "linearized_$(name).jld2");
            band=name, wavenumber_cm1=ν, wavelength_nm=λ_nm, toa, jacobian,
            parameter_labels=labels,
            n_atmosphere=layout.n_atmosphere,
            aerosol_params=layout.aerosol_params,
            n_aerosols=layout.n_aerosols, n_gases=layout.n_gases,
            n_surface=layout.n_surface)
    open(joinpath(ROOT, "output", "linearized_$(name).csv"), "w") do io
        println(io, join(vcat("wavelength_nm", labels), ','))
        writedlm(io, hcat(λ_nm, Array(result.toa_jacobian[1, 1, :, :])), ',')
    end
end
