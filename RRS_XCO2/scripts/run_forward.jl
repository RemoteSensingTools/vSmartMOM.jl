#!/usr/bin/env julia

include(joinpath(@__DIR__, "common.jl"))
using .RRSXCO2Common
using vSmartMOM
using JLD2

params = load_parameters()
model = model_from_parameters(params; external_solar=true)
mkpath(joinpath(ROOT, "output"))

for (iband, name) in enumerate(BAND_NAMES)
    @info "Running forward band" iband name
    toa_result = rt_run_toa(model; i_band=iband,
                            sources=sources_for_band(params, iband))
    ν = collect(params.spec_bands[iband])
    λ_nm = wavelengths_nm(ν)
    toa = Array(toa_result)
    boa = nothing # external-solar SFI intentionally computes TOA only
    jldsave(joinpath(ROOT, "output", "forward_$(name).jld2");
            band=name, wavenumber_cm1=ν, wavelength_nm=λ_nm, toa, boa)
end
