module RRSXCO2Common

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.SolarModel: planck_spectrum_wn

export ROOT, CONFIG, BAND_NAMES, configure_luts!, load_parameters, wavelengths_nm,
       sources_for_band

const ROOT = normpath(joinpath(@__DIR__, ".."))
const CONFIG = joinpath(ROOT, "config", "oco_grass_3aerosol.yaml")
const BAND_NAMES = ("o2a", "weak_co2", "strong_co2")

"Set portable LUT defaults; explicit environment values always win."
function configure_luts!()
    lut_dir = get(ENV, "VSMARTMOM_HITRAN_LUT_DIR",
                  joinpath(homedir(), "data", "HITRAN_LUTs"))
    defaults = Dict("O2_LUT" => "O2.jld2",
                    "H2O_LUT" => "H2O.jld2",
                    "CO2_LUT" => "CO2.jld2")
    for (key, filename) in defaults
        haskey(ENV, key) || (ENV[key] = joinpath(lut_dir, filename))
        isfile(ENV[key]) || error("Missing $key at $(ENV[key]). Set $key or VSMARTMOM_HITRAN_LUT_DIR.")
    end
    return nothing
end

function load_parameters()
    configure_luts!()
    params = parameters_from_yaml(CONFIG)
    length(params.scattering_params.rt_aerosols) == 3 ||
        error("Expected exactly three aerosol species")
    length(params.spec_bands) == 3 || error("Expected exactly three OCO bands")
    return params
end

wavelengths_nm(ν) = 1.0e7 ./ collect(ν)

"Physical solar beam for every band; retrievable SIF is present only in O2 A."
function sources_for_band(params, iband; SIF760=nothing, mSIF=nothing)
    FT = params.float_type
    ν = FT.(params.spec_bands[iband])
    nPol, nSpec = params.polarization_type.n, length(ν)
    solar = planck_spectrum_wn(FT(5777), collect(ν)) .* FT(2.1629e-5 * π)
    F₀ = zeros(FT, nPol, nSpec)
    @views F₀[1, :] .= solar
    beam = CoreRT.SolarBeam(F₀=F₀)
    iband == 1 || return beam
    state = sif_reference_state(total_sif=0.5, reference_wavelength_nm=760)
    SIF760 === nothing && (SIF760 = state.SIF760)
    mSIF === nothing && (mSIF = state.mSIF)
    return beam + CoreRT.SurfaceSIF(SIF760=FT(SIF760), mSIF=FT(mSIF),
        wavenumber_cm1=ν)
end

end
