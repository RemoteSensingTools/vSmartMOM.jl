using OCORamanWorkflow

function parse_float_list(value::AbstractString)
    isempty(strip(value)) && return Float64[]
    return parse.(Float64, split(value, ","))
end

wavelengths = parse_float_list(get(ENV, "WAVELENGTH_NM_LIST", "745,760,765,770,785"))
temperature = parse(Float64, get(ENV, "TEMPERATURE_K", "300"))
outdir = joinpath(output_root(), get(ENV, "OCORAMAN_RUN_ID", run_id("cabannes_sanity")))

summaries = Any[]
for wavelength in wavelengths
    summary = print_raman_air_summary(
        wavelength_nm=wavelength,
        temperature_K=temperature,
    )
    push!(summaries, Dict(string(k) => getfield(summary, k) for k in keys(summary)))
end

write_manifest(outdir;
    script=@__FILE__,
    extra=Dict(
        "wavelength_nm_list" => wavelengths,
        "temperature_K" => temperature,
    ))
write_json(joinpath(outdir, "cabannes_sanity.json"), summaries)
println("Wrote ", joinpath(outdir, "cabannes_sanity.json"))
