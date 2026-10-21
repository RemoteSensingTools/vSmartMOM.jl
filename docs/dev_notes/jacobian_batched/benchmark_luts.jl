# Optional caller-owned legacy HITRAN LUTs for the benchmark. Loading happens
# in fixture preparation, outside warmed model-construction/solve timings.
const benchmark_luts = Dict{String,Any}()

function benchmark_lut(molecule)
    get!(benchmark_luts,molecule) do
        path = joinpath(ENV["AUDIT_LUT_DIR"],molecule*".jld2")
        loaded = @timed vSmartMOM.Absorption.load_interpolation_model(path)
        lut = loaded.value
        @assert lut.mol == vSmartMOM.Absorption.mol_number(molecule)
        println("LUT_LOADED molecule=$molecule path=$path seconds=$(loaded.time) bytes=$(loaded.bytes) ",
                "ν=$(extrema(lut.ν_grid)) nν=$(length(lut.ν_grid)) ",
                "p=$(extrema(lut.p_grid)) T=$(extrema(lut.t_grid)) isotope=$(lut.iso)")
        flush(stdout)
        lut
    end
end

function set_benchmark_luts!(p,gases)
    haskey(ENV,"AUDIT_LUT_DIR") || return p
    isempty(gases) && error("AUDIT_LUT_DIR requires AUDIT_GASES")
    dry = filter(!=("H2O"),gases)
    p.absorption_params.luts = Any[[benchmark_lut(mol) for mol in dry]]
    "H2O" in gases && (p.absorption_params.h2o_lut[1] = benchmark_lut("H2O"))
    for mol in gases
        lut = benchmark_lut(mol)
        lo,hi = extrema(p.spec_bands[1])
        @assert first(lut.ν_grid) <= lo <= hi <= last(lut.ν_grid) "LUT spectral coverage: $mol"
    end
    return p
end
