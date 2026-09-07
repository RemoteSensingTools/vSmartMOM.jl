# Extract the two instrument mappings of each identical high-resolution output.
using JLD2, TOML, Statistics
root = ARGS[1]
metadata = TOML.parsefile(joinpath(root,"precision-frozen.toml"))
noise = sqrt.(JLD2.load(joinpath(root,
    "state035_corrected_siffalse-reference.jld2"),"variance")[1:934])
records = Dict[]
outputs = Dict()
for r in metadata["records"]
    file = r["file"]
    y = JLD2.load(joinpath(root,file),"y")
    mapped = JLD2.load(joinpath(root,replace(file,".jld2"=>"-instrument.jld2")),"instrument_y")
    @assert mapped[r["preparation"]] == y
    difference = (mapped["Float64"] - mapped["Float32"]) ./ noise
    push!(records,Dict("file"=>file,
        "max_noise_sigma"=>maximum(abs,difference),
        "rms_noise_sigma"=>sqrt(mean(abs2,difference))))
    outputs[(r["state"],r["preparation"],r["rt"])] = mapped
end
common = Dict[]
for state in ("reference","optimized"), instrument in ("Float32","Float64")
    a = outputs[(state,"Float32","Float64")][instrument]
    b = outputs[(state,"Float64","Float64")][instrument]
    difference = (b-a) ./ noise
    push!(common,Dict("state"=>state,"instrument"=>instrument,
        "max_noise_sigma"=>maximum(abs,difference),
        "rms_noise_sigma"=>sqrt(mean(abs2,difference))))
end
TOML.print(stdout,Dict("instrument_only"=>records,"preparation_common_instrument_RT64"=>common))
