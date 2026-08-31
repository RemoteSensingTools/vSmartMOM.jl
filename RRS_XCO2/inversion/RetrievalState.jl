module RetrievalState

using NCDatasets

export RetrievalPrior,
       ACTIVE_STATE_COUNT,
       load_retrieval_prior,
       retrieval_parameter_names

const ACTIVE_STATE_COUNT = 30
const DEFAULT_PRIOR_PATH = joinpath(
    @__DIR__, "retrieval_setup", "apriori_states.nc")

"""One surface-specific active retrieval prior and its full-state mapping."""
struct RetrievalPrior
    surface::Symbol
    xa::Vector{Float64}
    Sa::Matrix{Float64}
    active_to_full::Vector{Int}
    parameter_names::Vector{String}
end

function retrieval_parameter_names(dataset, active_to_full)
    full_names = split(String(dataset.attrib["parameter_names"]))
    maximum(active_to_full) <= length(full_names) || error(
        "active-state index exceeds the stored parameter-name list")
    return full_names[active_to_full]
end

"""Load the generated 30-element prior for one of the four surface classes."""
function load_retrieval_prior(surface::Symbol;
                              path::AbstractString=DEFAULT_PRIOR_PATH)
    isfile(path) || throw(ArgumentError(
        "missing generated prior $path; run retrieval_setup/build_apriori.jl"))
    return NCDataset(path) do dataset
        get(dataset.attrib, "apriori_complete", 0) == 1 || error(
            "prior file is not marked complete: $path")
        surfaces = Symbol.(split(String(dataset.attrib["surface_order"])))
        surface_index = findfirst(==(surface), surfaces)
        isnothing(surface_index) && throw(ArgumentError(
            "surface $surface is absent from $path"))
        active_to_full = Int.(dataset["active_parameter_index"][:])
        length(active_to_full) == ACTIVE_STATE_COUNT || error(
            "expected $ACTIVE_STATE_COUNT active parameters")
        xa_full = Float64.(dataset["xa"][:, surface_index])
        xa = xa_full[active_to_full]
        Sa = Float64.(dataset["Sa_active"][:, :, surface_index])
        names = retrieval_parameter_names(dataset, active_to_full)
        return RetrievalPrior(surface, xa, Sa, active_to_full, names)
    end
end

end # module RetrievalState
