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

"""
Load the generated active prior for one of the four surface classes.

Legacy retrieval priors contain 30 active coordinates.  New retrieval
experiments may deliberately fix additional entries and advertise their
solver dimension through the `active_state_count` global attribute.  Keeping
the active-to-full mapping in the prior file lets those experiments reuse the
same output and optimal-estimation machinery without pretending that a
zero-variance parameter is invertible.
"""
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
        expected_count = Int(get(
            dataset.attrib, "active_state_count", ACTIVE_STATE_COUNT))
        length(active_to_full) == expected_count || error(
            "prior advertises $expected_count active parameters but stores " *
            "$(length(active_to_full)) active indices")
        xa_full = Float64.(dataset["xa"][:, surface_index])
        xa = xa_full[active_to_full]
        Sa = Float64.(dataset["Sa_active"][:, :, surface_index])
        size(Sa) == (expected_count, expected_count) || error(
            "active covariance has size $(size(Sa)); expected " *
            "($expected_count,$expected_count)")
        names = retrieval_parameter_names(dataset, active_to_full)
        return RetrievalPrior(surface, xa, Sa, active_to_full, names)
    end
end

end # module RetrievalState
