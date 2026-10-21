module Round5FixedSIF

using ..OptimalEstimation: ForwardEvaluation

export Round5SIFMap,
       Round5ForwardEvaluator,
       CORE_STATE_COUNT,
       CORE_SIF760_INDEX,
       CORE_MSIF_INDEX,
       ROUND5_STATE_COUNT,
       NU_759_CM1,
       NU_760_CM1,
       DELTA_NU_760_MINUS_759_CM1,
       active_state_count,
       active_core_indices,
       expand_round5_state,
       reduce_core_jacobian

const WAVENUMBER_CONVERSION_NM_CM1 = 1.0e7
const NU_759_CM1 = WAVENUMBER_CONVERSION_NM_CM1 / 759.0
const NU_760_CM1 = WAVENUMBER_CONVERSION_NM_CM1 / 760.0
const DELTA_NU_760_MINUS_759_CM1 = NU_760_CM1 - NU_759_CM1

# The validated OCO_RRS_synth core order is fixed:
# 1:28 = pressure, active CO2, aerosol, and surface coordinates;
# 29 = SIF760; 30 = mSIF. Round 5 removes both SIF coordinates from the
# numerical state while retaining the unchanged core forward/Jacobian model.
const CORE_STATE_COUNT = 30
const CORE_SIF760_INDEX = 29
const CORE_MSIF_INDEX = 30
const ROUND5_STATE_COUNT = 28

"""
    Round5SIFMap(sif_on; Lnu759=0, mSIF=0)

Fixed SIF boundary map for retrieval round 5. For a SIF-on scene, both the
spectral radiance at 759 nm and the native wavenumber slope are prescribed.
The 760-nm core coefficient is derived exactly from the same linear source
model used by vSmartMOM:

```text
SIF760 = Lnu759 + mSIF * (nu760 - nu759).
```

For a SIF-off scene both coefficients must be exactly zero. Neither mode has
an active SIF coordinate.
"""
struct Round5SIFMap{T<:AbstractFloat}
    sif_on::Bool
    Lnu759::T
    mSIF::T
end

function Round5SIFMap(sif_on::Bool;
                      Lnu759::Real=0.0,
                      mSIF::Real=0.0)
    anchor = Float64(Lnu759)
    slope = Float64(mSIF)
    all(isfinite, (anchor, slope)) || throw(ArgumentError(
        "fixed round-5 SIF coefficients must be finite"))
    anchor >= 0 || throw(ArgumentError(
        "fixed round-5 Lnu759 must be nonnegative"))
    !sif_on && (!iszero(anchor) || !iszero(slope)) && throw(ArgumentError(
        "a SIF-off round-5 map requires Lnu759=mSIF=0"))
    return Round5SIFMap{Float64}(sif_on, anchor, slope)
end

active_state_count(::Round5SIFMap) = ROUND5_STATE_COUNT
active_core_indices(::Round5SIFMap) = collect(1:ROUND5_STATE_COUNT)

"""Expand a 28-coordinate round-5 state into the 30-coordinate core state."""
function expand_round5_state(map::Round5SIFMap,
                             state::AbstractVector{<:Real})
    length(state) == ROUND5_STATE_COUNT || throw(DimensionMismatch(
        "round-5 state has $(length(state)) entries; expected " *
        string(ROUND5_STATE_COUNT)))
    all(isfinite, state) || throw(ArgumentError(
        "round-5 state contains a non-finite value"))

    FT = promote_type(Float64, eltype(state), typeof(map.Lnu759))
    core = zeros(FT, CORE_STATE_COUNT)
    @views core[1:ROUND5_STATE_COUNT] .= state
    if map.sif_on
        core[CORE_SIF760_INDEX] = FT(map.Lnu759) + FT(map.mSIF) *
            FT(DELTA_NU_760_MINUS_759_CM1)
        core[CORE_MSIF_INDEX] = FT(map.mSIF)
    end
    return core
end

"""Drop the two fixed SIF columns from a core OCO_RRS_synth Jacobian."""
function reduce_core_jacobian(::Round5SIFMap,
                              jacobian::AbstractMatrix{<:Real})
    size(jacobian, 2) == CORE_STATE_COUNT || throw(DimensionMismatch(
        "core Jacobian has $(size(jacobian, 2)) columns; expected " *
        string(CORE_STATE_COUNT)))
    all(isfinite, jacobian) || throw(ArgumentError(
        "core Jacobian contains a non-finite value"))
    return Matrix(jacobian[:, 1:ROUND5_STATE_COUNT])
end

"""
Boundary wrapper around the validated 30-column OCO_RRS_synth evaluator.
The fixed SIF coefficients enter the forward calculation, but their two
Jacobian columns are omitted from the round-5 retrieval state.
"""
struct Round5ForwardEvaluator{E,M<:Round5SIFMap}
    core_evaluator::E
    map::M
end

function (evaluator::Round5ForwardEvaluator)(state::AbstractVector)
    core_state = expand_round5_state(evaluator.map, state)
    evaluation = evaluator.core_evaluator(core_state)
    evaluation isa ForwardEvaluation || throw(ArgumentError(
        "round-5 core evaluator must return ForwardEvaluation"))
    return ForwardEvaluation(
        evaluation.measurement,
        reduce_core_jacobian(evaluator.map, evaluation.jacobian),
        evaluation.band_ranges;
        timing=evaluation.timing)
end

end # module Round5FixedSIF
