module Round4KnownSIF

using ..OptimalEstimation: ForwardEvaluation

export Round4SIFMap,
       Round4ForwardEvaluator,
       CORE_STATE_COUNT,
       CORE_SIF760_INDEX,
       CORE_MSIF_INDEX,
       ROUND4_BASE_STATE_COUNT,
       NU_759_CM1,
       NU_760_CM1,
       DELTA_NU_760_MINUS_759_CM1,
       active_state_count,
       active_core_indices,
       expand_round4_state,
       reduce_core_jacobian

const WAVENUMBER_CONVERSION_NM_CM1 = 1.0e7
const NU_759_CM1 = WAVENUMBER_CONVERSION_NM_CM1 / 759.0
const NU_760_CM1 = WAVENUMBER_CONVERSION_NM_CM1 / 760.0
const DELTA_NU_760_MINUS_759_CM1 = NU_760_CM1 - NU_759_CM1

# `OCO_RRS_synth` currently exposes this stable global order:
# 1:28 = pressure, active CO2, aerosol, and surface coordinates;
# 29 = SIF760; 30 = mSIF.  Round 4 deliberately transforms only this
# retrieval boundary.  The already validated RT/source Jacobian remains in its
# original 30-column basis.
const CORE_STATE_COUNT = 30
const CORE_SIF760_INDEX = 29
const CORE_MSIF_INDEX = 30
const ROUND4_BASE_STATE_COUNT = 28

"""
    Round4SIFMap(sif_on; Lν759=0)

State-coordinate map for the round-4 retrieval in which the SIF spectral
radiance at 759 nm is known exactly. `Lν759` has the native vSmartMOM units
`mW m⁻² sr⁻¹ (cm⁻¹)⁻¹`.

For a SIF-on scene, the only active SIF coordinate is the native wavenumber
slope `mSIF = dLν/dν`. The core source still consumes `(SIF760,mSIF)`, so

```
SIF760 = Lν759 + mSIF * (ν760 - ν759).
```

For a SIF-off scene both core SIF coordinates are fixed to exact zero and the
solver has no active SIF coordinate.
"""
struct Round4SIFMap{T<:AbstractFloat}
    sif_on::Bool
    Lν759::T
end

function Round4SIFMap(sif_on::Bool; Lν759::Real=0.0)
    value = Float64(Lν759)
    isfinite(value) || throw(ArgumentError("known Lν759 must be finite"))
    value >= 0 || throw(ArgumentError("known Lν759 must be nonnegative"))
    !sif_on && !iszero(value) && throw(ArgumentError(
        "a SIF-off round-4 map requires Lν759=0"))
    return Round4SIFMap{Float64}(sif_on, value)
end

active_state_count(map::Round4SIFMap) =
    ROUND4_BASE_STATE_COUNT + Int(map.sif_on)

function active_core_indices(map::Round4SIFMap)
    indices = collect(1:ROUND4_BASE_STATE_COUNT)
    map.sif_on && push!(indices, CORE_MSIF_INDEX)
    return indices
end

"""Expand a 28/29-coordinate round-4 state into the core 30-coordinate state."""
function expand_round4_state(map::Round4SIFMap,
                             state::AbstractVector{<:Real})
    expected = active_state_count(map)
    length(state) == expected || throw(DimensionMismatch(
        "round-4 state has $(length(state)) entries; expected $expected"))
    all(isfinite, state) || throw(ArgumentError(
        "round-4 state contains a non-finite value"))

    FT = promote_type(Float64, eltype(state), typeof(map.Lν759))
    core = zeros(FT, CORE_STATE_COUNT)
    @views core[1:ROUND4_BASE_STATE_COUNT] .= state[1:ROUND4_BASE_STATE_COUNT]
    if map.sif_on
        mSIF = FT(state[end])
        core[CORE_SIF760_INDEX] = FT(map.Lν759) +
            mSIF * FT(DELTA_NU_760_MINUS_759_CM1)
        core[CORE_MSIF_INDEX] = mSIF
    end
    return core
end

"""
Transform an instrument-space core Jacobian to the round-4 active basis.

For SIF-on scenes the final column is the exact chain rule
`K_m = K_SIF760 * (ν760-ν759) + K_mSIF`. SIF-off scenes retain only
columns 1:28.
"""
function reduce_core_jacobian(map::Round4SIFMap,
                              jacobian::AbstractMatrix{<:Real})
    size(jacobian, 2) == CORE_STATE_COUNT || throw(DimensionMismatch(
        "core Jacobian has $(size(jacobian, 2)) columns; expected " *
        string(CORE_STATE_COUNT)))
    all(isfinite, jacobian) || throw(ArgumentError(
        "core Jacobian contains a non-finite value"))

    result = Matrix{promote_type(Float64, eltype(jacobian))}(
        undef, size(jacobian, 1), active_state_count(map))
    @views result[:, 1:ROUND4_BASE_STATE_COUNT] .=
        jacobian[:, 1:ROUND4_BASE_STATE_COUNT]
    if map.sif_on
        @views result[:, end] .=
            jacobian[:, CORE_SIF760_INDEX] .*
                DELTA_NU_760_MINUS_759_CM1 .+
            jacobian[:, CORE_MSIF_INDEX]
    end
    return result
end

"""
Callable boundary wrapper around the existing 30-column `OCO_RRS_synth`
evaluator. The wrapped evaluator remains unchanged; this object expands the
round-4 state before evaluation and applies the exact Jacobian chain rule on
return.
"""
struct Round4ForwardEvaluator{E,M<:Round4SIFMap}
    core_evaluator::E
    map::M
end

function (evaluator::Round4ForwardEvaluator)(state::AbstractVector)
    core_state = expand_round4_state(evaluator.map, state)
    evaluation = evaluator.core_evaluator(core_state)
    evaluation isa ForwardEvaluation || throw(ArgumentError(
        "round-4 core evaluator must return ForwardEvaluation"))
    jacobian = reduce_core_jacobian(evaluator.map, evaluation.jacobian)
    return ForwardEvaluation(
        evaluation.measurement, jacobian, evaluation.band_ranges;
        timing=evaluation.timing)
end

end # module Round4KnownSIF
