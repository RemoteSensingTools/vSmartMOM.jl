const _SS_MEASUREMENT_PATHS = (:total, :path1, :path2, :path3, :path4)

"""
    SSMeasurementSelector(; paths=(:total,), geometry_indices=:,
                            spectral_indices=:, stokes_indices=:)

Select and flatten StandaloneSS radiance outputs into a retrieval measurement
vector. Multiple paths are concatenated in the order provided. All Stokes
components are selected by default so future polarized outputs are not silently
reduced to Stokes I; pass `stokes_indices=1` for I-only retrieval vectors.
"""
struct SSMeasurementSelector{P,G,S,T}
    paths::P
    geometry_indices::G
    spectral_indices::S
    stokes_indices::T
end

function SSMeasurementSelector(; paths = (:total,),
                               geometry_indices = Colon(),
                               spectral_indices = Colon(),
                               stokes_indices = Colon())
    normalized_paths = _normalize_measurement_paths(paths)
    _validate_measurement_paths(normalized_paths)
    return SSMeasurementSelector(normalized_paths, geometry_indices,
                                 spectral_indices, stokes_indices)
end

_normalize_measurement_paths(path::Symbol) = (path,)
_normalize_measurement_paths(paths) = Tuple(paths)

function _validate_measurement_paths(paths)
    isempty(paths) &&
        throw(ArgumentError("SSMeasurementSelector requires at least one path"))
    for path in paths
        path in _SS_MEASUREMENT_PATHS ||
            throw(ArgumentError("unsupported StandaloneSS measurement path :$path; " *
                                "expected one of $_SS_MEASUREMENT_PATHS"))
    end
    return nothing
end

_slice_index(index::Integer) = index:index
_slice_index(index::Colon) = index
_slice_index(index) = index

function _selected_measurement_array(array, selector::SSMeasurementSelector)
    return array[_slice_index(selector.geometry_indices),
                 _slice_index(selector.stokes_indices),
                 _slice_index(selector.spectral_indices)]
end

function _selected_measurement_path(result, path::Symbol,
                                    selector::SSMeasurementSelector)
    hasproperty(result, path) ||
        throw(ArgumentError("result has no StandaloneSS path :$path"))
    return vec(_selected_measurement_array(getproperty(result, path),
                                           selector))
end

"""
    selected_measurements(result, selector=SSMeasurementSelector())

Return the selected StandaloneSS radiances as a one-dimensional measurement
vector. The vector order is geometry, then Stokes, then spectral point within
each selected path, with geometry varying fastest; multiple paths are
concatenated in selector order.
"""
function selected_measurements(result,
                               selector::SSMeasurementSelector =
                                   SSMeasurementSelector())
    chunks = map(path -> _selected_measurement_path(result, path, selector),
                 selector.paths)
    return length(chunks) == 1 ? first(chunks) : reduce(vcat, chunks)
end

selected_measurements(result; kwargs...) =
    selected_measurements(result, SSMeasurementSelector(; kwargs...))

function _validate_single_path_jacobian_selector(selector::SSMeasurementSelector)
    length(selector.paths) == 1 ||
        throw(ArgumentError("selected_measurement_jacobian for a single 4D " *
                            "Jacobian requires exactly one selected path"))
    return nothing
end

"""
    selected_measurement_jacobian(dL_dp, selector=SSMeasurementSelector())

Flatten a StandaloneSS Jacobian array with shape
`(nGeom, nStokes, nSpec, nParam)` into a retrieval Jacobian matrix with shape
`(nMeasurement, nParam)`, using the same measurement ordering as
`selected_measurements`.
"""
function selected_measurement_jacobian(
    dL_dp::AbstractArray{<:Any,4},
    selector::SSMeasurementSelector = SSMeasurementSelector(),
)
    _validate_single_path_jacobian_selector(selector)
    selected = dL_dp[_slice_index(selector.geometry_indices),
                     _slice_index(selector.stokes_indices),
                     _slice_index(selector.spectral_indices),
                     :]
    return reshape(selected, :, size(dL_dp, 4))
end

selected_measurement_jacobian(dL_dp; kwargs...) =
    selected_measurement_jacobian(dL_dp, SSMeasurementSelector(; kwargs...))
