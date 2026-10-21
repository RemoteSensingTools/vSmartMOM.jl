"""
    copy_parameters(params::vSmartMOM_Parameters; share_luts=false)

Copy a parameter template for an independent trial state. By default this is
equivalent to `deepcopy(params)`.

With `share_luts=true`, reuse the storage of loaded absorption tables in
`absorption_params.luts` and `absorption_params.h2o_lut`. The per-band LUT lists,
H₂O settings list, atmospheric profiles, VMR dictionary, aerosols, surfaces and
other parameter state are still deep-copied. Aliases within the copied state
are preserved, as with ordinary `deepcopy`.

Shared LUTs and everything they reference must remain read-only for the lifetime
of both templates. Replacing a copied list entry is safe; changing a shared
table's coefficients, grids or mutable metadata affects every owner. Model
construction and cross-section evaluation read the supplied tables without
modifying their storage. This option does not cache atmospheric optical depths
or change the behavior of `deepcopy` for any existing type.
"""
function copy_parameters(params::vSmartMOM_Parameters; share_luts::Bool=false)
    share_luts || return deepcopy(params)
    return deepcopy(_SharedLUTParameterCopy(params)).parameters
end

struct _SharedLUTParameterCopy{P}
    parameters::P
end

# Immutable LUT/interpolator wrappers are reconstructed by Base.deepcopy;
# memoizing only the outer LUT would not prevent their arrays being copied.
# Stop at mutable storage (especially large coefficient arrays), without
# walking its elements. Container lists are not passed to this helper.
function _share_lut_storage!(memo::IdDict, value)
    isbitstype(typeof(value)) && return
    if ismutable(value)
        memo[value] = value
    else
        for i in 1:fieldcount(typeof(value))
            isdefined(value,i) && _share_lut_storage!(memo,getfield(value,i))
        end
    end
    return nothing
end

# A private wrapper confines the custom deepcopy policy to this explicit
# operation. Base still copies the complete parameter graph with one memo,
# preserving shared references/cycles outside the read-only LUT payloads.
function Base.deepcopy_internal(wrapper::_SharedLUTParameterCopy, memo::IdDict)
    params = wrapper.parameters
    absorption = params.absorption_params
    if absorption !== nothing
        for band in absorption.luts, lut in band
            _share_lut_storage!(memo,lut)
        end
        for lut in absorption.h2o_lut
            _share_lut_storage!(memo,lut)
        end
    end
    return _SharedLUTParameterCopy(Base.deepcopy_internal(params,memo))
end
