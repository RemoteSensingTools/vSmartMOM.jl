using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering
using Statistics

const VLIDORT_BASELINE_DIR = @__DIR__
const REF_DIR     = joinpath(VLIDORT_BASELINE_DIR, "reference_data")
const CONFIG_DIR  = joinpath(VLIDORT_BASELINE_DIR, "configs")

"""
    rel_error(modeled, truth; near_zero_atol=1e-9)

Relative error matrix; positions where `|truth| < near_zero_atol` are masked
(returned as `NaN`) so they don't dominate the maximum.
"""
function rel_error(modeled::AbstractArray, truth::AbstractArray;
                   near_zero_atol::Real = 1e-9)
    size(modeled) == size(truth) || error("size mismatch: $(size(modeled)) vs $(size(truth))")
    out = similar(modeled, Float64)
    for i in eachindex(modeled, truth)
        if abs(truth[i]) < near_zero_atol
            out[i] = NaN
        else
            out[i] = abs(modeled[i] - truth[i]) / abs(truth[i])
        end
    end
    return out
end

"""
    summarize(label, modeled, truth; rtol, name="value")

Print summary statistics and return `(passes::Bool, max_rel)` for use with
`@test`. Skips NaN-masked entries.
"""
function summarize(label::AbstractString, modeled, truth; rtol::Real,
                   name::AbstractString = "value")
    re = rel_error(modeled, truth)
    valid = filter(!isnan, vec(re))
    if isempty(valid)
        @warn "$label: no comparable (non-zero-truth) entries"
        return (true, 0.0)
    end
    mx, mn = maximum(valid), median(valid)
    pass = mx < rtol
    flag = pass ? "PASS" : "FAIL"
    @info "$label [$flag]" name n=length(valid) median_rel=round(mn, sigdigits=3) max_rel=round(mx, sigdigits=3) rtol
    return (pass, mx)
end
