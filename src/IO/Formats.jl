# IO Formats and Sources registry

module Formats

using TOML
using YAML

export IOSource, FileSource, DictSource, load_config, register_format

abstract type IOSource end

struct FileSource <: IOSource
    path::String
end

struct DictSource <: IOSource
    data::Dict
end

# Simple format registry based on file extension
const _format_registry = Dict{String, Function}()

"""
    _format_error(msg)

Raise a stable `ArgumentError` for invalid IO format configuration.
"""
@inline _format_error(msg) = throw(ArgumentError(msg))

"""
    _require_format(cond, msg)

Validate IO format registry input and raise `ArgumentError` when invalid.
"""
@inline _require_format(cond, msg) = cond ? nothing : _format_error(msg)

function register_format(ext::AbstractString, parser::Function)
    _format_registry[string(lowercase(ext))] = parser
    return nothing
end

# TOML parses a mixed integer/float array such as `[60, 45.5]` as
# `Vector{Any}`; YAML gives `Vector{Real}`. Normalize all-numeric arrays the
# same way so that both formats pass the parser's `Array{<:Real}` checks.
_normalize_toml(x) = x
_normalize_toml(d::AbstractDict) = Dict{String,Any}(String(k) => _normalize_toml(v) for (k, v) in d)
# Homogeneous arrays keep their type; arrays of tables are normalized inside.
_normalize_toml(v::AbstractVector) =
    eltype(v) <: AbstractDict ? [_normalize_toml(x) for x in v] : v
function _normalize_toml(v::Vector{Any})
    w = Any[_normalize_toml(x) for x in v]
    return !isempty(w) && all(x -> x isa Real && !(x isa Bool), w) ? Real[w...] : w
end
_parse_toml_file(path) = _normalize_toml(TOML.parsefile(path))

# Built-in format registration
register_format(".toml", _parse_toml_file)
register_format(".yaml", YAML.load_file)
register_format(".yml", YAML.load_file)

"Load a configuration Dict from a source"
function load_config(src::FileSource)
    path = src.path
    ext = lowercase(splitext(path)[2])
    _require_format(haskey(_format_registry, ext), "No parser registered for extension: $(ext)")
    return _format_registry[ext](path)
end

load_config(src::DictSource) = src.data

load_config(src::IOSource) = _format_error("Unsupported IOSource: $(typeof(src))")

end # module
