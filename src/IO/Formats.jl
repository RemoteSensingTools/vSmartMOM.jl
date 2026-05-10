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

# Built-in YAML registration
register_format(".toml", TOML.parsefile)
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
