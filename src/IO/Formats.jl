# IO Formats and Sources registry

module Formats

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

function register_format(ext::AbstractString, parser::Function)
    _format_registry[string(lowercase(ext))] = parser
    return nothing
end

# Built-in YAML registration
register_format(".yaml", YAML.load_file)
register_format(".yml", YAML.load_file)

"Load a configuration Dict from a source"
function load_config(src::IOSource)
    if src isa FileSource
        path = (src::FileSource).path
        ext = lowercase(splitext(path)[2])
        @assert haskey(_format_registry, ext) "No parser registered for extension: $(ext)"
        return _format_registry[ext](path)
    elseif src isa DictSource
        return (src::DictSource).data
    else
        error("Unsupported IOSource: $(typeof(src))")
    end
end

end # module
