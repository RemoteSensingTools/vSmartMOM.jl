using vSmartMOM, JSON, CUDA

"Audit exported bindings and all registered docstrings in package modules."
function audit_api(output)
    modules = Module[]
    function visit(mod)
        mod in modules && return
        push!(modules, mod)
        for name in names(mod; all=true, imported=false)
            isdefined(mod, name) || continue
            value = getfield(mod, name)
            value isa Module || continue
            value !== mod && parentmodule(value) === mod && visit(value)
        end
    end
    visit(vSmartMOM)
    extensions = Dict{String,Bool}()
    for name in (:vSmartMOMCUDAExt, :vSmartMOMMetalExt)
        mod = Base.get_extension(vSmartMOM, name)
        extensions[string(name)] = mod !== nothing
        mod === nothing || visit(mod)
    end
    bindings = Any[]
    docstrings = Any[]
    failures = String[]
    for mod in sort(modules; by=string)
        for name in names(mod)
            name === nameof(mod) && continue
            defined = isdefined(mod, name)
            documented = defined && Base.Docs.doc(Base.Docs.Binding(mod, name)) !== nothing
            defined || push!(failures, "Undefined export: $mod.$name")
            documented || push!(failures, "Undocumented export: $mod.$name")
            push!(bindings, (; module_name=string(mod), name=string(name), defined,
                              documented, type=defined ? string(typeof(getfield(mod,name))) : "undefined"))
        end
        for (binding, multidoc) in Base.Docs.meta(mod; autoinit=false)
            for (signature, doc) in multidoc.docs
                # Parsing every registered docstring catches invalid attachment
                # data; strict Documenter separately resolves @ref and @docs.
                text = sprint(show, MIME"text/plain"(), Base.Docs.parsedoc(doc))
                isempty(strip(text)) && push!(failures, "Empty docstring: $binding $signature")
                path = string(get(doc.data, :path, ""))
                root = pkgdir(vSmartMOM)
                path = startswith(path, root) ? relpath(path, root) : path
                push!(docstrings, (; binding=string(binding), signature=string(signature),
                                   path, line=get(doc.data, :linenumber, 0)))
            end
        end
    end
    open(output, "w") do io
        JSON.print(io, (; modules=string.(modules), extensions, bindings, docstrings, failures), 2)
        println(io)
    end
    println("API audit: $(length(bindings)) exported bindings, $(length(docstrings)) registered docstrings, $(length(failures)) failures.")
    foreach(println, failures)
    isempty(failures) || error("API audit failed; see $output")
end

length(ARGS) == 1 || error("Usage: julia --project=test tools/audit_api.jl OUTPUT.json")
audit_api(only(ARGS))
