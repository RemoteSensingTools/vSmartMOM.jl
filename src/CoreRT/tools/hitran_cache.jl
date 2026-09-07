# Parsed line databases contain spectroscopy only, independent of the scene,
# spectral query, broadening settings and CPU/GPU backend. Treat their arrays
# as read-only; each solve still constructs its own LineByLineModel wrapper.
const _HITRAN_LINES_LOCK = ReentrantLock()
const _HITRAN_LINES = Dict{Tuple{String,DataType},Any}()

"""
    _hitran_lines(molecule, FT)

Load a HITRAN file once per process and floating-point type. The resolved path
separates artifacts and downloaded editions, so switching editions selects a
different key. Call `clear_spectroscopy_cache!` after replacing a downloaded
file in place; automatic reloads are deliberately avoided.
The lock covers the first parse as well as lookup: concurrent constructors do
not independently parse the same file. No spectral clipping is introduced;
line wings and pressure shifts retain the existing absorption-model semantics.
"""
function _hitran_lines(molecule, ::Type{FT}) where {FT<:AbstractFloat}
    return lock(_HITRAN_LINES_LOCK) do
        path = artifact(molecule)
        get!(_HITRAN_LINES, (path,FT)) do
            AtmosphericAbsorption.load_lines(AtmosphericAbsorption.HitranPort(path); FT)
        end::AtmosphericAbsorption.LineDatabase{FT}
    end
end

"""
    clear_spectroscopy_cache!()

Release the process-local parsed HITRAN line cache. Existing models remain
valid; later constructors parse their data again. This does not delete any
downloaded files or clear caller-owned LUTs/BatchContexts. Normal repeated
forward and linearized construction reuses cached data without this call.
Use it after replacing a downloaded spectroscopy file at the same path, or
to release retained databases when they are no longer needed.
"""
function clear_spectroscopy_cache!()
    lock(_HITRAN_LINES_LOCK) do
        empty!(_HITRAN_LINES)
    end
    return nothing
end
