#=
 
This file contains a read_hitran function to read through a HITRAN data file and 
produce a HitranTable struct. 
 
=#

"""
    read_hitran(filepath; mol=-1, iso=-1, ν_min=0, ν_max=Inf, min_strength=0)

Read and parse a HITRAN line-by-line data file into a HitranTable.

Filters lines by molecule, isotopologue, wavenumber range, and minimum line strength.
Uses fixed-width column parsing per the HITRAN format specification.

# Arguments
- `filepath::String`: Path to HITRAN .par or .hitran file
- `mol::Int=-1`: Filter by molecule ID (-1 = all)
- `iso::Int=-1`: Filter by isotopologue ID (-1 = all)
- `ν_min::Real=0`: Minimum wavenumber (cm⁻¹)
- `ν_max::Real=Inf`: Maximum wavenumber (cm⁻¹)
- `min_strength::Real=0`: Minimum line intensity at 296 K

# Returns
- `HitranTable`: Struct with line parameters (νᵢ, Sᵢ, γ_air, γ_self, E″, etc.)

# Throws
- `HitranEmptyError`: If no lines match the filter criteria
- `ArgumentError`: If a record is short, non-ASCII, or contains a malformed
  required numeric field. Blank upper/lower degeneracies are accepted as zero.
"""
function read_hitran(filepath::String; mol::Int=-1, iso::Int=-1, 
                     ν_min::Real=0, ν_max::Real=Inf, 
                     min_strength::Real=0)
    FT = typeof(float(ν_min))
    columns = (
        mol=Int[], iso=Int[], νᵢ=FT[], Sᵢ=FT[], Aᵢ=FT[], γ_air=FT[],
        γ_self=FT[], E″=FT[], n_air=FT[], δ_air=FT[],
        global_upper_quanta=String[], global_lower_quanta=String[],
        local_upper_quanta=String[], local_lower_quanta=String[],
        ierr=String[], iref=String[], line_mixing_flag=String[],
        g′=FT[], g″=FT[])

    open(filepath, "r") do file
        for (line_number, record) in enumerate(eachline(file))
            fields = _hitran_fields(record, filepath, line_number)
            molecule = _parse_hitran_number(Int, fields[1], filepath, line_number, "molec_id")
            isotopologue = _parse_hitran_isotopologue(
                fields[2], filepath, line_number)
            νᵢ = _parse_hitran_number(FT, fields[3], filepath, line_number, "nu")
            Sᵢ = _parse_hitran_number(FT, fields[4], filepath, line_number, "sw")

            (mol == -1 || molecule == mol) || continue
            (iso == -1 || isotopologue == iso) || continue
            ν_min <= νᵢ <= ν_max || continue
            Sᵢ >= min_strength || continue

            push!(columns.mol, molecule)
            push!(columns.iso, isotopologue)
            push!(columns.νᵢ, νᵢ)
            push!(columns.Sᵢ, Sᵢ)
            push!(columns.Aᵢ, _parse_hitran_number(FT, fields[5], filepath, line_number, "a"))
            push!(columns.γ_air, _parse_hitran_number(FT, fields[6], filepath, line_number, "gamma_air"))
            push!(columns.γ_self, _parse_hitran_number(FT, fields[7], filepath, line_number, "gamma_self"))
            push!(columns.E″, _parse_hitran_number(FT, fields[8], filepath, line_number, "elower"))
            push!(columns.n_air, _parse_hitran_number(FT, fields[9], filepath, line_number, "n_air"))
            push!(columns.δ_air, _parse_hitran_number(FT, fields[10], filepath, line_number, "delta_air"))
            push!(columns.global_upper_quanta, fields[11])
            push!(columns.global_lower_quanta, fields[12])
            push!(columns.local_upper_quanta, fields[13])
            push!(columns.local_lower_quanta, fields[14])
            push!(columns.ierr, fields[15])
            push!(columns.iref, fields[16])
            push!(columns.line_mixing_flag, fields[17])
            push!(columns.g′, _parse_hitran_number(FT, fields[18], filepath, line_number, "gp"; blank_is_zero=true))
            push!(columns.g″, _parse_hitran_number(FT, fields[19], filepath, line_number, "gpp"; blank_is_zero=true))
        end
    end

    isempty(columns.mol) && throw(HitranEmptyError())
    return HitranTable(; columns...)
end

const _HITRAN_FIELD_WIDTHS = (2, 1, 12, 10, 10, 5, 5, 10, 4, 8,
                              15, 15, 15, 15, 6, 12, 1, 7, 7)
const _HITRAN_FIELD_ENDS = Tuple(cumsum(collect(_HITRAN_FIELD_WIDTHS)))
const _HITRAN_RECORD_LENGTH = last(_HITRAN_FIELD_ENDS)

function _hitran_fields(record::AbstractString, filepath, line_number)
    isascii(record) || throw(ArgumentError(
        "$filepath:$line_number: HITRAN records must contain ASCII text"))
    ncodeunits(record) >= _HITRAN_RECORD_LENGTH || throw(ArgumentError(
        "$filepath:$line_number: short HITRAN record ($(ncodeunits(record)) bytes; " *
        "expected at least $_HITRAN_RECORD_LENGTH)"))
    bytes = codeunits(record)
    return ntuple(length(_HITRAN_FIELD_WIDTHS)) do index
        first_byte = index == 1 ? 1 : _HITRAN_FIELD_ENDS[index - 1] + 1
        String(view(bytes, first_byte:_HITRAN_FIELD_ENDS[index]))
    end
end

function _parse_hitran_number(::Type{T}, field, filepath, line_number, name;
                              blank_is_zero::Bool=false) where {T<:Real}
    text = strip(field)
    isempty(text) && blank_is_zero && return zero(T)
    value = tryparse(T, text)
    value === nothing && throw(ArgumentError(
        "$filepath:$line_number: invalid HITRAN $name field $(repr(field))"))
    return value
end

function _parse_hitran_isotopologue(field, filepath, line_number)
    text = strip(field)
    ncodeunits(text) == 1 || throw(ArgumentError(
        "$filepath:$line_number: invalid HITRAN local_iso_id field $(repr(field))"))
    code = only(codeunits(text))
    if UInt8('1') <= code <= UInt8('9')
        return Int(code - UInt8('0'))
    elseif code == UInt8('0')
        return 10
    elseif UInt8('A') <= code <= UInt8('Z')
        return 11 + Int(code - UInt8('A'))
    end
    throw(ArgumentError(
        "$filepath:$line_number: invalid HITRAN local_iso_id field $(repr(field))"))
end
