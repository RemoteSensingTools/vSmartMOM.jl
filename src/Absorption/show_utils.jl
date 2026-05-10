#=

Pretty-printing for Absorption module types

=#

function Base.show(io::IO, ::MIME"text/plain", x::HitranTable{FT}) where FT
    nlines = length(x.νᵢ)
    ν_lo = round(minimum(x.νᵢ), digits=2)
    ν_hi = round(maximum(x.νᵢ), digits=2)
    n_iso = length(unique(x.iso))
    n_mol = length(unique(x.mol))
    println(io, "HitranTable{$FT}")
    println(io, "├── lines: $nlines")
    println(io, "├── ν range: [$ν_lo, $ν_hi] cm⁻¹")
    println(io, "├── molecules: $n_mol unique ID(s)")
    print(io,   "└── isotopologues: $n_iso unique")
end

function Base.show(io::IO, x::HitranTable{FT}) where FT
    nlines = length(x.νᵢ)
    ν_lo = round(minimum(x.νᵢ), digits=1)
    ν_hi = round(maximum(x.νᵢ), digits=1)
    print(io, "HitranTable{$FT}($nlines lines, ν=[$ν_lo, $ν_hi] cm⁻¹)")
end
