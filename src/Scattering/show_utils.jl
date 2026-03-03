#=

Pretty-printing for Scattering module types (Oceananigans-inspired tree style)

=#

# Tree-drawing constants (matching CoreRT/tools/show_utils.jl)
const _S_TREE_MID   = "├── "
const _S_TREE_END   = "└── "
const _S_TREE_PIPE  = "│   "
const _S_TREE_SPACE = "    "

# ── Aerosol ───────────────────────────────────────────────────────────────

function Base.summary(io::IO, x::Aerosol)
    print(io, "Aerosol(", sprint(show, x.size_distribution), ", nᵣ=$(x.nᵣ), nᵢ=$(x.nᵢ))")
end

function Base.show(io::IO, ::MIME"text/plain", x::Aerosol)
    println(io, "Aerosol")
    println(io, _S_TREE_MID, "size_distribution: ", x.size_distribution)
    println(io, _S_TREE_MID, "nᵣ: ", x.nᵣ)
    print(io,   _S_TREE_END, "nᵢ: ", x.nᵢ)
end

function Base.show(io::IO, x::Aerosol)
    print(io, "Aerosol(nᵣ=$(x.nᵣ), nᵢ=$(x.nᵢ))")
end

# ── ScatteringMatrix ──────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", x::ScatteringMatrix)
    n = length(x.f₁₁)
    println(io, "ScatteringMatrix ($n angles)")
    println(io, _S_TREE_MID, "f₁₁: $n-element, range [$(round(minimum(x.f₁₁), sigdigits=3)), $(round(maximum(x.f₁₁), sigdigits=3))]")
    println(io, _S_TREE_MID, "f₁₂: $n-element")
    println(io, _S_TREE_MID, "f₂₂: $n-element")
    println(io, _S_TREE_MID, "f₃₃: $n-element")
    println(io, _S_TREE_MID, "f₃₄: $n-element")
    print(io,   _S_TREE_END, "f₄₄: $n-element")
end

function Base.show(io::IO, x::ScatteringMatrix)
    n = length(x.f₁₁)
    print(io, "ScatteringMatrix($n angles)")
end

# ── GreekCoefs ────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", x::GreekCoefs{FT}) where FT
    l_max = length(x.β) - 1
    println(io, "GreekCoefs{$FT} (l_max=$l_max)")
    println(io, _S_TREE_MID, "β: $(length(x.β))-element (phase function)")
    println(io, _S_TREE_MID, "α: $(length(x.α))-element")
    println(io, _S_TREE_MID, "γ: $(length(x.γ))-element")
    println(io, _S_TREE_MID, "δ: $(length(x.δ))-element")
    println(io, _S_TREE_MID, "ϵ: $(length(x.ϵ))-element")
    print(io,   _S_TREE_END, "ζ: $(length(x.ζ))-element")
end

function Base.show(io::IO, x::GreekCoefs{FT}) where FT
    l_max = length(x.β) - 1
    print(io, "GreekCoefs{$FT}(l_max=$l_max)")
end

# ── AerosolOptics ─────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", x::AerosolOptics{FT}) where FT
    l_max = length(x.greek_coefs.β) - 1
    println(io, "AerosolOptics{$FT}")
    println(io, _S_TREE_MID, "ω̃ (SSA): ", _scalar_or_range(x.ω̃))
    println(io, _S_TREE_MID, "k (extinction): ", _scalar_or_range(x.k))
    println(io, _S_TREE_MID, "fᵗ (truncation): ", _scalar_or_range(x.fᵗ))
    print(io,   _S_TREE_END, "greek_coefs: l_max=$l_max")
end

function Base.show(io::IO, x::AerosolOptics{FT}) where FT
    print(io, "AerosolOptics{$FT}(ω̃=", _scalar_or_range(x.ω̃), ")")
end

# ── MieModel ──────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", x::MieModel)
    comp_name = typeof(x.computation_type).name.name
    pol_name = typeof(x.polarization_type).name.name
    println(io, "MieModel ($comp_name)")
    println(io, _S_TREE_MID, "aerosol: ", x.aerosol)
    println(io, _S_TREE_MID, "λ: ", x.λ, " μm")
    println(io, _S_TREE_MID, "polarization: ", pol_name)
    println(io, _S_TREE_MID, "r_max: ", x.r_max, " μm")
    print(io,   _S_TREE_END, "nquad_radius: ", x.nquad_radius)
end

function Base.show(io::IO, x::MieModel)
    comp_name = typeof(x.computation_type).name.name
    print(io, "MieModel($comp_name, λ=$(x.λ) μm)")
end

# ── Polarization types ────────────────────────────────────────────────────

Base.show(io::IO, ::Stokes_I)    = print(io, "Stokes_I()")
Base.show(io::IO, ::Stokes_IQU)  = print(io, "Stokes_IQU()")
Base.show(io::IO, ::Stokes_IQUV) = print(io, "Stokes_IQUV()")

# ── Helpers ───────────────────────────────────────────────────────────────

function _scalar_or_range(x)
    if x isa Number
        return string(round(x, sigdigits=4))
    elseif x isa AbstractArray && length(x) > 0
        lo = round(minimum(x), sigdigits=4)
        hi = round(maximum(x), sigdigits=4)
        return lo == hi ? "$lo" : "[$lo, $hi] ($(length(x)) values)"
    else
        return string(x)
    end
end
