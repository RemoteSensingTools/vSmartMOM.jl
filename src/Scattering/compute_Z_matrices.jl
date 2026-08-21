# =====================================================================
# Phase-matrix Fourier moments Z⁺⁺(m), Z⁻⁺(m) for azimuthal mode m,
# evaluated on the discrete-ordinate stream cosines μ.
# =====================================================================
# Given the six Greek expansion coefficients α, β, γ, δ, ϵ, ζ that
# represent the polarized phase matrix (Hovenier convention — the four
# Stokes channels are I, Q, U, V), this routine returns the matrices
# Z⁺⁺ and Z⁻⁺ that multiply the upward / downward-going Stokes vectors
# in the vector RT equation.
#
#     Z(μ, μ′; m) = Σₗ ωₗ Pₗ(m, μ, μ′)         [Hovenier eq. 2.66, 2.69]
#
# where ωₗ are linear combinations of α…ζ and Pₗ are the generalised
# spherical functions (associated Legendre + Wigner-d). Z⁺⁺ takes
# (μ, μ′) on the same hemisphere; Z⁻⁺ flips the sign of one argument.
# The single-scattering albedo and τ are NOT folded in here — those
# enter elemental! / doubling!.
#
# Reference: Hovenier, van der Mee & Domke (2004), "Transfer of
# Polarized Light in Planetary Atmospheres", Ch. 2-3.
# =====================================================================
# ---------------------------------------------------------------------
# Per-solve precompute tables for compute_Z_moments.
#
# The expensive inputs of Z(μ, μ′; m) factor cleanly by what they depend on:
#
#   P_l^m, R_l^m, T_l^m  — generalised spherical functions; depend ONLY on
#                          (μ, l). One recurrence sweep produces the full
#                          (l, m) triangle for all Fourier moments at once.
#   Π_l^m                — assembled from P/R/T slices; depends on (l, m).
#   B_l                  — linear combinations of the Greek coefficients
#                          α…ζ; the ONLY factor that depends on the
#                          scatterer (Rayleigh vs each aerosol mode).
#
# The plain compute_Z_moments recomputes ALL of these per call, so an RT
# solve with nModes scatterers and m = 0:m_max pays the P/R/T recurrences
# 2·nModes·(m_max+1) times (±μ) and the Π assembly nModes·(m_max+1) times,
# even though both are scatterer-independent. ZMomentTables holds P/R/T
# once per solve; make_Π_lists shares Π across all scatterers of one m.
# ---------------------------------------------------------------------
"""
    ZMomentTables(μ, l_max)

Scatterer-independent precompute for [`compute_Z_moments`](@ref): the
generalised spherical functions ``P_l^m, R_l^m, T_l^m`` (Sanghavi 2014,
eq. 15) evaluated on the stream cosines ``±μ`` up to a global `l_max`
(the longest Greek-coefficient expansion among all scatterers in the run).
Build once per solve and pass via the `tables` keyword.
"""
struct ZMomentTables{FT}
    μ     :: Vector{FT}
    l_max :: Int
    # (nμ, l_max, l_max) arrays over (μ, l, m); the ⁻ set is evaluated at -μ.
    P  :: Array{FT,3};  R  :: Array{FT,3};  T  :: Array{FT,3}
    P⁻ :: Array{FT,3};  R⁻ :: Array{FT,3};  T⁻ :: Array{FT,3}
end

function ZMomentTables(μ::AbstractVector, l_max::Int)
    @assert all(0 .< μ .≤ 1) "all μ's within ZMomentTables have to be ∈ ]0,1]"
    P,  R,  T  = compute_associated_legendre_PRT(μ,  l_max)
    P⁻, R⁻, T⁻ = compute_associated_legendre_PRT(-μ, l_max)
    return ZMomentTables(collect(μ), l_max, P, R, T, P⁻, R⁻, T⁻)
end

"""
    make_Π_lists(mod, tables, m)

Assemble the Π matrices (Sanghavi 2014, eq. 15) for Fourier moment `m`
(0-indexed) from precomputed `tables`, for every l = m+1 … l_max (1-indexed).
The result is shared by ALL scatterers at this `m` — Π does not involve the
Greek coefficients. Returns `(Π, Π⁻)` indexed as `Π[l]`, or `nothing` when
`m ≥ l_max` (all Z moments vanish there).
"""
function make_Π_lists(mod::AbstractPolarizationType, t::ZMomentTables, m::Int)
    m1 = m + 1                       # 1-indexed Fourier moment
    m1 > t.l_max && return nothing
    Π1  = construct_Π_matrix(mod, t.P,  t.R,  t.T,  m1, m1)
    Π⁻1 = construct_Π_matrix(mod, t.P⁻, t.R⁻, t.T⁻, m1, m1)
    Π  = Vector{typeof(Π1)}(undef, t.l_max);   Π[m1]  = Π1
    Π⁻ = Vector{typeof(Π⁻1)}(undef, t.l_max);  Π⁻[m1] = Π⁻1
    for l in (m1 + 1):t.l_max
        Π[l]  = construct_Π_matrix(mod, t.P,  t.R,  t.T,  l, m1)
        Π⁻[l] = construct_Π_matrix(mod, t.P⁻, t.R⁻, t.T⁻, l, m1)
    end
    return (Π, Π⁻)
end

"""
    $(FUNCTIONNAME)(mod::AbstractPolarizationType, μ, α, β, γ, δ, ϵ, ζ, m::Int)
Compute moments of the phase matrix.

Keywords `tables`/`Π_pair` (see [`ZMomentTables`](@ref) / [`make_Π_lists`](@ref))
supply the scatterer-independent precompute; results are identical (the same
sums in the same order), only the redundant recomputation is skipped.
"""
function compute_Z_moments(mod::AbstractPolarizationType, μ, greek_coefs::GreekCoefs, m::Int ;
                           arr_type = Array,
                           tables::Union{Nothing, ZMomentTables} = nothing,
                           Π_pair = nothing)
    if tables !== nothing
        # The tables are only valid for the μ they were built on — a mismatch
        # would silently evaluate the phase matrix at the wrong stream angles.
        length(μ) == length(tables.μ) && all(μ .== tables.μ) ||
            throw(ArgumentError("compute_Z_moments: `tables` were built for a " *
                                "different μ than the one supplied"))
        # Callers that batch over scatterers pass the per-m Π lists explicitly
        # (built once, shared); a lone call may omit them — build them here
        # rather than silently skipping the accumulation (Z ≡ 0 bug).
        Π_pair === nothing && (Π_pair = make_Π_lists(mod, tables, m))
        return _compute_Z_moments_tabulated(mod, μ, greek_coefs, m,
                                            arr_type, tables, Π_pair)
    end
    (; α, β, γ, δ, ϵ, ζ) = greek_coefs
    FT = eltype(β)
    n = length(μ)

    # Set prefactor for moments (note 1-notation for `m` here):
    fact = (m == 0) ? 0.5 : 1.0

    # Change from 0-index to 1-index (i.e. the lowest m is 0 ), 
    # make more logical later to avoid confusion later (m=0 has a meaning!!)
    m = m+1
    
    # get l_max just from length of array:
    l_max = length(β)

    # Check that all μ are positive here ([0,1])
    # @show μ
    @assert all(0 .< μ .≤ 1) "all μ's within compute_Z_moments have to be ∈ ]0,1]"

    # Compute legendre Polynomials at μ and up to lmax
    P, R, T    = compute_associated_legendre_PRT(μ, l_max)
    P⁻, R⁻, T⁻ = compute_associated_legendre_PRT(-μ, l_max)
  
    # Pre-compute all required B matrices
    𝐁_all = [construct_B_matrix(mod, α, β, γ, δ, ϵ, ζ, i) for i in 1:l_max]
    # Get dimension of square matrix (easier for Scalar/Stokes dimensions)
    B_dim = Int(sqrt(length(𝐁_all[1])))
    
    # Create matrices:
    nb = B_dim * n
    𝐙⁺⁺, 𝐙⁻⁺ = (zeros(FT, nb, nb), zeros(FT, nb, nb))
    A⁺⁺, A⁻⁺ = (zeros(FT, B_dim, B_dim, n, n), zeros(FT, B_dim, B_dim, n, n))

    # Iterate over l
    for l = m:l_max

        # B matrix for l
        𝐁 = 𝐁_all[l];

        # Construct Π matrix for l,m pair (change to in place later!)
        # See eq. 15 in Sanghavi 2014, note that P,R,T are already normalized
        Π  = construct_Π_matrix(mod, P, R, T, l, m)
        Π⁻ = construct_Π_matrix(mod, P⁻, R⁻, T⁻, l, m)
        # Iterate over angles
        for j in eachindex(μ), i in eachindex(μ)
            if B_dim == 1
                A⁺⁺[B_dim,B_dim,i,j] += Π[i] * 𝐁 * Π[j]
                A⁻⁺[B_dim,B_dim,i,j] += Π[i] * 𝐁 * Π⁻[j]
            else
                A⁺⁺[:,:,i,j] += Π[i] * 𝐁 * Π[j]
                A⁻⁺[:,:,i,j] += Π[i] * 𝐁 * Π⁻[j]
            end
        end
    end

    # Now get to the Z part:
    for imu in eachindex(μ), jmu in eachindex(μ)
        
        # Indices adjusted for size of A
        ii, jj = ((imu - 1) * B_dim, (jmu - 1) * B_dim)
            
        # This is equivalent to Z̄ = 1/(1+δ) * C̄m+S̄m = 1/(1+δ) * (A+DAD+AD-DA) 
        # (see eq 11 in Sanghavi et al, 2013)
        for j in 1:B_dim, i in 1:B_dim
            𝐙⁺⁺[ii + i,jj + j] = 2fact * A⁺⁺[i,j,imu,jmu]
            if i <= 2 && j >= 3
                𝐙⁻⁺[ii + i,jj + j] = -2fact * A⁻⁺[i,j,imu,jmu]
            elseif i >= 3 && j <= 2
                𝐙⁻⁺[ii + i,jj + j] = -2fact * A⁻⁺[i,j,imu,jmu]
            else
                𝐙⁻⁺[ii + i,jj + j] = 2fact * A⁻⁺[i,j,imu,jmu]
            end
        end
    end

    # Return Z-moments
    return arr_type(𝐙⁺⁺), arr_type(𝐙⁻⁺)
end

# Tabulated body of compute_Z_moments: identical math —
#     Z(μ, μ′; m) = Σₗ Πₗ(μ) Bₗ Πₗ(μ′)      [Sanghavi 2014, eq. 14-16]
# with P/R/T and Π taken from the per-solve tables instead of being rebuilt,
# and the l-accumulation done on static blocks (no per-(i,j) slice copies).
# Same additions in the same order as the plain method ⇒ bitwise-equal Z.
function _compute_Z_moments_tabulated(mod::AbstractPolarizationType, μ,
                                      greek_coefs::GreekCoefs, m::Int,
                                      arr_type, t::ZMomentTables, Π_pair)
    (; α, β, γ, δ, ϵ, ζ) = greek_coefs
    FT = eltype(β)
    n  = length(μ)

    # Prefactor from the azimuthal expansion: the m = 0 term carries 1/2.
    fact = (m == 0) ? 0.5 : 1.0
    m1 = m + 1                                   # 1-indexed Fourier moment
    l_max = min(length(β), t.l_max)

    # B matrices — the only scatterer-dependent factor (Sanghavi eq. 16).
    𝐁_all = [construct_B_matrix(mod, α, β, γ, δ, ϵ, ζ, l) for l in 1:l_max]
    B_dim = Int(sqrt(length(𝐁_all[1])))

    nb = B_dim * n
    𝐙⁺⁺, 𝐙⁻⁺ = (zeros(FT, nb, nb), zeros(FT, nb, nb))
    A⁺⁺, A⁻⁺ = (zeros(FT, B_dim, B_dim, n, n), zeros(FT, B_dim, B_dim, n, n))

    if Π_pair !== nothing && m1 ≤ l_max
        Π_all, Π⁻_all = Π_pair
        # A(μᵢ, μⱼ) = Σₗ Πₗ(μᵢ) Bₗ Πₗ(±μⱼ)   — Sanghavi eq. 14
        @inbounds for l in m1:l_max
            𝐁 = 𝐁_all[l]
            Π, Π⁻ = Π_all[l], Π⁻_all[l]
            for j in 1:n, i in 1:n
                if B_dim == 1
                    A⁺⁺[1, 1, i, j] += Π[i] * 𝐁 * Π[j]
                    A⁻⁺[1, 1, i, j] += Π[i] * 𝐁 * Π⁻[j]
                else
                    S⁺⁺ = Π[i] * 𝐁 * Π[j]        # static-matrix products,
                    S⁻⁺ = Π[i] * 𝐁 * Π⁻[j]       # no slice allocation
                    for b2 in 1:B_dim, b1 in 1:B_dim
                        A⁺⁺[b1, b2, i, j] += S⁺⁺[b1, b2]
                        A⁻⁺[b1, b2, i, j] += S⁻⁺[b1, b2]
                    end
                end
            end
        end
    end
    # (m ≥ l_max ⇒ the sum is empty and Z ≡ 0, as in the plain method.)

    # Z̄ = 1/(1+δ)·(C̄m + S̄m); the sign flips on the (I,Q)×(U,V) off-diagonal
    # blocks of Z⁻⁺ implement the D-matrix mirror symmetry (Sanghavi eq. 11).
    @inbounds for jmu in 1:n, imu in 1:n
        ii, jj = ((imu - 1) * B_dim, (jmu - 1) * B_dim)
        for j in 1:B_dim, i in 1:B_dim
            𝐙⁺⁺[ii + i, jj + j] = 2fact * A⁺⁺[i, j, imu, jmu]
            if (i <= 2 && j >= 3) || (i >= 3 && j <= 2)
                𝐙⁻⁺[ii + i, jj + j] = -2fact * A⁻⁺[i, j, imu, jmu]
            else
                𝐙⁻⁺[ii + i, jj + j] = 2fact * A⁻⁺[i, j, imu, jmu]
            end
        end
    end

    return arr_type(𝐙⁺⁺), arr_type(𝐙⁻⁺)
end

