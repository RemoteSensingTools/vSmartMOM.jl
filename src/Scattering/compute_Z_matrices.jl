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
"""
    $(FUNCTIONNAME)(mod::AbstractPolarizationType, μ, α, β, γ, δ, ϵ, ζ, m::Int)
Compute moments of the phase matrix
"""
function compute_Z_moments(mod::AbstractPolarizationType, μ, greek_coefs::GreekCoefs, m::Int ; arr_type = Array)
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

