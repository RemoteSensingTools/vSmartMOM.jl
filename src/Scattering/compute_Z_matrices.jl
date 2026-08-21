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
    if ndims(β) == 2
        n_spec = size(β, 2)
        moments = map(1:n_spec) do iν
            gc = GreekCoefs(collect(@view(α[:, iν])), collect(@view(β[:, iν])),
                            collect(@view(γ[:, iν])), collect(@view(δ[:, iν])),
                            collect(@view(ϵ[:, iν])), collect(@view(ζ[:, iν])))
            compute_Z_moments(mod, μ, gc, m; arr_type=Array)
        end
        Z⁺⁺ = cat((z[1] for z in moments)...; dims=3)
        Z⁻⁺ = cat((z[2] for z in moments)...; dims=3)
        return arr_type(Z⁺⁺), arr_type(Z⁻⁺)
    elseif ndims(β) != 1
        throw(DimensionMismatch("Greek coefficients must be (l,) or (l,nSpec); got $(size(β))"))
    end
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

"""
    compute_Z_source_moments(mod, μ_out, μ₀, greek_coefs, m; arr_type=Array)

Evaluate the two Fourier phase blocks coupling one collimated incident
direction `μ₀` to the diffuse output directions `μ_out`, without constructing
an augmented square angular matrix. The angular shape is `length(μ_out) × 1`;
after Stokes expansion the returned arrays have shape
`(length(μ_out)*nStokes, nStokes[, nSpec])`.

Keeping the complete incident Stokes block preserves support for polarized
collimated sources. For the usual unpolarized solar beam only column one is
used, so its contraction is an `NquadN × 1` source vector.
"""
function compute_Z_source_moments(mod::AbstractPolarizationType, μ_out,
                                  μ₀::Real, greek_coefs::GreekCoefs, m::Int;
                                  arr_type=Array)
    (; α, β, γ, δ, ϵ, ζ) = greek_coefs
    if ndims(β) == 2
        moments = map(axes(β, 2)) do iν
            gc = GreekCoefs(collect(@view(α[:, iν])), collect(@view(β[:, iν])),
                            collect(@view(γ[:, iν])), collect(@view(δ[:, iν])),
                            collect(@view(ϵ[:, iν])), collect(@view(ζ[:, iν])))
            compute_Z_source_moments(mod, μ_out, μ₀, gc, m; arr_type=Array)
        end
        return arr_type(cat((z[1] for z in moments)...; dims=3)),
               arr_type(cat((z[2] for z in moments)...; dims=3))
    elseif ndims(β) != 1
        throw(DimensionMismatch(
            "Greek coefficients must be (l,) or (l,nSpec); got $(size(β))"))
    end

    FT = promote_type(eltype(β), eltype(μ_out), typeof(μ₀))
    μ = FT.(μ_out)
    μ0 = FT(μ₀)
    @assert all(zero(FT) .< μ .<= one(FT))
    @assert zero(FT) < μ0 <= one(FT)

    fact = m == 0 ? FT(0.5) : one(FT)
    m1 = m + 1
    l_max = length(β)
    # Evaluate the output and solar directions in one recurrence. Besides
    # sharing work, this is bit-identical to the historical augmented-angle
    # evaluation and avoids endpoint-sensitive singleton recurrence behavior.
    μ_eval = vcat(μ, μ0)
    P, R, T = compute_associated_legendre_PRT(μ_eval, l_max)
    Pm, Rm, Tm = compute_associated_legendre_PRT(-μ_eval, l_max)
    B_all = [construct_B_matrix(mod, α, β, γ, δ, ϵ, ζ, l)
             for l in 1:l_max]
    nstokes = Int(sqrt(length(B_all[1])))
    Aplus = zeros(FT, nstokes, nstokes, length(μ))
    Aminus = zeros(FT, nstokes, nstokes, length(μ))

    for l in m1:l_max
        B = B_all[l]
        Π = construct_Π_matrix(mod, P, R, T, l, m1)
        Πm = construct_Π_matrix(mod, Pm, Rm, Tm, l, m1)
        Π0 = Π[end]
        Π0m = Πm[end]
        for iμ in eachindex(μ)
            if nstokes == 1
                Aplus[1, 1, iμ] += Π[iμ] * B * Π0
                Aminus[1, 1, iμ] += Π[iμ] * B * Π0m
            else
                Aplus[:, :, iμ] .+= Π[iμ] * B * Π0
                Aminus[:, :, iμ] .+= Π[iμ] * B * Π0m
            end
        end
    end

    Z0plus = zeros(FT, nstokes * length(μ), nstokes)
    Z0minus = similar(Z0plus)
    for iμ in eachindex(μ), j in 1:nstokes, i in 1:nstokes
        row = (iμ - 1) * nstokes + i
        Z0plus[row, j] = 2fact * Aplus[i, j, iμ]
        cross_qu = (i <= 2 && j >= 3) || (i >= 3 && j <= 2)
        Z0minus[row, j] = (cross_qu ? -2fact : 2fact) * Aminus[i, j, iμ]
    end
    return arr_type(Z0plus), arr_type(Z0minus)
end
