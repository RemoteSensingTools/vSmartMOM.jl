"""
    $(FUNCTIONNAME)(mod::AbstractPolarizationType, μ, α, β, γ, δ, ϵ, ζ, m::Int)
Compute moments of the phase matrix 
"""
function compute_Z_moments(mod::AbstractPolarizationType, μ, greek_coefs::GreekCoefs, m::Int ; arr_type = Array)
    @unpack α, β, γ, δ, ϵ, ζ = greek_coefs
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
    P, R, T    = compute_associated_legendre_PRT( μ, l_max)
    P⁻, R⁻, T⁻ = compute_associated_legendre_PRT(-μ, l_max)
  
    # Pre-compute all required B matrices
    𝐁_all = [construct_B_matrix(mod, α, β, γ, δ, ϵ, ζ, i) for i in 1:l_max]
#@show 𝐁_all 
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
        #@show Π, Π⁻
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
        #=for i in eachindex(μ)
                @show m,l
                @show A⁺⁺[:,:,i,i] 
                @show A⁻⁺[:,:,i,i] 
        end=#
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
    $(FUNCTIONNAME)(mod::AbstractPolarizationType, μ,  μ₀, greek_coefs::GreekCoefs, m::Int ; arr_type = Array)
Compute moments of the phase matrix 
"""
function compute_Z_moments(mod::AbstractPolarizationType, μ, μ₀, greek_coefs::GreekCoefs, m::Int ; arr_type = Array)
    @unpack α, β, γ, δ, ϵ, ζ = greek_coefs
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
    P, R, T    = Scattering.compute_associated_legendre_PRT(μ, l_max)
    P⁻, R⁻, T⁻ = Scattering.compute_associated_legendre_PRT(-μ, l_max)
    
    # For incoming direction only:
    μ₀P, μ₀R, μ₀T    = Scattering.compute_associated_legendre_PRT(μ₀, l_max)
    μ₀P⁻, μ₀R⁻, μ₀T⁻ = Scattering.compute_associated_legendre_PRT(-μ₀, l_max)
  
    # Pre-compute all required B matrices
    𝐁_all = [construct_B_matrix(mod, α, β, γ, δ, ϵ, ζ, i) for i in 1:l_max]

    # Get dimension of square matrix (easier for Scalar/Stokes dimensions)
    B_dim = Int(sqrt(length(𝐁_all[1])))
    
    # Create matrices:
    nb = B_dim * n
    𝐙⁺⁺, 𝐙⁻⁺ = (zeros(FT, nb, nb), zeros(FT, nb, nb));
    A⁺⁺, A⁻⁺ = (zeros(FT, B_dim, B_dim, n, n), zeros(FT, B_dim, B_dim, n, n));

    μ₀𝐙⁺⁺, μ₀𝐙⁻⁺ = (zeros(FT, nb), zeros(FT, nb));
    μ₀A⁺⁺, μ₀A⁻⁺ = (zeros(FT, B_dim, B_dim, n), zeros(FT, B_dim, B_dim, n));
    # Iterate over l
    for l = m:l_max

        # B matrix for l
        𝐁 = 𝐁_all[l];

        # Construct Π matrix for l,m pair (change to in place later!)
        # See eq. 15 in Sanghavi 2014, note that P,R,T are already normalized
        Π    = construct_Π_matrix(mod, P, R, T, l, m)
        Π⁻   = construct_Π_matrix(mod, P⁻, R⁻, T⁻, l, m)
        #μ₀Π  = construct_Π_matrix(mod, μ₀P, μ₀R, μ₀T, l, m)[1]
        #μ₀Π⁻ = construct_Π_matrix(mod, μ₀P⁻, μ₀R⁻, μ₀T⁻, l, m)
        
        i = 1; j=3
        
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
        #=for j in eachindex(μ)
            if B_dim == 1
                #@show μ₀Π * 𝐁 * Π[j]
                μ₀A⁺⁺[B_dim,B_dim,i] += μ₀Π * 𝐁 * Π[j]
                μ₀A⁻⁺[B_dim,B_dim,i] += μ₀Π * 𝐁 * Π⁻[j]  
            else
                #@show size((μ₀Π * 𝐁 * Π[j])*mod.I₀)
                μ₀A⁺⁺[:,:,j] += μ₀Π * 𝐁 * Π[j]
                μ₀A⁻⁺[:,:,j] += μ₀Π * 𝐁 * Π⁻[j]
            end
        end=#
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

    # for μ₀
    # Now get to the Z part:
    for jmu in eachindex(μ)
        imu = 1
        # Indices adjusted for size of A
        ii, jj = ((imu - 1) * B_dim, (jmu - 1) * B_dim)
            
        # This is equivalent to Z̄ = 1/(1+δ) * C̄m+S̄m = 1/(1+δ) * (A+DAD+AD-DA) 
        # (see eq 11 in Sanghavi et al, 2013)
        #=@inbounds for j in 1:B_dim
            @show size((2fact * μ₀A⁺⁺[:,j,jmu])' * mod.I₀)
            μ₀𝐙⁺⁺[jj + j] = (2fact * μ₀A⁺⁺[:,j,jmu])' * mod.I₀
            if i <= 2 && j >= 3
                μ₀𝐙⁻⁺[jj + j] = (-2fact * μ₀A⁻⁺[:,j,jmu])' * mod.I₀
            elseif i >= 3 && j <= 2
                μ₀𝐙⁻⁺[jj + j] = (-2fact * μ₀A⁻⁺[:,j,jmu])' * mod.I₀
            else
                μ₀𝐙⁻⁺[jj + j] = (2fact * μ₀A⁻⁺[:,j,jmu])' * mod.I₀
            end
        end=#
    end


    # Return Z-moments
    return arr_type(𝐙⁺⁺), arr_type(𝐙⁻⁺)
end