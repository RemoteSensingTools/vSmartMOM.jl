using Revise
using Plots
using Distributions
using RadiativeTransfer
using RadiativeTransfer.PhaseFunction

# Eqn. 1
function compute_C_scatt(k, nmax, an, bn)
    return (2π/k^2) * sum(n->(2n+1)*(abs2(an[n]) + abs2(bn[n])), 1:nmax)
end

function compute_avg_C_scatt(k, an,bn,w)
    n_ = collect(1:size(an)[1]);
    n_ = 2n_ .+ 1
    return 2π/k^2 * n_' * (w' * (abs2.(an) + abs2.(bn))')'
end

function compute_avg_anbn(an,bn,w,Nmax)
    FT = eltype(an)
    mat_anam = UpperTriangular(zeros(FT,Nmax,Nmax))
    mat_bnbm = UpperTriangular(zeros(FT,Nmax,Nmax))
    mat_anbm = UpperTriangular(zeros(FT,Nmax,Nmax))
    mat_bnam = UpperTriangular(zeros(FT,Nmax,Nmax))
    @inbounds for n=1:Nmax
        @inbounds for m=n:Nmax
            #anam = 0;
            #bnbm = 0;
            #anbm = 0;
            #bnam = 0;
            @inbounds for i = 1:size(an)[2]
                #anam += w[i] * an[n,i]' * an[m,i]
                #bnbm += w[i] * bn[n,i]' * bn[m,i]
                #anbm += w[i] * an[n,i]' * bn[m,i]
                #bnam += w[i] * bn[n,i]' * an[m,i]
                mat_anam[n,m] += w[i] * an[n,i]' * an[m,i]
                mat_bnbm[n,m] += w[i] * bn[n,i]' * bn[m,i]
                mat_anbm[n,m] += w[i] * an[n,i]' * bn[m,i]
                mat_bnam[n,m] += w[i] * bn[n,i]' * an[m,i]
            end 
            #@inbounds mat_anam[n,m] = anam
        end
    end
    return mat_anam, mat_bnam, mat_anbm, mat_bnam
end

# <|a_n|^2 + |b_n|^2> averaged over size distribution
function compute_avg_C_scatt2(k, ans_bns,w)
    a = 0
    for n=1:size(ans_bns)[2]
        for m=1:size(ans_bns)[3]
            a += (2n+1)* (w[m] * (abs2(ans_bns[1,n,m]) + abs2(ans_bns[2,n,m])))
        end
    end
    return (2π/k^2) * a
end

# Eqn. 22
function compute_Sl_νν(l, ν, k, N_max, ans_bns, wigner_A, wigner_B,w)
    ll = l-1
    coef = (2ll+1)*π / k^2
    #@show ll
    first_term = 0
    second_term = 0

    for n = 1:N_max

        m_star = max(ll-n,n+1)
        m_max = min(ll+n,N_max)
        for m = m_star:m_max

            real_avg = real(w' *( (adjoint.(ans_bns[1,n,:]) + adjoint.(ans_bns[2,n,:])) .* (ans_bns[1,m,:] + ans_bns[2,m,:])) )
            #@show real_avg
            #real_avg = 1.3;
            A_lnm2 = 2 * (2m + 1)  * (2n + 1)
            first_term += (real_avg * (wigner_A[m, n, l])^2 * A_lnm2)
        end

        A_lnn2 = (2n + 1)^2
        second_term += w' *(abs2.(ans_bns[1,n,:] .+ ans_bns[2,n,:])) * (wigner_A[n, n, l])^2 * A_lnn2

    end

    return coef * (first_term + second_term)
end

# Eqn. 22
function compute_Sl_νmν(l, ν, k, N_max, ans_bns, wigner_A, wigner_B,w)
    ll = l-1
    coef = (2ll+1)*π / k^2

    first_term = 0
    second_term = 0

    for n = 1:N_max

        m_star  = max(ll-n,n+1)
        m_max = min(ll+n,N_max)
        for m = m_star:m_max

            real_avg = real(w' *( (adjoint.(ans_bns[1,n,:]) - adjoint.(ans_bns[2,n,:])) .* (ans_bns[1,m,:] - ans_bns[2,m,:])) )
            #real_avg = 1.2;
            A_lnm2 = ((-1) ^ (ll + n + m)) * 2 * (2m + 1)  * (2n + 1)
            first_term += (real_avg * (wigner_A[m, n, l])^2 * A_lnm2)
        end

        A_lnn2 = ((-1)^ll) * (2n + 1)^2
        second_term += w' *(abs2.(ans_bns[1,n,:] .- ans_bns[2,n,:])) * (wigner_A[n, n, l])^2 * A_lnn2

    end

    return coef * (first_term + second_term)
    
end

function compute_B2(aerosol::UnivariateAerosol, wigner_A, wigner_B, wl, radius,w)
    #@show w
    # Find overall N_max from the maximum radius
    N_max = PhaseFunction.get_n_max(2 * π * aerosol.r_max/ wl)
    #@show N_max
    # Where to store an, bn, computed over size distribution
    ans_bns = zeros(Complex{Float64}, 2, 2 * N_max + 1, aerosol.nquad_radius)

    #Dn = zeros(Complex{FT}, N_max)

    # Loop over the size distribution, and compute an, bn, for each size
    for i in 1:aerosol.nquad_radius

        r = radius[i] 
        size_param = 2 * π * r / wl
        # Pre-allocate Dn:
        y = size_param * (aerosol.nᵣ-aerosol.nᵢ);
        nmx = round(Int, max(N_max, abs(y))+51 )
        Dn = zeros(Complex{FT},nmx)

        PhaseFunction.compute_mie_ab!(size_param, aerosol.nᵣ + aerosol.nᵢ * im, 
                                      view(ans_bns, 1, :, i), 
                                      view(ans_bns, 2, :, i), Dn)
    end

    # Compute the average cross-sectional scattering
    k = 2 * π / wl
    avg_C_scatt = compute_avg_C_scatt(k, ans_bns,w)

    # Only do these l's for now
    ls = 1:10

    # Where to store the β values
    β_ls = zeros(size(ls, 1))

    # For each l
    for l in ls

        # Compute β_l

        #println(l)
        Sl_00  =  compute_Sl_νν(l, 0, k, N_max, ans_bns, wigner_A, wigner_B,w)
        Sl_0m0 = compute_Sl_νmν(l, 0, k, N_max, ans_bns, wigner_A, wigner_B,w)
        println("Sl_00: ", Sl_00)
        println("Sl_0m0: ",Sl_0m0)
        println("avg_C_scatt: ",avg_C_scatt)
        β_l = (1/avg_C_scatt) * (Sl_00 + Sl_0m0)
        println("beta: ", β_l)
        β_ls[l] = β_l
    end

    return β_ls

end


function compute_abns(aerosol::UnivariateAerosol, wigner_A, wigner_B,wl,radius)

    # Find overall N_max from the maximum radius
    N_max = PhaseFunction.get_n_max(2 * π * aerosol.r_max/ wl)

    # Where to store an, bn, computed over size distribution
    ans_bns = zeros(Complex{Float32}, 2, N_max, aerosol.nquad_radius)

    #Dn = zeros(Complex{FT}, N_max)

    # Loop over the size distribution, and compute an, bn, for each size
    for i in 1:aerosol.nquad_radius

        r = radius[i]
        size_param = 2 * π * r / wl
        n_max = PhaseFunction.get_n_max(size_param)
        y = size_param * (aerosol.nᵣ-aerosol.nᵢ);
        nmx = round(Int, max(n_max, abs(y))+51 )
        Dn = zeros(Complex{FT},nmx)
        
        PhaseFunction.compute_mie_ab!(size_param, aerosol.nᵣ + aerosol.nᵢ * im, 
                                      view(ans_bns, 1, 1:n_max, i), 
                                      view(ans_bns, 2, 1:n_max, i), Dn)
    end

    # Compute the average cross-sectional scattering
    k = 2 * π / wl
    #avg_C_scatt = compute_avg_C_scatt(k, ans_bns)
    return ans_bns;
    

end


# Constants

const μ  = 0.3
const σ  = 6.82
const wl = 0.55

const FT = Float64

size_distribution = LogNormal(log(μ), log(σ))

# Generate aerosol:
aero = PhaseFunction.UnivariateAerosol(size_distribution, 30.0, 2500, 1.3, 0.0)
r, wᵣ = PhaseFunction.gauleg(aero.nquad_radius, 0.0, aero.r_max ; norm=true)
wₓ = pdf.(aero.size_distribution,r)
# pre multiply with wᵣ to get proper means eventually:
wₓ .*= wᵣ
# normalize (could apply a check whether cdf.(aero.size_distribution,r_max) is larger than 0.99:
wₓ /= sum(wₓ)

N_max = PhaseFunction.get_n_max(2 * π * aero.r_max/ wl)
ls = 1:10
ms = 1:(2 * N_max + 1)
ns = 1:(N_max + 1)

# Pre-compute Wigner symbols

# Note: these only have the wigner symbols, not the added coefficients like 2l+1, 2m+1, etc. 

wigner_A = zeros(size(ms, 1), size(ns, 1), size(ls, 1))
wigner_B = zeros(size(ms, 1), size(ns, 1), size(ls, 1))

for l in ls
    for m in ms
        for n in ns 
            wigner_A[m, n, l] = wigner!(m, n, l-1, -1, 1, 0)
            wigner_B[m, n, l] = wigner!(m, n, l-1, -1, -1, 2)
        end
    end
end

# Compute B matrix (just betas for now)
β_ls = compute_B2(aero, wigner_A, wigner_B,wl,r,wₓ)