using Revise
using Plots
using Distributions
using RadiativeTransfer
using RadiativeTransfer.PhaseFunction

# Eqn. 1
function compute_C_scatt(k, nmax, an, bn)
    return (2π/k^2) * sum(n->(2n+1)*(abs2(an[n]) + abs2(bn[n])), 1:nmax)
end

# <|a_n|^2 + |b_n|^2> averaged over size distribution
function compute_avg_C_scatt(k, ans_bns)
    return (2π/k^2) * sum(n->(2n+1)*mean(abs2.(ans_bns[1,n,:]) + abs2.(ans_bns[2,n,:])), 1:size(ans_bns)[2])
end

# Eqn. 22
function compute_Sl_νν(l, ν, k, N_max, ans_bns, wigner_A, wigner_B)
    
    coef = π / k^2

    first_term = 0
    second_term = 0

    for n = 1:N_max

        m_star = max(abs(l-n), n+1)

        for m = m_star:(n+l)

            real_avg = mean( real( (adjoint.(ans_bns[1,n,:]) + adjoint.(ans_bns[2,n,:])) .* (ans_bns[1,m,:] + ans_bns[1,m,:])) )
            A_lnm2 = 2 * (2 * l + 1) * (2 * m + 1)  * (2 * n + 1)
            first_term += (real_avg * (wigner_A[m, n, l])^2 * A_lnm2)
        end

        A_lnn2 = 2 * (2 * l + 1) * (2 * n + 1)^2
        second_term += mean(abs2.(ans_bns[1,n,:] .+ ans_bns[2,n,:])) * (wigner_A[n, n, l])^2 * A_lnn2

    end

    return coef * (first_term + 0.5 * second_term)
end

# Eqn. 22
function compute_Sl_νmν(l, ν, k, N_max, ans_bns, wigner_A, wigner_B)

    coef = π / k^2

    first_term = 0
    second_term = 0

    for n = 1:N_max

        m_star = max(abs(l-n), n+1)

        for m = m_star:(n+l)

            real_avg = mean( real( (adjoint.(ans_bns[1,n,:]) - adjoint.(ans_bns[2,n,:])) .* (ans_bns[1,m,:] - ans_bns[1,m,:])) )
            A_lnm2 = ((-1) ^ (l + n + m)) * 2 * (2 * l + 1) * (2 * m + 1)  * (2 * n + 1)
            first_term += (real_avg * (wigner_A[m, n, l])^2 * A_lnm2)
        end

        A_lnn2 = ((-1)^l) * 2 * (2 * l + 1) * (2 * n + 1)^2
        second_term += mean(abs2.(ans_bns[1,n,:] .- ans_bns[2,n,:])) * (wigner_A[n, n, l])^2 * A_lnn2

    end

    return coef * (first_term + 0.5 * second_term)
    
end

function compute_B(aerosol::UnivariateAerosol, wigner_A, wigner_B)

    # Find overall N_max from the maximum radius
    N_max = PhaseFunction.get_n_max(2 * π * aerosol.r_max/ wl)

    # Where to store an, bn, computed over size distribution
    ans_bns = zeros(Complex, 2, 2 * N_max + 1, aerosol.nquad_radius)

    Dn = zeros(Complex{FT}, N_max)

    # Loop over the size distribution, and compute an, bn, for each size
    for i in 1:aerosol.nquad_radius

        r = i * aerosol.r_max/aerosol.nquad_radius
        size_param = 2 * π * r

        PhaseFunction.compute_mie_ab!(size_param, aerosol.nᵣ - aerosol.nᵢ * im, 
                                      view(ans_bns, 1, :, i), 
                                      view(ans_bns, 2, :, i), Dn)
    end

    # Compute the average cross-sectional scattering
    k = 2 * π / wl
    avg_C_scatt = compute_avg_C_scatt(k, ans_bns)

    # Only do these l's for now
    ls = 1:10

    # Where to store the β values
    β_ls = zeros(size(ls, 1))

    # For each l
    for l in ls

        # Compute β_l

        println(l)
        Sl_00 = compute_Sl_νν(l, 0, k, N_max, ans_bns, wigner_A, wigner_B)
        Sl_0m0 = compute_Sl_νmν(l, 0, k, N_max, ans_bns, wigner_A, wigner_B)
        println("Sl_00: ", Sl_00)
        println("Sl_0m0: ",Sl_0m0)
        println("avg_C_scatt: ",avg_C_scatt)
        β_l = (1/avg_C_scatt) * (Sl_00 + Sl_0m0)
        println("beta: ", β_l)
        β_ls[l] = β_l
    end

    return β_ls

end

# Constants

μ  = 0.3
σ  = 6.82
wl = 0.55

FT = Float64

size_distribution = LogNormal(log(μ), log(σ))

# Generate aerosol:
aero1 = PhaseFunction.UnivariateAerosol(size_distribution, 30.0, 10000, 1.3, 0.0)


N_max = PhaseFunction.get_n_max(2 * π * aero1.r_max/ wl)
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
            wigner_A[m, n, l] = wigner!(m, n, l, -1, 1, 0)
            wigner_B[m, n, l] = wigner!(m, n, l, -1, -1, 2)
        end
    end
end

# Compute B matrix (just betas for now)
β_ls = compute_B(aero1, wigner_A, wigner_B)