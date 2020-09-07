using Revise
using Plots
using Distributions
using RadiativeTransfer
using RadiativeTransfer.PhaseFunction
using KernelAbstractions

# Eqn. 1
function compute_C_scatt(k, nmax, an, bn)
    return (2π/k^2) * sum(n->(2n+1)*(abs2(an[n]) + abs2(bn[n])), 1:nmax)
end

function compute_avg_C_scatt(k, an,bn,w)
    n_ = collect(1:size(an)[1]);
    n_ = 2n_ .+ 1
    return 2π/k^2 * n_' * (w' * (abs2.(an) + abs2.(bn))')'
end

@kernel function compute_Sl_νν!(@Const(wignerA),@Const(wignerB), an, bn, Sνν, @Const(lMax))
    # Get length of an:
    san = size(an)[1];
    # Indices over n and m
    n, m = @index(Global, NTuple)
    
    # Outer loop over l
    for l = 1:lMax
        if max(l-n,n) <= m <= min(san,n+l)
            #println((an[n]' + bn[n]') * (an[m] + bn[m]) * wignerA[l,n,m]^2)
            Sνν[l] += real((an[n]' + bn[n]') * (an[m] + bn[m]) * wignerA[l,n,m]^2)
            Sνν[l] += 1/2 *  real(abs2(an[n] + bn[n]) * wignerA[l,n,m]^2)
        end
    end
    
end

# This can add another flag/number to avoid multiplications by 0 (e.g. where an,bn is 0)
@kernel function avg_anbn!(@Const(an), @Const(bn),mat_anam,mat_bnbm,mat_anbm,mat_bnam,@Const(w),@Const(nMax))
    FT = eltype(an)
    # Indices over n and m
    m, n, i = @index(Global, NTuple)
    if m>=n && m<nMax[i] && n < nMax[i]
        @inbounds mat_anam[n,m] += (w[i] * (an[n,i]' * an[m,i]));
        @inbounds mat_bnbm[n,m] += (w[i] * (bn[n,i]' * bn[m,i]));
        @inbounds mat_anbm[n,m] += (w[i] * (an[n,i]' * bn[m,i]));
        @inbounds mat_bnam[n,m] += (w[i] * (bn[n,i]' * an[m,i]));
 
    end
end

function compute_avg_anbn!(an,bn,w,Nmax,N_max_, mat_anam,mat_bnbm,mat_anbm,mat_bnam)
    FT2 = eltype(an)
    fill!(mat_anam,0)
    fill!(mat_bnbm,0)
    fill!(mat_anbm,0)
    fill!(mat_bnam,0)
    @inbounds for n=1:Nmax
        @inbounds for m=n:Nmax
            anam = FT2(0);
            bnbm = FT2(0);
            anbm = FT2(0);
            bnam = FT2(0);
            @inbounds for i = 1:size(an)[2]
                if m < N_max_[i] && n < N_max_[i]
                    anam += w[i] * an[n,i]' * an[m,i]
                    bnbm += w[i] * bn[n,i]' * bn[m,i]
                    anbm += w[i] * an[n,i]' * bn[m,i]
                    bnam += w[i] * bn[n,i]' * an[m,i]
                end
            end 
            @inbounds mat_anam[n,m] = anam;
            @inbounds mat_bnbm[n,m] = bnbm;
            @inbounds mat_anbm[n,m] = anbm;
            @inbounds mat_bnam[n,m] = bnam;
        end
    end
    return nothing
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
function compute_Sl_νν(l, ν, k, N_max, an,bn, wigner_A, wigner_B,w)
    
    ll = l-1
    coef = (2ll+1)*π / k^2
    #@show ll
    first_term = 0
    second_term = 0

    for n = 1:N_max

        m_star = max(ll-n,n+1)
        m_max = min(ll+n,N_max)
        for m = m_star:m_max
            anam,bnbm,anbm,bnam = get_abnabm(an,bn,n,m,w)
            # real_avg = real(w' *( (adjoint.(ans_bns[1,n,:]) + adjoint.(ans_bns[2,n,:])) .* (ans_bns[1,m,:] + ans_bns[2,m,:])) )
            real_avg = real(anam+anbm+bnam+bnbm);
            #@show real_avg
            #real_avg = 1.3;
            A_lnm2 = 2 * (2m + 1)  * (2n + 1)
            first_term += (real_avg * (wigner_A[m, n, l])^2 * A_lnm2)
        end

        A_lnn2 = (2n + 1)^2
        second_term += w' *(abs2.(an[n,:] .+ bn[n,:])) * (wigner_A[n, n, l])^2 * A_lnn2

    end

    return coef * (first_term + second_term)
end

# Eqn. 22
function compute_Sl_νmν(l, ν, k, N_max, an,bn, wigner_A, wigner_B,w)
    ll = l-1
    coef = (2ll+1)*π / k^2

    first_term = 0
    second_term = 0

    for n = 1:N_max

        m_star  = max(ll-n,n+1)
        m_max = min(ll+n,N_max)
        for m = m_star:m_max
            anam,bnbm,anbm,bnam = get_abnabm(an,bn,n,m,w)
            real_avg = real(anam-anbm-bnam+bnbm);
            #real_avg = real(w' *( (adjoint.(ans_bns[1,n,:]) - adjoint.(ans_bns[2,n,:])) .* (ans_bns[1,m,:] - ans_bns[2,m,:])) )
            #real_avg = 1.2;
            A_lnm2 = ((-1) ^ (ll + n + m)) * 2 * (2m + 1)  * (2n + 1)
            first_term += (real_avg * (wigner_A[m, n, l])^2 * A_lnm2)
        end

        A_lnn2 = ((-1)^ll) * (2n + 1)^2
        second_term += w' *(abs2.(an[n,:] .- bn[n,:])) * (wigner_A[n, n, l])^2 * A_lnn2

    end

    return coef * (first_term + second_term)
    
end

function compute_B2(aerosol::UnivariateAerosol, wigner_A, wigner_B, wl, radius,w)
    #@show w
    # Find overall N_max from the maximum radius
    N_max = PhaseFunction.get_n_max(2 * π * aerosol.r_max/ wl)
    #@show N_max
    # Where to store an, bn, computed over size distribution
    an = zeros(Complex{Float64}, N_max, aerosol.nquad_radius)
    bn = zeros(Complex{Float64}, N_max, aerosol.nquad_radius)
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
                                      view(an, :, i), 
                                      view(bn, :, i), Dn)
    end

    # Compute the average cross-sectional scattering
    k = 2 * π / wl
    avg_C_scatt = compute_avg_C_scatt(k, an,bn,w)

    # Only do these l's for now
    ls = 1:10

    # Where to store the β values
    β_ls = zeros(size(ls, 1))

    # For each l
    for l in ls

        # Compute β_l

        #println(l)
        Sl_00  =  compute_Sl_νν(l, 0, k, N_max, an,bn, wigner_A, wigner_B,w)
        Sl_0m0 = compute_Sl_νmν(l, 0, k, N_max, an,bn, wigner_A, wigner_B,w)
        println("Sl_00: ", Sl_00)
        println("Sl_0m0: ",Sl_0m0)
        println("avg_C_scatt: ",avg_C_scatt)
        β_l = (1/avg_C_scatt) * (Sl_00 + Sl_0m0)
        println("beta: ", β_l)
        β_ls[l] = β_l
    end

    return β_ls

end


function compute_abns(aerosol::UnivariateAerosol,wl,radius)
    FT = eltype(radius)
    # Find overall N_max from the maximum radius
    N_max = PhaseFunction.get_n_max(2 * π * aerosol.r_max/ wl)

    # Where to store an, bn, computed over size distribution
    ans_bns = zeros(Complex{FT}, 2, N_max, aerosol.nquad_radius)

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

function goCPU()
    fill!(mat_anam,0)
    fill!(mat_bnbm,0)
    fill!(mat_bnam,0)
    fill!(mat_anbm,0)
    kernel! = avg_anbn!(CPU())
    event = kernel!(an,bn,mat_anam, mat_bnbm,mat_anbm,mat_bnam,wₓ,N_max_, ndrange=(Nmax,Nmax,length(wₓ))); 
    wait(CPU(), event) 
    @show  mat_anam[10,12]       
end



function goCUDA!()
    fill!(mat_anamC,0);
    fill!(mat_bnbmC,0);
    fill!(mat_bnamC,0);
    fill!(mat_anbmC,0);
    kernel! = avg_anbn!(CUDADevice())
    event = kernel!(anC,bnC,mat_anamC, mat_bnbmC,mat_anbmC,mat_bnamC,wₓC,N_max_C, ndrange=(Nmax,Nmax,length(wₓ))); 
    wait(CUDADevice(), event)    ;
    #@show  Array(mat_anamC)[10,12]  
    #return mat_anamC
end

Nmax = N_max

using LinearAlgebra
N_max_ = PhaseFunction.get_n_max.(2π * r/ wl)

anC = CuArray(an)
bnC = CuArray(bn)
mat_anamC = CuArray(mat_anam)
mat_bnbmC = CuArray(mat_bnbm)
mat_bnamC = CuArray(mat_bnam)
mat_anbmC = CuArray(mat_anbm)
wₓC = CuArray(wₓ)
N_max_C = CuArray(N_max_)
FT2 = Complex{Float64}
mat_anam = UpperTriangular(zeros(FT2,Nmax,Nmax));
mat_bnbm = UpperTriangular(zeros(FT2,Nmax,Nmax));
mat_anbm = UpperTriangular(zeros(FT2,Nmax,Nmax));
mat_bnam = UpperTriangular(zeros(FT2,Nmax,Nmax));

ans_bns = compute_abns(aero, wl, r);
an = ans_bns[1,:,:];
bn = ans_bns[2,:,:];

@time goCPU()

function get_abnabm(an,bn,n,m,w)
    FT = eltype(an)
    anam = FT(0);
    bnbm = FT(0);
    anbm = FT(0);
    bnam = FT(0);
    @inbounds for i=1:size(an)[2]
        anam += w[i] * an[n,i]' * an[m,i]
        bnbm += w[i] * bn[n,i]' * bn[m,i]
        anbm += w[i] * an[n,i]' * bn[m,i]
        bnam += w[i] * bn[n,i]' * an[m,i]
    end
    return anam,bnbm,anbm,bnam
end