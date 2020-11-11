#= Main file for RTM =#
using FastGaussQuadrature   
function rt_set_streams(Ltrunc, quadtype, sza, vza)
    #Ltrunc + 1 = number of spherical coefficients considered (+1 for l=0)
    #quadtype = 'r' for Radau (with DNI), 'g' for Gauss (no DNI)
    #sza = single solar zenith angle
    #lza = array of or single line-of-sight zenith angle
    Nquad=ceil(0.5*(Ltrunc+1))
    if quadtype='g'
        qp_μ, wt_μ = gausslegendre(Nquad, μ1, μ2) #quadrature limits are 0.0-1.0
        return Nquad, qp_μ, wt_μ
    elseif quadtype='r'
        tqp_μ0, twt_μ0 = gaussradau(Nquad)
        qp_μ0 = -1*reverse(tqp_μ0)
        wt_μ0 = reverse(twt_μ0)
        μ0 = cos(sza) #check for degree to radian conversion
        for i in 1:Nquad
            qp_μ[i] = 0.5 * μ0 * qp_μ0[i]
            wt_μ[i] = 0.5 * μ0 * wt_μ0[i]
            qp_μ[Nquad+i] = 0.5 * ((1 + μ0) + (1 - μ0) * qp_μ0[i])
            wt_μ[Nquad+i] = 0.5 * (1-μ0) * wt_μ0[i]
        end
        Ncam = length(vza)
        μ = cos(vza) #check for degree to radian conversion

        # Screen out duplicate camera zenith angles
        txμ = zeros(Ncam)
        idx_vza = zeros(Ncam)
        txμ = 99.99
        txμ[1] = μ[1]
        xctr=1
        for icam in 2:Ncam
            ix=1
            while ((ix<=xctr)&(abs(txμ[ix]-μ[icam])>=1.e-8))
                ix = ix+1
            end
            if (ix>xctr)
                xctr = xctr+1 #Number of unique angles
                txμ[xctr] = μ[icam]
            end
        end 
        xμ = sort(txμ[1:xctr]) #Sorting xμ in increasing order

        # Find camera zenith cosines that are identical to the quadrature cosines   
        insert = zeros(xctr)
        for icam in 1:xctr
            for iμ = 1:2*Nquad
                if abs(qp_μ[iμ] - xμ[icam] > 1.e-8) # distinct cosines
                    insert[icam] = 1
                else #camera cosine already exists
                    insert[icam] = 0
                end
            end
        end
        # Attach unique viewing to quadrature rule
        iμ = 2*Nquad+1
        for icam in 1:xctr
            if (insert[icam]==1)
                qp_μ[iμ] = xμ[icam]
                wt_μ[iμ] = 0.0
                iμ = iμ + 1
            end
        end
        Nquad = Nquad*2 + sum(insert)
        return Nquad, qp_μ, wt_μ
    end   
end

# Finds index i of f_array (i) which is nearest point to f
function nearest_point(Npts, f_array, f)
    d0 = 999.9
    for i in 1:Npts
        d = abs(f_array[i]-f)
        if d<d0
            d0=d
            index=i
        end     
    end
    return index
end

# Assuming a vertical array of SSPs of molecular 
# (Rayleigh) scattering and <naer> aerosol species, compute
# composite single scattering properties τ_comp[1:nz], 
# ϖ_comp[1:nz] and Z matrices      
struct params
    Nz::Int #Number of vertical layers

end
struct atmos

end
