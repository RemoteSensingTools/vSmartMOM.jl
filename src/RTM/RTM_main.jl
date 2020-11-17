#= Main file for RTM =#


function rt_set_streams(model::GaussQuadHemisphere, Ltrunc::Int, sza, vza)
    Nquad=convert(Int,ceil((Ltrunc+1)/2))
    qp_μ, wt_μ = gauleg(Nquad, 0.0, 1.0) #quadrature limits are 0.0-1.0
    return Nquad, qp_μ, wt_μ
end

function rt_set_streams(model::GaussQuadFullSphere, Ltrunc::Int, sza, vza)
    Nquad=convert(Int,ceil((Ltrunc+1)/2))
    qp_μ, wt_μ = gausslegendre(2Nquad) #quadrature limits are 0.0-1.0
    return Nquad, qp_μ[Nquad+1:end], wt_μ[Nquad+1:end]
end

function rt_set_streams(model::RadauQuad, Ltrunc::Int, sza, vza)
    FT = eltype(sza)
    #Ltrunc + 1 = number of spherical coefficients considered (+1 for l=0)
    #quadtype = 'r' for Radau (with DNI), 'g' for Gauss (no DNI)
    #sza = single solar zenith angle
    #lza = array of or single line-of-sight zenith angle
    Nquad=convert(Int,ceil((Ltrunc+1)/2))
    
    tqp_μ₀, twt_μ₀ = gaussradau(Nquad)
    qp_μ₀ = -reverse(tqp_μ₀)
    wt_μ₀ =  reverse(twt_μ₀)
    μ₀ = cosd(sza) #check for degree to radian conversion

    qp_μ = zeros(FT,2Nquad)
    wt_μ = zeros(FT,2Nquad)
    for i = 1:Nquad
        qp_μ[i] = (μ₀ + μ₀ * qp_μ₀[i])/2
        wt_μ[i] = μ₀ * wt_μ₀[i]/2
        qp_μ[Nquad+i] = ((1+μ₀) + (1-μ₀) * qp_μ₀[i])/2
        wt_μ[Nquad+i] = (1-μ₀) * wt_μ₀[i]/2
    end
    Ncam = length(vza)
    μ = cosd.(vza) #check for degree to radian conversion
    
    # Screen out duplicate camera zenith angles
    qp_μ = unique([qp_μ; cosd.(vza)])
    # Assign zero-weights to remaining camera zenith angles
    wt_μ = [wt_μ; zeros(length(qp_μ)-length(wt_μ))]
    
    Nquad = length(qp_μ)
    return Nquad, qp_μ, wt_μ   
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
#struct params
#    Nz::Int #Number of vertical layers#
#
#end
#struct atmos

#end
