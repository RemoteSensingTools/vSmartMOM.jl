
using NCDatasets
using ProgressMeter

using Distributions
using Interpolations
using Polynomials
using DelimitedFiles




"Struct for an atmospheric profile"
struct AtmosphericProfile{FT}
    lat::FT
    lon::FT
    psurf::FT
    T::Array{FT,1}
    q::Array{FT,1}
    p::Array{FT,1}
    p_levels::Array{FT,1}
    vmr_h2o::Array{FT,1}
    vcd_dry::Array{FT,1}
    vcd_h2o::Array{FT,1}
end;


"Read atmospheric profile (just works for our file, can be generalized"
function read_atmos_profile(file::String, lat::Real, lon::Real, timeIndex; g₀=9.8196)
    @assert 1 <= timeIndex <= 4

    ds = Dataset(file)

    # See how easy it is to actually extract data? Note the [:] in the end reads in ALL the data in one step
    lat_   = ds["YDim"][:]
    lon_   = ds["XDim"][:]
    
    FT = eltype(lat_)
    lat = FT(lat)
    lon = FT(lon)
    
    # Find index (nearest neighbor, one could envision interpolation in space and time!):
    iLat = argmin(abs.(lat_ .- lat))
    iLon = argmin(abs.(lon_ .- lon))

    # Temperature profile
    T    = convert(Array{FT,1}, ds["T"][iLat, iLon, :, timeIndex])
    # specific humidity profile
    q    = convert(Array{FT,1}, ds["QV"][iLat, iLon, :, timeIndex])
    
    # Surafce pressure
    psurf = convert(FT, ds["PS"][iLat, iLon, timeIndex])
    
    # AK and BK global attributes (important to calculate pressure half-levels)
    ak = ds.attrib["HDF_GLOBAL.ak"][:]
    bk = ds.attrib["HDF_GLOBAL.bk"][:]

    p_half = (ak + bk * psurf)
    p_full = (p_half[2:end] + p_half[1:end - 1]) / 2
    close(ds)
    
    # Avogradro's number:
    Na = 6.0221415e23;
    # Dry and wet mass
    dryMass = 28.9647e-3  / Na  # in kg/molec, weighted average for N2 and O2
    wetMass = 18.01528e-3 / Na  # just H2O
    ratio = dryMass / wetMass 
    n_layers = length(T)
    # also get a VMR vector of H2O (volumetric!)
    vmr_h2o = zeros(FT, n_layers, )
    vcd_dry = zeros(FT, n_layers, )
    vcd_h2o = zeros(FT, n_layers, )

    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        vmr_h2o[i] = q[i] * ratio
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dryMass + vmr_h2o[i] * wetMass
        vcd_dry[i] = vmr_dry * Δp / (M * g₀ * 100.0^2)   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * Δp / (M * g₀ * 100^2)
    end
    return AtmosphericProfile(lat, lon, psurf, T, q, p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o)
end;

#for terrestrial atmospheres 
# psurf in hPa, λ in μm 
function getRayleighLayerOptProp(psurf, λ, depol_fct, vcd_dry) 
    # Total vertical Rayleigh scattering optical thickness 
    tau_scat = 0.00864 * (psurf/1013.25) * λ^(-3.916 - 0.074*λ - 0.05/λ) 
    tau_scat = tau_scat*(6.0+3.0*depol_fct)/(6.0-7.0*depol_fct)
    Nz = length(vcd_dry)
    τRayl = zeros(Nz)
    k = tau_scat/sum(vcd_dry)
    for i = 1:Nz
        τRayl[i] = k * vcd_dry[i]
    end
    return τRayl
end

#Gaussian distribution on a pressure grid
function getAerosolLayerOptProp(total_τ, p₀, σp, p_half)
    Nz = length(p_half)
    ρ = zeros(Nz)
    for i = 2:Nz
        dp = p_half[i]-p_half[i-1]
        ρ[i] = (1/(σp*sqrt(2π)))*exp(-(p_half[i]-p₀)^2/(2σp^2))*dp
         #@show ρ[i]  
    end
    Norm = sum(ρ)
    #@show Norm
    τAer  =  (total_τ/Norm) * ρ
    return τAer
end
#computes the composite single scattering parameters (τ, ϖ, Z⁺⁺, Z⁻⁺) for a given atmospheric layer iz for a given Fourier component m
function construct_atm_layer(τRayl, τAer, ϖRayl, ϖAer, fᵗ, Zpp_Rayl, Zmp_Rayl, Zpp_Aer, Zmp_Aer)
    FT = eltype(τRayl)
    @assert length(τAer) == length(ϖAer) == length(fᵗ) "Sizes don't match"
    
    #@show τRayl , sum(τAer)

    τ = FT(0)
    ϖ = FT(0)
    A = FT(0)
    Z⁺⁺ = similar(Zpp_Rayl); 
    Z⁻⁺ = similar(Zpp_Rayl);
    
    if (τRayl + sum(τAer)) < eps(FT)
        fill!(Z⁺⁺,0); fill!(Z⁻⁺,0);
        return FT(0), FT(1), Z⁺⁺, Z⁻⁺
    end
    
    τ += τRayl
    ϖ += τRayl * ϖRayl
    A += τRayl * ϖRayl

    Z⁺⁺ = τRayl * ϖRayl * Zpp_Rayl
    Z⁻⁺ = τRayl * ϖRayl * Zmp_Rayl

    for i = 1:length(τAer)
        τ += τAer[i]
        ϖ += τAer[i] * ϖAer[i]
        A += τAer[i] * ϖAer[i] * (1-fᵗ[i])
        Z⁺⁺ += τAer[i] * ϖAer[i] * (1-fᵗ[i]) * Zpp_Aer[i]
        Z⁻⁺ += τAer[i] * ϖAer[i] * (1-fᵗ[i]) * Zmp_Aer[i]
    end
    
    Z⁺⁺ /= A
    Z⁻⁺ /= A
    A /= ϖ
    ϖ /= τ
    
    #Rescaling composite SSPs according to Eqs. A.3 of Sanghavi et al. (2013) or Eqs.(8) of Sanghavi & Stephens (2015)
    τ *= (FT(1)-(FT(1)-A)*ϖ)
    ϖ *= ϖ*A / (1-(1-A)*ϖ)
    return τ, ϖ, Z⁺⁺, Z⁻⁺  
end


