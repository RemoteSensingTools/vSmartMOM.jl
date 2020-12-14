# Types for describing atmospheric parameters #
abstract type AbstractObsGeometry end

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
end

"Observation Geometry (basics)" 
struct ObsGeometry{FT} <: AbstractObsGeometry
    "altitude of observer `[Pa]`"
    obs_alt::FT
    "solar zenith angle `[Degree]`"
    sza::FT
    "viewing zenith angle(s) `[Degree]`" 
    vza::Array{FT,1}
    "viewing azimuth angle(s) `[Degree]`" 
    vaz::Array{FT,1}
end

"Quadrature Types for RT streams"
abstract type AbstractQuadratureType end

struct RadauQuad <:AbstractQuadratureType end
struct GaussQuadHemisphere <:AbstractQuadratureType end
struct GaussQuadFullSphere <:AbstractQuadratureType end


