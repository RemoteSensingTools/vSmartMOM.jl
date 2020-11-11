# Types for describing atmospheric parameters #
abstract type AbstractAtmParamType end

# Reflected/transmitted stellar irradiation over a given geolocation on a planet 
struct LocalReflectedStarlight{FT,FT2} <: AbstractAtmParamType
    # latitude of ground pixel
    lat::FT
    # longitude of ground pixel
    lon::FT
    # altitude of ground pixel
    alt::FT
    #time of year
    month::FT

    # solar zenith angle
    sza::FT
    # number of viewing angles 
    Ncam::Int 
    # viewing zenith angles
    vza::FT2
    # viewing azimuth angles
    vaz::FT2

    # wavelength in micron
    lambda::FT

    #RT type (redun?)
    # Needs to be adapted later to have a specific Stokes_ type? TBD
    RT_type::AbstractPolarizationType #0: scalar, 1: vector (IQU), 2: vector (IQUV)

    #Surface type (TBD AbstractSurfaceType needs to be defined still)
    Surf_type::AbstractSurfaceType #0: Lambertian, 2: Cox-Munk, 3: RPV
end
# Reflected/transmitted stellar irradiation by the visible disk of a planet 
struct DiskReflectedStarlight{FT, FT2} <: AbstractAtmParamType
    #Geo-quadrature number
    quadlat::Int #Number of quadrature points for latitude
    quadlon::Int #Number of quadrature points for longitude
    # latitude grid of ground pixels
    lat::FT2
    # longitude of ground pixels
    lon::FT2
    # altitude of ground pixels
    alt::FT2 #(unlikely to be used often)
    
    # planetary phase
    ϕ_pl::FT
    # Inclination angle 
    θ_incl::FT
    
    ## <move this to a new struct called RTGeom>
    # local viewing zenith angles
    # vza::FT2
    # local viewing azimuth angles
    # vaz::FT2

    # wavelength in micron
    lambda::FT

    #RT type
    RT_type::Int #0: scalar, 1: vector (IQU), 2: vector (IQUV)

    #Surface type
    Surf_type::Int #0: Lambertian, 2: Cox-Munk, 3: RPV
end
# Stellar/planetary thermal infrared radiation over a given geolocation 
struct LocalTIR{FT,FT2} <: AbstractAtmParamType
    # latitude of ground pixel
    lat::FT
    # longitude of ground pixel
    lon::FT
    # altitude of ground pixel
    alt::FT
    #time of year
    month::FT
    
    # solar zenith angle
    #sza::FT
    # number of viewing angles 
    Ncam::Int 
    # viewing zenith angles
    vza::FT2
    # viewing azimuth angles
    vaz::FT2

    # wavelength in micron
    lambda::FT

    #RT type
    RT_type::Int #0: scalar, 1: vector (IQU), 2: vector (IQUV)

    #Surface type
    Surf_type::Int #0: Lambertian, 2: Cox-Munk, 3: RPV
end
# Stellar/planetary thermal infrared radiation over its visible disk
struct DiskTIR{FT, FT2} <: AbstractAtmParamType
    #Geo-quadrature number
    quadlat::Int #Number of quadrature points for latitude
    quadlon::Int #Number of quadrature points for longitude
    # latitude grid of ground pixels
    lat::FT2
    # longitude of ground pixels
    lon::FT2
    # altitude of ground pixels
    alt::FT2 #(unlikely to be used often)
    
    # Inclination angle 
    θ_incl::FT
    
    ## <move this to a new struct called RTGeom>
    # local viewing zenith angles
    # vza::FT2
    # local viewing azimuth angles
    # vaz::FT2

    # wavelength in micron
    lambda::FT

    #RT type
    RT_type::Int #0: scalar, 1: vector (IQU), 2: vector (IQUV)

    #Surface type
    Surf_type::Int #0: Lambertian, 2: Cox-Munk, 3: RPV
end

abstract type AbstractRTGeomType end
