
# Overload the show method for vSmartMOM_Parameters
function Base.show(io::IO, x::vSmartMOM_Parameters)

    println(io, "\n------------------")
    println(io, "Radiative Transfer")
    println(io, "------------------")
    println(io, "\tSpectral Bands (ν₁:n:νₙ, cm⁻¹):")
    for band in x.spec_bands
        println(io, "\t\t- $(length(band))-length vector from $(round(band[1], digits=2)) to $(round(band[end], digits=2))")
    end
    println(io, "\tSurface BRDFs: ")
    for surface in x.brdf
        println(io, "\t\t- $(surface)")
    end
    println(io, "\tQuadrature Type: $(x.quadrature_type)")
    println(io, "\tPolarization Type: $(x.polarization_type)")
    println(io, "\tmax_m: $(x.max_m)")
    println(io, "\tΔ_angle: $(x.Δ_angle)°")
    println(io, "\tl_trunc: $(x.l_trunc)")
    println(io, "\tDepolarization: $(x.depol)")
    println(io, "\tFloating-point type: $(x.float_type)")
    println(io, "\tArchitecture: $(x.architecture)")

    println(io, "\n--------")
    println(io, "Geometry")
    println(io, "--------")
    println(io, "\tSZA (deg): $(x.sza)")
    println(io, "\tVZA (deg): $(x.vza)")
    println(io, "\tVAZ (deg): $(x.vaz)")
    println(io, "\tObservation Altitude: $(x.obs_alt)")

    println(io, "\n-------------------")
    println(io, "Atmospheric Profile")
    println(io, "-------------------")
    println(io, "\tTemperature: $(length(x.T))-length array")
    println(io, "\tPressure: $(length(x.p))-length array")
    println(io, "\tSpecific humidity: $(length(x.q))-length array")
    if x.profile_reduction_n == -1
        println(io, "\tProfile Reduction: None")
    else
        println(io, "\tProfile Reduction: $(x.profile_reduction_n) layers")
    end

    println(io, "\n----------")
    println(io, "Absorption")
    println(io, "----------")
    if (!isnothing(x.absorption_params))
        println(io, "\tMolecules")
        for molecule_band in x.absorption_params.molecules
            println(io, "\t\t- $(molecule_band)")
        end
        println(io, "\tVMR:")
        for mol in keys(x.absorption_params.vmr)
            println(io, "\t\t- $(mol): $(x.absorption_params.vmr[mol])")
        end
        println(io, "\tBroadening: $(x.absorption_params.broadening_function)")
        println(io, "\tCEF: $(x.absorption_params.CEF)")
        println(io, "\tWing Cutoff: $(x.absorption_params.wing_cutoff) cm⁻¹")
    else
        println(io, "(No Absorption Parameters Specified)")
    end
    
    println(io, "\n----------")
    println(io, "Scattering")
    println(io, "----------")
    if (!isnothing(x.scattering_params))
        println(io, "\tAerosols:")
        for i in 1:length(x.scattering_params.aerosols)
            println(io, "\t\t- Aerosol #$(i)")
            println(io, "\t\t  τ_ref: $(x.scattering_params.aerosols[i].τ_ref)")
            println(io, "\t\t  μ: $(x.scattering_params.aerosols[i].μ) μm")
            println(io, "\t\t  σ: $(x.scattering_params.aerosols[i].σ) μm")
            println(io, "\t\t  nᵣ: $(x.scattering_params.aerosols[i].nᵣ)")
            println(io, "\t\t  nᵢ: $(x.scattering_params.aerosols[i].nᵢ)")
            println(io, "\t\t  p₀: $(x.scattering_params.aerosols[i].p₀) Pa") 
            println(io, "\t\t  σp: $(x.scattering_params.aerosols[i].σp) Pa")
        end
        println(io, "\tr_max: $(x.scattering_params.r_max) μm")
        println(io, "\tnquad_radius: $(x.scattering_params.nquad_radius)")
        println(io, "\tλ_ref: $(x.scattering_params.λ_ref) μm")
        print(io, "\tDecomposition Type: $(x.scattering_params.decomp_type)")
    else
        print(io, "(No Scattering Parameters Specified)")
    end
    
end

# Overload the show method for vSmartMOM_Model
function Base.show(io::IO, x::vSmartMOM_Model)
    
    println(io, "        params| ::vSmartMOM_Parameters")
    println(io, "aerosol_optics| [i_band][i_aer] ($(length(x.aerosol_optics)) x $(length(x.aerosol_optics[1])))")
    println(io, "greek_rayleigh| ::GreekCoefs")
    println(io, "   quad_points| ::QuadPoints")
    println(io, "         τ_abs| [i_band][i_spec,i_z] ($(length(x.τ_abs)) x (n_spec x $(size(x.τ_abs[1])[2])))")
    println(io, "        τ_rayl| [i_band][i_z] ($(length(x.τ_rayl)) x $(length(x.τ_rayl[1])))")
    println(io, "         τ_aer| [i_band][i_aer,i_z] ($(length(x.τ_aer)) x ($(size(x.τ_aer[1])[1]) x $(size(x.τ_aer[1])[2])))")
    println(io, "      obs_geom| ::ObsGeometry")
    println(io, "       profile| ::AtmosphericProfile")
end