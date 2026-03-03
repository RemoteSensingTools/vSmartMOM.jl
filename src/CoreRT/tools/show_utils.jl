#=
 
This file specifies how to pretty-print vSmartMOM_Parameters and RTModel types
 
=#

function Base.show(io::IO,::MIME"text/plain", x::vSmartMOM_Parameters)

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
        for i in 1:length(x.scattering_params.rt_aerosols)
            println(io, "\t\t- Aerosol #$(i)")
            curr_rt_aerosol = x.scattering_params.rt_aerosols[i]
            println(io, "\t\t  τ_ref: $(curr_rt_aerosol.τ_ref)")
            println(io, "\t\t  Size Distribution: $(curr_rt_aerosol.aerosol.size_distribution)")
            #println(io, "\t\t  σ: $(curr_rt_aerosol.aerosol.size_distribution.σ) μm; geometric standard deviation")
            println(io, "\t\t  nᵣ: $(curr_rt_aerosol.aerosol.nᵣ)")
            println(io, "\t\t  nᵢ: $(curr_rt_aerosol.aerosol.nᵢ)")
            
            println(io, "\t\t  Vertical profile in p: $(curr_rt_aerosol.profile)") 
        end
        println(io, "\tr_max: $(x.scattering_params.r_max) μm")
        println(io, "\tnquad_radius: $(x.scattering_params.nquad_radius)")
        println(io, "\tλ_ref: $(x.scattering_params.λ_ref) μm")
        println(io, "\treference n: $(x.scattering_params.n_ref)")
        print(io, "\tDecomposition Type: $(x.scattering_params.decomp_type)")
    else
        print(io, "(No Scattering Parameters Specified)")
    end
end

# Pretty-print RTModel
function Base.show(io::IO, x::RTModel)
    nBands = length(x.atmosphere.spec_bands)
    nAer = n_aerosols(x)
    println(io, "RTModel{$(typeof(x.architecture)), $(float_type(x))}")
    println(io, "  architecture| $(x.architecture)")
    println(io, "        solver| $(typeof(x.solver.polarization_type)), Nquad=$(x.quad_points.Nquad), max_m=$(x.solver.max_m)")
    println(io, "      geometry| SZA=$(x.geometry.sza)°, VZA=$(x.geometry.vza)°")
    println(io, "    atmosphere| $(nBands) band(s), $(length(x.atmosphere.profile.p_full)) levels")
    println(io, "       optics | τ_abs: $(nBands) band(s), $(nAer) aerosol(s)")
    println(io, "     surfaces | $(nBands) BRDF(s): $(join([typeof(s).name.name for s in x.surfaces], ", "))")
end