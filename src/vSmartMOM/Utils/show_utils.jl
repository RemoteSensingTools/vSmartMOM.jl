
# Overload the show method for vSmartMOM_Parameters
function Base.show(::IO, x::vSmartMOM_Parameters)

    println("----------")
    println("Absorption")
    println("----------")
    println("\tBroadening: $(x.broadening_function)")
    println("\tCEF: $(x.CEF)")
    println("\tWing Cutoff: $(x.wing_cutoff) cm⁻¹")
    bands = map(i-> "$(round(x.spec_grid_start[i], digits=2)):$(x.spec_grid_n[i]):$(round(x.spec_grid_end[i], digits=2))", 1:length(x.spec_grid_start))
    println("\tSpectral Bands (ν₁:n:νₙ, cm⁻¹):")
    for band in bands
        println("\t\t- $(band)")
    end
    println("\tMolecules")
    for molecule_band in x.molecules
        println("\t\t- $(molecule_band)")
    end

    println("\n----------")
    println("Scattering")
    println("----------")
    println("\tAerosols:")
    for i in 1:x.nAer
        println("\t\t- Aerosol #$(i)")
        println("\t\t  τ_ref: $(x.τAer_ref[i])")
        println("\t\t  μ: $(x.μ[i]) μm")
        println("\t\t  σ: $(x.σ[i]) μm")
        println("\t\t  nᵣ: $(x.nᵣ[i])")
        println("\t\t  nᵢ: $(x.nᵢ[i])")
        println("\t\t  p₀: $(x.p₀[i]) Pa") 
        println("\t\t  σp: $(x.σp[i]) Pa")
    end
    println("\tr_max: $(x.r_max) μm")
    println("\tnquad_radius: $(x.CEF)")
    println("\tλ_bands:")
    for λ in x.λ_band
        println("\t\t- $(λ) μm")
    end
    println("\tλ_ref: $(x.λ_ref) μm")
    println("\tDepolarization: $(x.depol)")
    println("\tPolarization Type: $(x.polarization_type)")
    println("\tDecomposition Type: $(x.decomp_type)")

    println("\n--------")
    println("Geometry")
    println("--------")
    println("\tVZA (deg): $(x.vza)")
    println("\tVAZ (deg): $(x.vaz)")
    println("\tSZA (deg): $(x.sza)")
    println("\tObservation Altitude: $(x.obs_alt)")

    println("\n-------------------")
    println("Atmospheric Profile")
    println("-------------------")
    println("\tFile Path: $(x.file)")
    println("\tProfile Reduction: $(x.profile_reduction_n) layers")

    println("\n--------")
    println("Surfaces")
    println("--------")
    println("\tBRDFs: ")
    for surface in x.brdf
        println("\t\t- $(surface)")
    end

    println("\n---------------")
    println("Truncation Type")
    println("---------------")
    println("\tl_trunc: $(x.l_trunc)")
    println("\tΔ_angle: $(x.Δ_angle)°")

    println("\n---------")
    println("vSmartMOM")
    println("---------")
    println("\tquadrature_type: $(x.quadrature_type)")
    println("\tmax_m: $(x.max_m)")
    println("\tarchitecture: $(x.architecture)")
    println("\tfloat_type: $(x.float_type)")
    
end

# Overload the show method for vSmartMOM_Model
function Base.show(::IO, x::vSmartMOM_Model)
    
    println("        params| ::vSmartMOM_Parameters")
    println("aerosol_optics| [iBand][iAer] ($(length(x.aerosol_optics)) x $(length(x.aerosol_optics[1])))")
    println("greek_rayleigh| ::GreekCoefs")
    println("   quad_points| ::QuadPoints")
    println("         τ_abs| [iBand][iSpec,iZ] ($(length(x.τ_abs)) x (nSpec x $(size(x.τ_abs[1])[2])))")
    println("         τRayl| [iBand][iZ] ($(length(x.τRayl)) x $(length(x.τRayl[1])))")
    println("          τAer| [iBand][iAer,iZ] ($(length(x.τAer)) x ($(size(x.τAer[1])[1]) x $(size(x.τAer[1])[2])))")
    println("      obs_geom| ::ObsGeometry")
    println("       profile| ::AtmosphericProfile")
end