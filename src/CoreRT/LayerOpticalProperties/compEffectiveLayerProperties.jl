function constructCoreOpticalProperties(RS_type::AbstractRamanType{FT}, 
                                        iBand, m, model) where FT
    @unpack τ_rayl, τ_aer, τ_abs, aerosol_optics, greek_rayleigh  = model
    #@show typeof(τ_rayl[1]), typeof(τ_aer[1]), typeof(τ_abs[1])
    @assert all(iBand .≤ length(τ_rayl)) "iBand exceeded number of bands"

    arr_type = array_type(model.params.architecture)

    pol_type = model.params.polarization_type
    
    # Quadrature points:
    μ = Array(model.quad_points.qp_μ )
    N = length(model.quad_points.qp_μN)
    # Number of Aerosols:
    nAero = size(τ_aer[iBand[1]],1)
    nZ    = size(τ_rayl[1],2)
    # Rayleigh Z matrix:
    Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, μ, 
                                                    greek_rayleigh, m, 
                                                    arr_type = arr_type);
    
    band_layer_props    = [];
    band_fScattRayleigh = [];
    # @show arr_type
    for iB in iBand
        # Define G here as ones:
        #G = arr_type(ones(FT,N)
        rayl =  [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]),RS_type.ϖ_Cabannes[iB], 
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺) for i=1:nZ]
        
        # Initiate combined properties with rayleigh
        combo = rayl

        # Loop over all aerosol types:
        for i=1:nAero
            # Precomute Z matrices per type (constant per layer)
            #@show typeof(aerosol_optics[iB][i].greek_coefs), typeof(pol_type), typeof(μ)
            AerZ⁺⁺, AerZ⁻⁺ = Scattering.compute_Z_moments(
                                pol_type, μ, 
                                aerosol_optics[iB][i].greek_coefs, 
                                m, arr_type=arr_type)
            #dAerZ⁺⁺ = zeros(4, size(AerZ⁺⁺,1), size(AerZ⁺⁺,2))                    
            #dAerZ⁻⁺ = zeros(4, size(AerZ⁺⁺,1), size(AerZ⁺⁺,2))                    
            
            #for ctr=1:4
            #    dAerZ⁺⁺[ctr,:,:], dAerZ⁻⁺[ctr,:,:] = 
            #                    Scattering.compute_Z_moments(
            #                    pol_type, μ, 
            #                    aerosol_optics[iB][i].greek_coefs, 
            #                    m, arr_type=arr_type)
            #end
            #@show typeof(AerZ⁺⁺), typeof(aerosol_optics[iB][i]), typeof(FT.(τ_aer[iB][i,:]))
            # Generate Core optical properties for Aerosols i
            aer   = createAero.(τ_aer[iB][i,:], 
                                [aerosol_optics[iB][i]], 
                                [AerZ⁺⁺], [AerZ⁻⁺])
            #@show typeof(aer), typeof(combo)
            # Mix with previous Core Optical Properties
            combo = combo .+ aer
        end
        #@show typeof(combo)
        # TODO Type check τ_abs, τ_aer, rayl[i].τ  ./ combo[i].τ
        # Somewhere here we can add canopy later as well!
        ###

        # fScattRayleigh:
        # Assume ϖ of 1 for Rayleight here:
        #@show size(combo)
        #fScattRayleigh = [FT.(Array(rayl[i].τ  ./ combo[i].τ)) for i=1:nZ]
        fScattRayleigh = [Array(rayl[i].τ  ./ combo[i].τ) for i=1:nZ]
        #@show fScattRayleigh, rayl[1].τ, combo[1].τ
        # Create Core Optical Properties merged with trace gas absorptions:
        #@show typeof(combo.+ 
        #[CoreAbsorptionOpticalProperties(arr_type((τ_abs[iB][:,i]))) for i=1:nZ])
        push!(band_layer_props,
                combo .+ 
                [CoreAbsorptionOpticalProperties(arr_type((τ_abs[iB][:,i]))) for i=1:nZ])
        push!(band_fScattRayleigh,fScattRayleigh)
    end

    layer_opt = []
    fscat_opt = []
    for iz = 1:nZ
        push!(layer_opt, prod([band_layer_props[i][iz] for i=1:length(iBand)]));
        push!(fscat_opt, [band_fScattRayleigh[i][iz] for i=1:length(iBand)]);
    end
    # For now just one band_fScattRayleigh
    return layer_opt, fscat_opt # Suniti: this needs to be modified because Rayleigh scattering fraction varies dramatically with wavelength
end

function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺)
    @unpack fᵗ, ω̃ = aerosol_optics
    τ_mod = (1-fᵗ * ω̃ ) * τAer;
    ϖ_mod = (1-fᵗ) * ω̃/(1-fᵗ * ω̃)
    CoreScatteringOpticalProperties(τ_mod, ϖ_mod,AerZ⁺⁺, AerZ⁻⁺)
end

# Extract scattering definitions and integrated absorptions for the source function!
function extractEffectiveProps(
            lods::Array,#{CoreScatteringOpticalProperties{FT},1}
            quad_points::QuadPoints{FT}
            ) where FT

    #FT    = eltype(lods[1].τ)
    nSpec = length(lods[1].τ)
    nZ    = length(lods)
    # First the Scattering Interfaces:
    scattering_interface = ScatteringInterface_00()
    scattering_interfaces_all = []
    τ_sum_all = similar(lods[1].τ,(nSpec,nZ+1))
    τ_sum_all[:,1] .= 0
    #@show FT
    for iz =1:nZ
        # Need to check max entries in Z matrices here as well later!
        scatter = maximum(lods[iz].τ .* lods[iz].ϖ) > 2eps(FT)
        scattering_interface = get_scattering_interface(scattering_interface, scatter, iz)
        push!(scattering_interfaces_all, scattering_interface)
        #@show typeof(τ_sum_all[:,iz]), typeof(lods[iz].τ)
        @views τ_sum_all[:,iz+1] = τ_sum_all[:,iz] + getG_atSun(lods[iz], quad_points) * lods[iz].τ 
    end
    return scattering_interfaces_all, τ_sum_all
end

function getG_atSun(lod::CoreScatteringOpticalProperties,quad_points::QuadPoints{FT}) where FT
    FT(1)
end

function getG_atSun(lod::CoreDirectionalScatteringOpticalProperties,quad_points::QuadPoints{FT}) where FT
    @unpack iμ₀ = quad_points
    gfct = Array(lod.G)[iμ₀]
    return gfct
end


function expandOpticalProperties(in::CoreScatteringOpticalProperties, arr_type)
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = in 
    @assert length(τ) == length(ϖ) "τ and ϖ sizes need to match"
    if size(Z⁺⁺,3) == 1
        Z⁺⁺ = _repeat(Z⁺⁺,1,1,length(τ))
        Z⁻⁺ = _repeat(Z⁻⁺,1,1,length(τ))
        return CoreScatteringOpticalProperties(arr_type(τ), arr_type(ϖ), arr_type(Z⁺⁺), arr_type(Z⁻⁺)) 
    else
        @assert size(Z⁺⁺,3) ==  length(τ) "Z and τ dimensions need to match "
        CoreScatteringOpticalProperties(arr_type(τ), arr_type(ϖ), arr_type(Z⁺⁺), arr_type(Z⁻⁺)) 
    end
end

function expandOpticalProperties(in::CoreDirectionalScatteringOpticalProperties, arr_type)
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺, G = in 
    @assert length(τ) == length(ϖ) "τ and ϖ sizes need to match"
    if size(Z⁺⁺,3) == 1
        Z⁺⁺ = _repeat(Z⁺⁺,1,1,length(τ))
        Z⁻⁺ = _repeat(Z⁻⁺,1,1,length(τ))
        return CoreDirectionalScatteringOpticalProperties(arr_type(τ), arr_type(ϖ), arr_type(Z⁺⁺), arr_type(Z⁻⁺), arr_type(G)) 
    else
        @assert size(Z⁺⁺,3) ==  length(τ) "Z and τ dimensions need to match "
        CoreDirectionalScatteringOpticalProperties(arr_type(τ), arr_type(ϖ), arr_type(Z⁺⁺), arr_type(Z⁻⁺), arr_type(G)) 
    end
end

function expandBandScalars(RS_type, x)
    #test = [length(RS_type.bandSpecLim[iB]) for iB in RS_type.iBand]
    #@show test, sum(test), size(x[1])
    #@show eltype(x[1]),sum([length(RS_type.bandSpecLim[iB]) for iB in RS_type.iBand])
    out = zeros(eltype(x[1]),sum([length(RS_type.bandSpecLim[iB]) for iB in RS_type.iBand]))
    for iB in RS_type.iBand
        out[RS_type.bandSpecLim[iB]] .= expandScalar(x[iB],length(RS_type.bandSpecLim[iB]))
    end
    return out
end

expandScalar(x,n) = x.*ones(n);