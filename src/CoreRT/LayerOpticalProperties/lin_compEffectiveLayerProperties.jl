# Construct τ, ϖ, Z⁺⁺, Z⁻⁺ and
# lin_τ, lin_ϖ, linZ⁺⁺, linZ⁻⁺ wrt all state vector parameters
# for each atmospheric layer

function constructCoreOpticalProperties(
        RS_type::AbstractRamanType{FT}, 
        iBand, m, model, lin) where FT
    @unpack τ_rayl, τ_aer, τ_abs, aerosol_optics, greek_rayleigh  = model
    @unpack lin_τ_abs, lin_τ_rayl, lin_τ_aer_psurf, lin_τ_aer_z₀, lin_τ_aer_σz, lin_aerosol_optics = lin

    #@show typeof(τ_rayl[1]), typeof(τ_aer[1]), typeof(τ_abs[1])
    @assert all(iBand .≤ length(τ_rayl)) "iBand exceeded number of bands"

    arr_type = array_type(model.params.architecture)

    pol_type = model.params.polarization_type
    #nparams = 7*Naer + 6 + 2 + 1 + 2*n_bands + 2 
    #  7 for each aerosol type ,
    #  6 for CO2 VMR (VMR(z)=a₁f₁(z)+a₂f₂(z)+a₃f₃(z)),
    #  2 for H2O VMR (VMR(z)=b₁f₁(z)+b₂f₂(z)),
    #  1 for p_surf,
    #  nbands for ρ₀ (RPV),
    #  nbands for ρc (RPV),
    #  1 for Θ (RPV),
    #  1 for k (RPV) 
    nparams = 23 #excluding surface for now
    # Quadrature points:
    μ = Array(model.quad_points.qp_μ )
    N = length(model.quad_points.qp_μN)
    # Number of Aerosols:
    nAero = size(τ_aer[iBand[1]],1)
    Nz    = size(τ_rayl[1],2)
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
        #rayl =  [CoreScatteringOpticalProperties(
        #    arr_type(τ_rayl[iB][:,i]),RS_type.ϖ_Cabannes[iB], 
        #    Rayl𝐙⁺⁺, Rayl𝐙⁻⁺) for i=1:Nz]

        #@show size(τ_rayl[iB][:,end]), size(rayl[end].τ) 
        #@show RS_type.ϖ_Cabannes[iB], rayl[end].ϖ
        #@show size(Rayl𝐙⁺⁺), size(rayl[end].Z⁺⁺)     
        dtmp_τ = [zeros(FT, nparams, size(lin_τ_rayl[1],1)) for i=1:Nz]
        for iz = 1:Nz
            @show size(lin_τ_rayl[iB],1)
            dtmp_τ[iz][1,:] .= lin_τ_rayl[iB][:,iz]
        end
        dtmp_ϖ = zeros(FT, nparams)
        dtmp_Z⁺⁺ = zeros(FT, nparams, size(Rayl𝐙⁺⁺,1), size(Rayl𝐙⁺⁺,1))
        dtmp_Z⁻⁺ = zeros(FT, nparams, size(Rayl𝐙⁺⁺,1), size(Rayl𝐙⁺⁺,1))


        #lin_rayl =  [linCoreScatteringOpticalProperties( 
        #        arr_type(dtmp_τ[iz]), 
        #        dtmp_ϖ, 
        #        dtmp_Z⁺⁺, dtmp_Z⁻⁺) 
        #        for iz=1:Nz]
        rayl =  [scatt_paired_xdx(CoreScatteringOpticalProperties(
                            arr_type(τ_rayl[iB][:,i]),RS_type.ϖ_Cabannes[iB], 
                            Rayl𝐙⁺⁺, Rayl𝐙⁻⁺),
                linCoreScatteringOpticalProperties( 
                            arr_type(dtmp_τ[i]), 
                            dtmp_ϖ, 
                            dtmp_Z⁺⁺, dtmp_Z⁻⁺))
                for i=1:Nz]
        @show size(rayl)
        
        #@show size(dtmp_τ[end]), size(lin_rayl[end].lin_τ) 
        #@show dtmp_ϖ, lin_rayl[end].lin_ϖ
        #@show size(dtmp_Z⁺⁺), size(lin_rayl[end].lin_Z⁺⁺)           
        # Initiate combined properties with rayleigh
        combo = rayl #scatt_paired_xdx(rayl, lin_rayl)
        @show size(combo[end].x.τ),size(combo[end].dx.lin_τ)    
        #combo, lin_combo = rayl, lin_rayl

        #@show size(combo.x[1].τ)
        #@show size(combo.dx[1].lin_τ)
        #lin_combo = []

        # Loop over all aerosol types:
        for i=1:nAero
        # Precomute Z matrices per type (constant per layer)
        #@show typeof(aerosol_optics[iB][i].greek_coefs), typeof(pol_type), typeof(μ)
            AerZ⁺⁺, AerZ⁻⁺ = Scattering.compute_Z_moments(
                pol_type, μ, 
                aerosol_optics[iB][i].greek_coefs, 
                m, arr_type=arr_type)
            dAerZ⁺⁺ = zeros(4, size(AerZ⁺⁺,1), size(AerZ⁺⁺,2))                    
            dAerZ⁻⁺ = zeros(4, size(AerZ⁺⁺,1), size(AerZ⁺⁺,2))                    

            for ctr=1:4 
                dAerZ⁺⁺[ctr,:,:], dAerZ⁻⁺[ctr,:,:] = 
                        Scattering.compute_Z_moments(
                        pol_type, μ, 
                        lin_aerosol_optics[iB][i].d_greek_coefs[ctr], 
                        m, arr_type=arr_type)
            end

            #@show typeof(AerZ⁺⁺), typeof(aerosol_optics[iB][i]), typeof(FT.(τ_aer[iB][i,:]))
            # Generate Core optical properties for Aerosols i
            #aer, lin_aer   = 
            aer = createAero(τ_aer[iB][i,:], 
                (aerosol_optics[iB][i]), 
                AerZ⁺⁺, AerZ⁻⁺, 
                nparams, Nz, i,
                lin_τ_aer_psurf[iB][i,:], 
                lin_τ_aer_z₀[iB][i,:], 
                lin_τ_aer_σz[iB][i,:], 
                lin_aerosol_optics[iB][i], 
                dAerZ⁺⁺, dAerZ⁻⁺)

            #combo_aer = scatt_paired_xdx(aer, lin_aer)
            #@show size(τ_aer[iB][i,end])
            #@show size(aer), size(aer[end].τ) 
            #@show aerosol_optics[iB][i].ω̃, aer[end].ϖ
            #@show size(AerZ⁺⁺), size(aer[end].Z⁺⁺) 
            
            #@show size(dtmp_τ[end]), size(lin_aer[end].lin_τ) 
            #@show lin_aerosol_optics[iB][i].dω̃, lin_aer[end].lin_ϖ
            #@show size(dAerZ⁺⁺), size(lin_aer[end].lin_Z⁺⁺)   
            # aer -> CoreScatteringOpticalProperties(τ_mod, ϖ_mod,AerZ⁺⁺, AerZ⁻⁺),
            # lin_aer -> linCoreScatteringOpticalProperties(lin_τ_mod, lin_ϖ_mod, dAerZ⁺⁺, dAerZ⁻⁺)
            
                #@show typeof(aer), typeof(combo)
            # Mix with previous Core Optical Properties
            
            #lin_combo = lin_combo .+ lin_aer
            combo = combo .+ aer
            #combo = combo + combo_aer
            #(combo, lin_combo) = (combo, lin_combo) .+ (aer, lin_aer)
        end
        #@show typeof(combo)
        ###

        # fScattRayleigh:
        # Assume ϖ of 1 for Rayleight here:
        #@show size(combo)
        #fScattRayleigh = [FT.(Array(rayl[i].τ  ./ combo[i].τ)) for i=1:nZ]
        fScattRayleigh = [Array(rayl[i].τ  ./ combo[i].τ) for i=1:nZ]
        #lin_fScattRayleigh = [Array(-rayl[i].τ .* lin_combo[i].lin_τ  ./ combo[i].τ^2) for i=1:nZ]

        #@show fScattRayleigh, rayl[1].τ, combo[1].τ
        # Create Core Optical Properties merged with trace gas absorptions:
        #@show typeof(combo.+ 
        #[CoreAbsorptionOpticalProperties(arr_type((τ_abs[iB][:,i]))) for i=1:nZ])
        
        # Gaseous Absorption
        gabs = [abs_paired_xdx(CoreAbsorptionOpticalProperties(
                        arr_type(τ_abs[iB][:,i])), 
                linCoreAbsorptionOpticalProperties(
                        arr_type(lin_τ_abs[iB][:,:,i])))        
        for i=1:nZ]
        #lin_gabs = [linCoreAbsorptionOpticalProperties(
        #    arr_type(lin_τ_abs[iB][:,:,i])) for i=1:nZ]
        #combo_gabs = abs_paired_xdx(gabs, lin_gabs)
        #combo = combo .+ combo_gabs
        combo = combo .+ gabs
        #(combo, lin_combo) = (combo, lin_combo) .+ (gabs, lin_gabs)
        
        # Initiate combined properties with rayleigh
        #(combo, lin_combo) = (combo, lin_combo) .+ (gabs, lin_gabs)
        
        push!(band_layer_props, combo.x)
        push!(lin_band_layer_props, combo.dx)
        #push!(lin_band_layer_props, lin_combo)
        push!(band_fScattRayleigh,fScattRayleigh)
        #push!(lin_band_fScattRayleigh,lin_fScattRayleigh)
    end

    layer_opt = []
    fscat_opt = []
    lin_layer_opt = []
    lin_fscat_opt = []
    for iz = 1:nZ
        push!(layer_opt, prod([band_layer_props[i][iz] for i=1:length(iBand)]));
        push!(fscat_opt, [band_fScattRayleigh[i][iz] for i=1:length(iBand)]);
        push!(lin_layer_opt, prod([lin_band_layer_props[i][iz] for i=1:length(iBand)]));
        #push!(lin_fscat_opt, [lin_band_fScattRayleigh[i][:,iz] for i=1:length(iBand)]);
    end
    # For now just one band_fScattRayleigh
    return layer_opt, fscat_opt, lin_layer_opt#, lin_fscat_opt
end

function createAero(τAer, 
    aerosol_optics, 
    AerZ⁺⁺, AerZ⁻⁺, 
    nparams, Nz, iAer,
    lin_τ_aer_psurf, 
    lin_τ_aer_z₀,
    lin_τ_aer_σz,
    lin_aerosol_optics, 
    dAerZ⁺⁺, dAerZ⁻⁺) 

    @unpack k_ref, k, fᵗ, ω̃ = aerosol_optics
    @unpack dk_ref, dk, dfᵗ, dω̃ = lin_aerosol_optics
    τ_mod = (1-fᵗ * ω̃ ) * τAer;
    ϖ_mod = (1-fᵗ) * ω̃/(1-fᵗ * ω̃)
    FT = eltype(ϖ_mod)
    lin_τ_mod = zeros(5,length(τ_mod))
    lin_ϖ_mod = zeros(5)
    lin_τ_mod[1,:] = (τ_mod./τAer)*(k/k_ref)
    lin_ϖ_mod[1] = 0
    lin_τ_aer_psurf .*= (τ_mod./τAer)
    lin_τ_aer_z₀ .*= (τ_mod./τAer)
    lin_τ_aer_σz .*= (τ_mod./τAer)
    for ctr=2:5 #ctr=1 corresponds to the derivative with respect to τ_ref, the rest for the microphysical parameters nᵣ, nᵢ, r₀. and σ₀
        mctr = ctr-1
        lin_τ_mod[ctr,:] = τ_mod * 
            (dk[mctr]/k - dk_ref[mctr]/k_ref)
                - τAer * (fᵗ*dω̃[mctr] + dfᵗ[mctr]*ω̃)
        lin_ϖ_mod[ctr] = (dω̃[mctr] - 
            (dfᵗ[mctr]*ω̃ + fᵗ*dω̃[mctr])*(1-ϖ_mod))/
            (1-fᵗ * ω̃)
    end
    dtmp_τ = [zeros(FT, nparams) for i=1:Nz]
    dtmp_ϖ = zeros(FT, nparams)
    dtmp_Z⁺⁺ = zeros(FT, nparams, size(AerZ⁺⁺,1), size(AerZ⁺⁺,1))
    dtmp_Z⁻⁺ = zeros(FT, nparams, size(AerZ⁺⁺,1), size(AerZ⁺⁺,1))
    @show size(τ_mod[end]), length(τ_mod), size(ϖ_mod), size(lin_τ_mod), size(lin_ϖ_mod)
    @show nparams, Nz
    for iz = 1:Nz            
        dtmp_τ[iz][1,:] .= lin_τ_aer_psurf[iz]
        ctr = 9 + (iAer-1)*7 
        dtmp_τ[iz][ctr+6,:] .= lin_τ_aer_z₀[iz]
        dtmp_τ[iz][ctr+7,:] .= lin_τ_aer_σz[iz]
        dtmp_τ[iz][ctr+1,:] .= lin_τ_mod[1,iz]
        dtmp_ϖ[ctr+1] = lin_ϖ_mod[1]
        dtmp_Z⁺⁺[ctr+1,:,:] .= 0
        dtmp_Z⁻⁺[ctr+1,:,:] .= 0
        for i=2:5
            mctr=i-1
            dtmp_τ[iz][ctr+i,:] .= lin_τ_mod[mctr,iz]
            dtmp_ϖ[ctr+i] = lin_ϖ_mod[mctr]
            dtmp_Z⁺⁺[ctr+i,:,:] .= dAerZ⁺⁺[mctr,:,:]
            dtmp_Z⁻⁺[ctr+i,:,:] .= dAerZ⁻⁺[mctr,:,:]
        end
    end

    aer =  [scatt_paired_xdx(CoreScatteringOpticalProperties(
                    (τ_mod[i]), ϖ_mod, 
                    AerZ⁺⁺, AerZ⁻⁺), 
            linCoreScatteringOpticalProperties( 
                    (dtmp_τ[i]), 
                    dtmp_ϖ, 
                    dtmp_Z⁺⁺, dtmp_Z⁻⁺))        
            for i=1:Nz]

    #lin_aer =  [linCoreScatteringOpticalProperties( 
    #    (dtmp_τ[iz]), 
    #    dtmp_ϖ, 
    #    dtmp_Z⁺⁺, dtmp_Z⁻⁺) for iz=1:Nz]

    return aer #, lin_aer
    #CoreScatteringOpticalProperties(τ_mod, ϖ_mod, AerZ⁺⁺, AerZ⁻⁺),
    #linCoreScatteringOpticalProperties(dtmp_τ, dtmp_ϖ, dtmp_Z⁺⁺, dtmp_Z⁻⁺)

    #CoreScatteringOpticalProperties(τ_mod, ϖ_mod, AerZ⁺⁺, AerZ⁻⁺),
    #linCoreScatteringOpticalProperties(lin_τ_mod, lin_ϖ_mod, dAerZ⁺⁺, dAerZ⁻⁺)
end

# Extract scattering definitions and integrated absorptions for the source function!
function extractEffectiveProps(
    lods::Array,#{CoreScatteringOpticalProperties{FT},1}
    lin_lods::Array,
    quad_points::QuadPoints{FT}
    ) where FT

    #FT    = eltype(lods[1].τ)
    nSpec = length(lods[1].τ)
    Nz    = length(lods)
    nparams = size(lin_lods[1].lin_τ, 1)
    # First the Scattering Interfaces:
    scattering_interface = ScatteringInterface_00()
    scattering_interfaces_all = []
    τ_sum_all = similar(lods[1].τ,(nSpec,Nz+1)) #??
    lin_τ_sum_all = similar(lin_lods[1].lin_τ,(nparams,nSpec,Nz+1)) #??
    τ_sum_all[:,1] .= 0
    lin_τ_sum_all[:,:,1] .= 0
    #@show FT
    for iz =1:Nz
        # Need to check max entries in Z matrices here as well later!
        scatter = maximum(lods[iz].τ .* lods[iz].ϖ) > 2eps(FT)
        scattering_interface = get_scattering_interface(scattering_interface, 
                                scatter, iz)
        push!(scattering_interfaces_all, scattering_interface)
        #@show typeof(τ_sum_all[:,iz]), typeof(lods[iz].τ)
        @views τ_sum_all[:,iz+1] = τ_sum_all[:,iz] + 
            getG_atSun(lods[iz], quad_points) * lods[iz].τ 
        for ctr = 1:nparams
            @views lin_τ_sum_all[ctr,:,iz+1] = lin_τ_sum_all[ctr,:,iz] #+ getG_atSun(lods[iz], quad_points) * lods[iz].lin_τ[ctr,:] 
        end
    end
    return scattering_interfaces_all, τ_sum_all, lin_τ_sum_all
end
#=
function getG_atSun(lod::CoreScatteringOpticalProperties,quad_points::QuadPoints{FT}) where FT
    FT(1)
end

function getG_atSun(lod::CoreDirectionalScatteringOpticalProperties,quad_points::QuadPoints{FT}) where FT
    @unpack iμ₀ = quad_points
    gfct = Array(lod.G)[iμ₀]
    return gfct
end
=#

function expandOpticalProperties(in::linCoreScatteringOpticalProperties, arr_type)
    @unpack lin_τ, lin_ϖ, lin_Z⁺⁺, lin_Z⁻⁺ = in 
    @assert size(lin_τ) == size(lin_ϖ) "τ and ϖ sizes need to match"
    if size(lin_Z⁺⁺,4) == 1
        lin_Z⁺⁺ = _repeat(lin_Z⁺⁺,1,1,1,length(lin_τ[1,:]))
        lin_Z⁻⁺ = _repeat(lin_Z⁻⁺,1,1,1,length(lin_τ[1,:]))
        return linCoreScatteringOpticalProperties(arr_type(τ), arr_type(ϖ), arr_type(Z⁺⁺), arr_type(Z⁻⁺)) 
    else
        @assert size(lin_Z⁺⁺,4) ==  length(lin_τ[1,:]) "Z and τ dimensions need to match "
        linCoreScatteringOpticalProperties(arr_type(lin_τ), arr_type(lin_ϖ), arr_type(lin_Z⁺⁺), arr_type(lin_Z⁻⁺)) 
    end
end

function expandOpticalProperties_lin(in::linCoreDirectionalScatteringOpticalProperties, arr_type)
    @unpack lin_τ, lin_ϖ, lin_Z⁺⁺, lin_Z⁻⁺, lin_G = in 
    @assert size(lin_τ) == size(lin_ϖ) "τ and ϖ sizes need to match"
    if size(lin_Z⁺⁺,5) == 1
        lin_Z⁺⁺ = _repeat(lin_Z⁺⁺,1,1,1,1,length(lin_τ[1,:]))
        lin_Z⁻⁺ = _repeat(lin_Z⁻⁺,1,1,1,1,length(lin_τ[1,:]))
        return linCoreDirectionalScatteringOpticalProperties(
            arr_type(lin_τ), arr_type(lin_ϖ), 
            arr_type(lin_Z⁺⁺), arr_type(lin_Z⁻⁺), 
            arr_type(lin_G)) 
    else
        @assert size(lin_Z⁺⁺,5) ==  length(lin_τ[1,:]) "Z and τ dimensions need to match "
        linCoreDirectionalScatteringOpticalProperties(
            arr_type(lin_τ), arr_type(lin_ϖ), 
            arr_type(lin_Z⁺⁺), arr_type(lin_Z⁻⁺), 
            arr_type(lin_G)) 
    end
end

#=
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
=#
#expandScalar(x,n) = x.*ones(n);