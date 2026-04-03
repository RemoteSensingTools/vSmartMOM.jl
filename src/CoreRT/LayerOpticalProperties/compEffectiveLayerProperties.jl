function constructCoreOpticalProperties(RS_type, iBand, m, model)
    @unpack τ_rayl, τ_aer, τ_abs, aerosol_optics, 
            greek_rayleigh, greek_cabannes, ϖ_Cabannes = model
    
    @assert all(iBand .≤ length(τ_rayl)) "iBand exceeded number of bands"
    FT = eltype(τ_rayl[1])
    arr_type = array_type(model.params.architecture)

    pol_type = model.params.polarization_type
    # Convert pol_type to match FT if needed
    if FT != eltype(pol_type.D)
        PT = Base.typename(typeof(pol_type)).wrapper
        pol_type = PT{FT}(n=pol_type.n, D=FT.(pol_type.D), I₀=FT.(pol_type.I₀))
    end
    # Do this in CPU space only first:
    
    # Quadrature points:
    μ = Array(model.quad_points.qp_μ )
    # Number of Aerosols:
    nAero = size(τ_aer[iBand[1]],1)
    nZ    = size(τ_rayl[1],2)
    #@show greek_rayleigh
    # Rayleigh Z matrix:
    
                                                        #@show Rayl𝐙⁺⁺

    band_layer_props    = [];
    band_fScattRayleigh = [];
    # @show arr_type
    for iB in iBand
        if (typeof(RS_type)<:noRS) #!(typeof(RS_type)<:Union{RRS,RRS_plus})
            Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, μ,
                                                            greek_rayleigh[iB], m,
                                                            arr_type = arr_type);
        else
            Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, μ,
                                                            greek_cabannes[iB], m,
                                                            arr_type = arr_type);
            #Rayl2𝐙⁺⁺, Rayl2𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, μ, 
            #                                                greek_rayleigh[iB], m, 
            #                                                arr_type = arr_type);
        end
        # Convert Z-matrices to FT if needed (compute_Z_moments may return Float64)
        if eltype(Rayl𝐙⁺⁺) != FT
            Rayl𝐙⁺⁺ = FT.(Rayl𝐙⁺⁺)
            Rayl𝐙⁻⁺ = FT.(Rayl𝐙⁻⁺)
        end

        if (typeof(RS_type)<:noRS) #if !(typeof(RS_type)<:Union{RRS,RRS_plus})
            rayl =  [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]), FT(1),
                (Rayl𝐙⁺⁺), (Rayl𝐙⁻⁺)) for i=1:nZ]
        else
            rayl =  [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]), FT(ϖ_Cabannes[iB]),
                (Rayl𝐙⁺⁺), (Rayl𝐙⁻⁺)) for i=1:nZ]
            #@show τ_rayl[iB][1,i]
            #rayl2 =  [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]), 1.0, 
            #    (Rayl2𝐙⁺⁺), (Rayl2𝐙⁻⁺)) for i=1:nZ]
        end
        #@show τ_rayl[iB][1,1], τ_rayl[iB][1,end]
        #@show τ_aer[iB][1,1,1], τ_aer[iB][1,1,end]
        #CoreScatteringOpticalProperties.(
        #        τ_rayl[iB], 
        #        [RS_type.ϖ_Cabannes[iB]], 
        #        [Rayl𝐙⁺⁺], [Rayl𝐙⁻⁺])
        
        #@show size(rayl)
        # Initiate combined properties with rayleigh
        combo = rayl
        # test:
        # combo = combo .+ rayl
        # this throws the following error:
        # ERROR: MethodError: Cannot `convert` an object of type 
        #  vSmartMOM.CoreRT.CoreScatteringOpticalProperties{CuArray{Float64{},1,CUDA.Mem.DeviceBuffer{}},CuArray{Float64, 1, CUDA.Mem.DeviceBuffer},CuArray{Float64{},3,CUDA.Mem.DeviceBuffer{}}} to an object of type 
        #  vSmartMOM.CoreRT.CoreScatteringOpticalProperties{CuArray{Float64{},1,CUDA.Mem.DeviceBuffer{}},Float64,CuArray{Float64{},2,CUDA.Mem.DeviceBuffer{}}}
        # Closest candidates are:
        #  convert(::Type{T}, ::T) where T
        #   @ Base Base.jl:64
        #  (::Type{vSmartMOM.CoreRT.CoreScatteringOpticalProperties{FT, FT2, FT3}} where {FT, FT2, FT3})(::Any, ::Any, ::Any, ::Any)
        #   @ vSmartMOM ~/code/github/vSmartMOM.jl/src/CoreRT/types.jl:605

        #@show combo[1].τ[1], combo[1].τ[end]
        #@show combo[1].ϖ
        #@show RS_type.ϖ_Cabannes
        # Loop over all aerosol types:
        for iaer=1:nAero
            # Precomute Z matrices per type (constant per layer)
            #@show iB,i
            AerZ⁺⁺, AerZ⁻⁺ = Scattering.compute_Z_moments(
                                pol_type, μ,
                                aerosol_optics[iB][iaer].greek_coefs,
                                m, arr_type=arr_type)
            # Convert to FT if needed
            if eltype(AerZ⁺⁺) != FT
                AerZ⁺⁺ = FT.(AerZ⁺⁺)
                AerZ⁻⁺ = FT.(AerZ⁻⁺)
            end
            # Generate Core optical properties for Aerosols iaer
            #@show size(τ_aer[iB][iaer,:,:])
            #aer = Vector{CoreScatteringOpticalProperties}
            #aer =  [CoreScatteringOpticalProperties(zeros(length(τ_rayl[iB][:,1])), zeros(length(τ_rayl[iB][:,1])), 
            #    zeros(size(Rayl𝐙⁺⁺)), zeros(size(Rayl𝐙⁻⁺))) for i=1:nZ]
            #for i=1:nZ   
                #aer[i]   = createAero(τ_aer[iB][iaer,:,i], 
                #                aerosol_optics[iB][iaer], 
                #                AerZ⁺⁺, AerZ⁻⁺)
            #    push!(aer, createAero(τ_aer[iB][iaer,:,i], 
            #                    aerosol_optics[iB][iaer], 
            #                    AerZ⁺⁺, AerZ⁻⁺))                
            #end
            #aer =  [createAero(arr_type(τ_aer[iB][iaer,:,i]), 
            #            aerosol_optics[iB][iaer], 
            #            AerZ⁺⁺, AerZ⁻⁺ ) for i=1:nZ]
            aer =  [createAero(arr_type(τ_aer[iB][iaer,:,i]), 
                        aerosol_optics[iB][iaer], 
                        AerZ⁺⁺, AerZ⁻⁺, arr_type) for i=1:nZ]
                
            #@show aer[1].τ[1], aer[1].τ[end]
            #@show size(aer[end].τ), aer[end].τ[1], aer[end].τ[end]
            #@show size(aer[end].ϖ), aer[end].ϖ[1], aer[end].ϖ[end]
            #@show τ_aer[iB][iaer,:,:]
            # Mix with previous Core Optical Properties
            #@show combo[1].ϖ   , aer[1].ϖ
            #@show typeof(combo)
            #@show typeof(aer)
            combo = combo .+ aer

            #@show combo[1].ϖ   , aer[1].ϖ
        end

        # Somewhere here we can add canopy later as well!
        ###

        # fScattRayleigh:
        #@show rayl[1].τ * rayl[1].ϖ, combo[1].τ
        # Assume ϖ of 1 for Rayleight here:
        #@show size(combo)
        #fScattRayleigh = [Array(rayl[i].τ  ./ combo[i].τ) for i=1:nZ]
        #@show fScattRayleigh, rayl[1].τ, combo[1].τ
        # Create Core Optical Properties merged with trace gas absorptions:
        #@show size(combo)
        
        #@show size(fScattRayleigh)
        #@show size(combo[1].τ), size(τ_abs[iB][:,1])

        combo2 = combo .+ [CoreAbsorptionOpticalProperties(arr_type(FT.(τ_abs[iB][:,i]))) for i=1:nZ]
        # Ensure all layer properties are in FT precision
        if eltype(combo2[1].τ) != FT
            combo2 = [CoreScatteringOpticalProperties{FT}(FT.(c.τ), FT.(c.ϖ), FT.(c.Z⁺⁺), FT.(c.Z⁻⁺)) for c in combo2]
        end
        fScattRayleigh = [Array(rayl[i].τ  ./ combo2[i].τ) for i=1:nZ]
        #@show fScattRayleigh[1]
        #for i=1:nZ
        #    @show i, rayl[i].τ, combo[1].τ#,combo2[1].τ
        #end
        push!(band_layer_props,combo2 )
        push!(band_fScattRayleigh,fScattRayleigh)
        #aType = array_type(model.params.architecture)
        #combo2 = [CoreScatteringOpticalProperties(aType(combo[i].τ),aType(combo[i].ϖ), aType(combo[i].Z⁺⁺), aType(combo[i].Z⁻⁺)) for i in eachindex(combo)]
        # Need to check how to convert to GPU later as well!
        #return combo,fScattRayleigh
        #@show rayl[1].τ 
        #@show rayl[1].ϖ
        #@show rayl[1].Z⁺⁺
        #@show typeof(rayl[1].τ)
        #@show Array(rayl[1].τ)[1] * rayl[1].ϖ * Array(rayl[1].Z⁺⁺)
        #@show Array(rayl[1].τ)[1] * sum(RS_type.ϖ_λ₁λ₀) * Array(RS_type.Z⁺⁺_λ₁λ₀) 
        #@show Array(rayl2[1].τ)[1] * rayl2[1].ϖ * Array(rayl2[1].Z⁺⁺)
    
        #=@show sum(Array(rayl[1].τ)[1] * rayl[1].ϖ * Array(rayl[1].Z⁺⁺) + 
        Array(rayl[1].τ)[1] * sum(RS_type.ϖ_λ₁λ₀) * Array(RS_type.Z⁺⁺_λ₁λ₀) - 
        Array(rayl2[1].τ)[1] * rayl2[1].ϖ * Array(rayl2[1].Z⁺⁺), dims=1)
        @show sum(Array(rayl[1].τ)[1] * rayl[1].ϖ * Array(rayl[1].Z⁻⁺) + 
        Array(rayl[1].τ)[1] * sum(RS_type.ϖ_λ₁λ₀) * Array(RS_type.Z⁻⁺_λ₁λ₀) - 
        Array(rayl2[1].τ)[1] * rayl2[1].ϖ * Array(rayl2[1].Z⁻⁺), dims=1)
        =#
        #@show rayl2[1].Z⁺⁺[:,:,1] #.==0
        #@show rayl[1].Z⁺⁺[:,:,1]
        #@show RS_type.Z⁺⁺_λ₁λ₀[:,:,1] #.==0

    end
    #bla
    #@show RS_type.bandSpecLim[1]
    #@show RS_type.iBand
    layer_opt = []
    fscat_opt = []
    for iz = 1:nZ
        push!(layer_opt, prod([band_layer_props[i][iz] for i=1:length(iBand)]));
        #push!(fscat_opt, expandBandScalars(RS_type,[band_fScattRayleigh[i][iz] for i=1:length(iBand)]));
        push!(fscat_opt, [band_fScattRayleigh[i][iz] for i=1:length(iBand)]);
    end
    # For now just one band_fScattRayleigh
    #@show typeof(layer_opt[1].τ)
    return layer_opt, fscat_opt # Suniti: this needs to be modified because Rayleigh scattering fraction varies dramatically with wavelength
end

function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺)
    @unpack fᵗ, ω̃ = aerosol_optics
    #τ_mod = (1-fᵗ * ω̃ ) * τAer;
    #ϖ_mod = (1-fᵗ) * ω̃/(1-fᵗ * ω̃)
    @show typeof(fᵗ), typeof(ω̃)
    τ_mod = (1 .- fᵗ * ω̃ ) .* τAer;
    ϖ_mod = (1 .- fᵗ) .* ω̃ ./ (1 .- fᵗ * ω̃)
    CoreScatteringOpticalProperties(τ_mod, ϖ_mod,AerZ⁺⁺, AerZ⁻⁺)
end

function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺, arr_type)
    @unpack fᵗ, ω̃ = aerosol_optics
    if ω̃ isa Number
        nothing
    else
        ω̃ = arr_type(aerosol_optics.ω̃) 
    end
    #@show typeof(ω̃), typeof(fᵗ)
    #@show size(fᵗ)
    #@show size(ω̃)
    #@show size(τAer), τAer[1], τAer[end]
    #τ_mod = zeros(size(τAer,1), size(τAer,2))
    #for iz = 1:size(τAer,1)
    τ_mod = (1 .- fᵗ * ω̃ ) .* τAer;
    #@show τ_mod[1], τ_mod[end]  
    #end
    ϖ_mod = (1 .- fᵗ) .* ω̃ ./ (1 .- fᵗ * ω̃)
    CoreScatteringOpticalProperties(τ_mod, ϖ_mod, AerZ⁺⁺, AerZ⁻⁺)
end

# Extract scattering definitions and integrated absorptions for the source function!
function extractEffectiveProps(
                lods::Array#{CoreScatteringOpticalProperties{FT},1}
                ) #where FT

    FT    = eltype(lods[1].τ)
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
        @views τ_sum_all[:,iz+1] = τ_sum_all[:,iz] + lods[iz].τ 
    end
    return scattering_interfaces_all, τ_sum_all
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