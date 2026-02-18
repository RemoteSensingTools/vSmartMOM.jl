"""
    constructCoreOpticalProperties(RS_type, iBand, m, model, lin_model)

Construct the combined effective layer optical properties and their linearized derivatives
for all atmospheric layers, merging Rayleigh scattering, aerosol scattering (with δ-M
truncation), and trace gas absorption.

For each atmospheric layer, this function:
1. Creates Rayleigh scattering properties (τ_rayl, ϖ=1, Z_rayl).
2. For each aerosol type, computes δ-M scaled optical properties via [`createAero`](@ref),
   combining Mie scattering with the aerosol vertical profile.
3. Adds Rayleigh and aerosol properties using the `+` operator on
   `UmbrellaCoreScatteringOpticalProperties`, which correctly propagates derivatives through
   the mixing formulas for combined ϖ and Z.
4. Adds gas absorption via `UmbrellaCoreAbsorptionOpticalProperties`, updating ϖ and τ.

The output derivative dimension has `Nparams = 7×NAer + NGas` entries per layer
(surface parameters are handled separately in the interaction step).

# Returns
- `layer_opt`: Array of `CoreScatteringOpticalProperties` (one per layer).
- `layer_opt_lin`: Array of `CoreScatteringOpticalPropertiesLin` (derivatives, one per layer).
- `fscat_opt`: Rayleigh scattering fraction per layer (for inelastic scattering weight).
"""
function constructCoreOpticalProperties(RS_type, iBand, m, model, lin_model) #where {FT<:Real}
    @unpack τ_rayl, τ_aer, τ_abs, aerosol_optics, 
            greek_rayleigh, greek_cabannes, ϖ_Cabannes = model
    @unpack τ̇_aer, τ̇_abs, lin_aerosol_optics = lin_model
    @assert all(iBand .≤ length(τ_rayl)) "iBand exceeded number of bands"
    FT = eltype(τ_rayl[1])
    
    # Debug: Check what architecture and array_type we're getting
    arr_type = array_type(model.params.architecture)
    
    pol_type = model.params.polarization_type
    # Do this in CPU space only first:
    
    # Quadrature points:
    μ = collect(model.quad_points.qp_μ )
    # Number of Aerosols:
    nAero = size(τ_aer[iBand[1]],1)
    nZ    = size(τ_rayl[1],2)
    #@show greek_rayleigh
    # Rayleigh Z matrix:
    
                                                        #@show Rayl𝐙⁺⁺

    band_layer_props     = [];
    band_layer_props_lin = [];
    band_fScattRayleigh  = [];
    # @show arr_type
    for iB in iBand
        #@show iB
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

        #if (typeof(RS_type)<:noRS) #if !(typeof(RS_type)<:Union{RRS,RRS_plus})
        
        # Debug the exact line that's causing the error
        #@show arr_type
        #@show typeof(arr_type)
        #@show τ_rayl[iB][:,1]  # Check one element to see what we're passing
        #@show typeof(τ_rayl[iB][:,1]), size(τ_rayl)
        
        #@show Cucollect(τ_rayl[iB][:,1])  # Convert to CuArray if needed
        #@show typeof(Rayl𝐙⁺⁺), size(Rayl𝐙⁺⁺)
        CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,1]), FT(1.0), 
                (Rayl𝐙⁺⁺), (Rayl𝐙⁻⁺))

              
        rayl =  [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]), 1.0, 
                (Rayl𝐙⁺⁺), (Rayl𝐙⁻⁺)) for i=1:nZ]
            #rayl_lin = [CoreScatteringOpticalPropertiesLin(arr_type(τ̇_rayl[iB][:,iz]), 0.0, 
            #(0.0.*Rayl𝐙⁺⁺), (0.0.*Rayl𝐙⁻⁺)) for iz=1:nZ]    
        #else
        #    @error("Cannot linearize Raman computations")
            #rayl =  [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]), ϖ_Cabannes[iB], 
            #    (Rayl𝐙⁺⁺), (Rayl𝐙⁻⁺)) for i=1:nZ]
            #@show τ_rayl[iB][1,i]
            #rayl2 =  [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]), 1.0, 
            #    (Rayl2𝐙⁺⁺), (Rayl2𝐙⁻⁺)) for i=1:nZ]
        #end
        #@show τ_rayl[iB][1,1], τ_rayl[iB][1,end]
        #@show τ_aer[iB][1,1,1], τ_aer[iB][1,1,end]
        #CoreScatteringOpticalProperties.(
        #        τ_rayl[iB], 
        #        [RS_type.ϖ_Cabannes[iB]], 
        #        [Rayl𝐙⁺⁺], [Rayl𝐙⁻⁺])
        
        #@show size(rayl)
        # Initiate combined properties with rayleigh
        #combo = rayl
        combrella = [UmbrellaCoreScatteringOpticalProperties(rayl[i],nothing) for i=1:nZ]
        #@show 1
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
            AerZ⁺⁺, AerZ⁻⁺, AerŻ⁺⁺, AerŻ⁻⁺ = Scattering.compute_Z_moments(
                                pol_type, μ, 
                                aerosol_optics[iB][iaer].greek_coefs, 
                                lin_aerosol_optics[iB][iaer].lin_greek_coefs, 
                                m, arr_type=arr_type)
            #@show AerŻ⁺⁺, size(AerŻ⁺⁺)
            #@show AerŻ⁻⁺, size(AerŻ⁻⁺)
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
            aer = []
            lin_aer = []
            for iz=1:nZ
                #@show size(τ_aer[iB][iaer,:,iz])
                #@show size(τ̇_aer[iB][iaer,:,:,iz])
                #@show size(aerosol_optics[iB][iaer].fᵗ), size(aerosol_optics[iB][iaer].ω̃)
                #@show size(lin_aerosol_optics[iB][iaer].ḟᵗ), size(lin_aerosol_optics[iB][iaer].ω̃̇)   
                t_aer, t_lin_aer =  createAero(
                            arr_type(τ_aer[iB][iaer,:,iz]), 
                            aerosol_optics[iB][iaer], 
                            arr_type(AerZ⁺⁺), arr_type(AerZ⁻⁺), 
                            arr_type(τ̇_aer[iB][iaer,:,:,iz]), 
                            lin_aerosol_optics[iB][iaer], 
                            arr_type(AerŻ⁺⁺), arr_type(AerŻ⁻⁺), 
                            arr_type) 
                #@show iz, size(t_aer.τ)
                #@show iz, size(t_aer.ϖ)

                #@show iz, size(t_lin_aer.ϖ̇)
                #@show iz, size(t_lin_aer.τ̇)
                push!(aer, t_aer)
                push!(lin_aer, t_lin_aer)
            end 
            aer_combrella = [UmbrellaCoreScatteringOpticalProperties(aer[i],lin_aer[i]) for i=1:nZ]   
            #@show 2
            #@show aer[1].τ[1], aer[1].τ[end]
            #@show size(aer[end].τ), aer[end].τ[1], aer[end].τ[end]
            #@show size(aer[end].ϖ), aer[end].ϖ[1], aer[end].ϖ[end]
            #@show τ_aer[iB][iaer,:,:]
            # Mix with previous Core Optical Properties
            #@show combo[1].ϖ   , aer[1].ϖ
            #@show typeof(combo)
            #@show typeof(aer)
            #combo = combo .+ aer
            #combrella = combrella .+ aer_combrella
            tmp = combrella[1]+aer_combrella[1]
            #@show 2.1
            tmp = [combrella[i]+aer_combrella[i] for i=1:nZ]
            #@show 3
            combrella = tmp
            #@show 4
            #@show combo[1].ϖ   , aer[1].ϖ
        end

        # Somewhere here we can add canopy later as well!
        ###

        # fScattRayleigh:
        #@show rayl[1].τ * rayl[1].ϖ, combo[1].τ
        # Assume ϖ of 1 for Rayleight here:
        #@show size(combo)
        #fScattRayleigh = [collect(rayl[i].τ  ./ combo[i].τ) for i=1:nZ]
        #@show fScattRayleigh, rayl[1].τ, combo[1].τ
        # Create Core Optical Properties merged with trace gas absorptions:
        #@show size(combo)
        
        #@show size(fScattRayleigh)
        #@show size(combo[1].τ), size(τ_abs[iB][:,1])
        gas = [CoreAbsorptionOpticalProperties(arr_type(τ_abs[iB][:,iz])) for iz=1:nZ]
        lin_gas = [CoreAbsorptionOpticalPropertiesLin(arr_type(τ̇_abs[iB][:,:,iz])) for iz=1:nZ] 
        #combo2 = combo .+ gas
        gas_combrella = [UmbrellaCoreAbsorptionOpticalProperties(gas[iz],lin_gas[iz]) for iz=1:nZ]   
        tmp = combrella .+ gas_combrella            
        combrella = tmp
        #@show size(combo2[1].τ)
        #fScattRayleigh = [collect(rayl[iz].τ  ./ combo2[iz].τ) for iz=1:nZ] 
        fScattRayleigh = [collect(rayl[iz].τ  ./ combrella[iz].fwd.τ) for iz=1:nZ] 
        
        #@show fScattRayleigh[1]
        #for i=1:nZ
        #    @show i, rayl[i].τ, combo[1].τ#,combo2[1].τ
        #end
        #combo_lin = [include_rayl!(combo[iz], combo_lin[iz], rayl[iz], rayl_lin[iz]) for iz=1:nZ]
        #for iaer=1:nAero
        #    aer_lin =  [createAeroLin(arr_type(τ̇_aer[iB][iaer,1:7,:,iz]), 
        #                aerosol_optics_lin[iB][iaer], 
        #                AerŻ⁺⁺, AerŻ⁻⁺, arr_type) for iz=1:nZ]
        #    combo_lin = [include_aer!(iaer, combo[iz], combo_lin[iz], aer[iz], aer_lin[iz]) for i=1:nZ]
        #end
        # Use the following two lines for every gas to be included in the state vector
        #igas=1
        #combo2_lin = [include_gas(nAero, igas[iz], combo[iz], combo_lin[iz], gas_lin[iz]) for i=1:nZ]
        combo =     [combrella[iz].fwd for iz=1:nZ]
        combo_lin = [combrella[iz].lin for iz=1:nZ]

        push!(band_layer_props, combo)
        push!(band_layer_props_lin, combo_lin)
        push!(band_fScattRayleigh,fScattRayleigh)
        #aType = array_type(model.params.architecture)
        #combo2 = [CoreScatteringOpticalProperties(aType(combo[i].τ),aType(combo[i].ϖ), aType(combo[i].Z⁺⁺), aType(combo[i].Z⁻⁺)) for i in eachindex(combo)]
        # Need to check how to convert to GPU later as well!
        #return combo,fScattRayleigh
        #@show rayl[1].τ 
        #@show rayl[1].ϖ
        #@show rayl[1].Z⁺⁺
        #@show typeof(rayl[1].τ)
        #@show collect(rayl[1].τ)[1] * rayl[1].ϖ * collect(rayl[1].Z⁺⁺)
        #@show collect(rayl[1].τ)[1] * sum(RS_type.ϖ_λ₁λ₀) * collect(RS_type.Z⁺⁺_λ₁λ₀) 
        #@show collect(rayl2[1].τ)[1] * rayl2[1].ϖ * collect(rayl2[1].Z⁺⁺)
    
        #@show rayl2[1].Z⁺⁺[:,:,1] #.==0
        #@show rayl[1].Z⁺⁺[:,:,1]
        #@show RS_type.Z⁺⁺_λ₁λ₀[:,:,1] #.==0

    end
    #bla
    #@show RS_type.bandSpecLim[1]
    #@show RS_type.iBand
    layer_opt = []
    layer_opt_lin = []
    fscat_opt = []
    for iz = 1:nZ
        push!(layer_opt, prod([band_layer_props[i][iz] for i=1:length(iBand)])); #Is this intentional?? Why prod?
        push!(layer_opt_lin, prod([band_layer_props_lin[i][iz] for i=1:length(iBand)])); #Is this intentional?? Why prod?
        #push!(fscat_opt, expandBandScalars(RS_type,[band_fScattRayleigh[i][iz] for i=1:length(iBand)]));
        push!(fscat_opt, [band_fScattRayleigh[i][iz] for i=1:length(iBand)]);
    end
    # For now just one band_fScattRayleigh
    #@show typeof(layer_opt[1].τ)
    return layer_opt, layer_opt_lin, fscat_opt # Suniti: this needs to be modified because Rayleigh scattering fraction varies dramatically with wavelength
end
 
function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺,
                    τ̇Aer, lin_aerosol_optics, AerŻ⁺⁺, AerŻ⁻⁺)
    @unpack fᵗ, ω̃ = aerosol_optics
    @unpack ḟᵗ, ω̃̇ = lin_aerosol_optics
    #τ_mod = (1-fᵗ * ω̃ ) * τAer;
    #ϖ_mod = (1-fᵗ) * ω̃/(1-fᵗω̃)
    #τ̇_mod = (1-fᵗ * ω̃ ) * τ̇Aer - (ḟᵗϖ+fᵗϖ̇) * τAer;
    #ϖ̇_mod = [ϖ̇(1-fᵗ) - ḟᵗϖ(1-ϖ)]/(1-fᵗω̃)²
    #@show typeof(fᵗ), typeof(ω̃)
    
    τ_mod = (1 .- fᵗ * ω̃ ) .* τAer;
    ϖ_mod = (1 .- fᵗ) .* ω̃ ./ (1 .- fᵗ * ω̃)
    τ̇_mod = similar(7, length(τAer))
    ϖ̇_mod = similar(7, length(ω̃))
    Ż⁺⁺_mod = similar(7, size(AerZ⁺⁺,1), size(AerZ⁺⁺,2))
    Ż⁻⁺_mod = similar(7, size(AerZ⁻⁺,1), size(AerZ⁻⁺,2))
    #Derivatives with respect to τAer
    τ̇_mod[1,:] = (1 .- fᵗ * ω̃ )
    for iparam=1:4
        tmp = fᵗ.*ω̃̇[iparam,:] .+ ḟᵗ[iparam,:].*ω̃

        τ̇_mod[1+iparam,:] = (1 .- fᵗ .* ω̃ ) .* τ̇Aer[iparam,:];
        τ̇_mod[1+iparam,:] .-= tmp .* τAer
        ϖ̇_mod[1+iparam,:] = (ω̃̇[iparam,:].*(1 .- fᵗ) .- ḟᵗ[iparam,:].*ω̃.*(1 .- ω̃))
        ϖ̇_mod[1+iparam,:] ./= (1 .- fᵗ * ω̃).^2
        Ż⁺⁺_mod[1+iparam,:,:] .= AerŻ⁺⁺[iparam,:,:]
        Ż⁻⁺_mod[1+iparam,:,:] .= AerŻ⁻⁺[iparam,:,:]
    end
    for iparam=5:6
        #tmp = 0 #fᵗ.*ω̃̇[iparam,:] .+ ḟᵗ[iparam,:].*ω̃
        τ̇_mod[1+iparam,:] = (1 .- fᵗ .* ω̃ ) .* τ̇Aer[iparam,:];
        ϖ̇_mod[1+iparam,:] .= 0
        Ż⁺⁺_mod[1+iparam,:,:] .= 0
        Ż⁻⁺_mod[1+iparam,:,:] .= 0
        #τ̇_mod[1+iparam,:] .-= tmp .* τAer
        #ϖ̇_mod[1+iparam,:] = (ω̃̇[iparam,:].*(1 .- fᵗ) .- ḟᵗ[iparam,:].*ω̃.*(1 .- ω̃))
        #ϖ̇_mod[1+iparam,:] ./= (1 .- fᵗ * ω̃).^2
    end
    return CoreScatteringOpticalProperties(τ_mod, ϖ_mod, AerZ⁺⁺, AerZ⁻⁺),
        CoreScatteringOpticalPropertiesLin(τ̇_mod, ϖ̇_mod, Ż⁺⁺_mod, Ż⁻⁺_mod)
end

"""
    createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺, τ̇Aer, lin_aerosol_optics, AerŻ⁺⁺, AerŻ⁻⁺, arr_type)

Compute **δ-M scaled** aerosol optical properties and their derivatives for one aerosol
type in one atmospheric layer.

Applies the δ-M truncation correction (Nakajima & Tanaka, 1988; Sanghavi & Stephens, 2013)
to the aerosol optical depth and single-scattering albedo, and computes derivatives with
respect to 7 aerosol sub-parameters: `[τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]`.

# δ-M scaling formulas

```math
\\tau_\\text{mod} = (1 - f^t \\tilde{\\omega}) \\cdot \\tau_\\text{aer}
```
```math
\\varpi_\\text{mod} = \\frac{(1 - f^t) \\tilde{\\omega}}{1 - f^t \\tilde{\\omega}}
```

# Derivative chain rule

For each physical parameter ``p_j``:
```math
\\frac{\\partial \\tau_\\text{mod}}{\\partial p_j} = 
  (1 - f^t \\tilde{\\omega}) \\frac{\\partial \\tau_\\text{aer}}{\\partial p_j}
  - \\tau_\\text{aer} \\left(f^t \\frac{\\partial \\tilde{\\omega}}{\\partial p_j} + 
    \\tilde{\\omega} \\frac{\\partial f^t}{\\partial p_j}\\right)
```

For `τ_ref, p₀, σp`: only the `τ_aer` chain contributes (Mie properties are independent).

!!! note "Bug 18 fix"
    Corrected index mapping: `τ̇Aer[k,:]` now correctly maps to parameter `k`
    (was previously off-by-one, mixing τ_ref with Mie parameter derivatives).

# Returns
- `CoreScatteringOpticalProperties`: Forward δ-M scaled properties.
- `CoreScatteringOpticalPropertiesLin`: Linearized properties (7 sub-params).
"""
function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺, 
                    τ̇Aer, lin_aerosol_optics, AerŻ⁺⁺, AerŻ⁻⁺,
                    arr_type)

    @unpack fᵗ, ω̃ = aerosol_optics
    @unpack ḟᵗ, ω̃̇ = lin_aerosol_optics

    n  = size(τAer,1)
    #fᵗ = arr_type(fᵗ)
    #ω̃  = arr_type(ω̃)
    #ḟᵗ = arr_type(ḟᵗ)
    #ω̃̇  = arr_type(ω̃̇)

    sz = size(AerŻ⁺⁺)
    tmpŻ⁺⁺ = AerŻ⁺⁺
    tmpŻ⁻⁺ = AerŻ⁻⁺
    AerŻ⁺⁺ = arr_type(zeros(7, size(AerZ⁺⁺,1), size(AerZ⁺⁺,2)))#, n)
    AerŻ⁻⁺ = arr_type(zeros(7, size(AerZ⁻⁺,1), size(AerZ⁻⁺,2)))#, n)
    AerŻ⁺⁺[2:5,:,:] .= tmpŻ⁺⁺
    AerŻ⁻⁺[2:5,:,:] .= tmpŻ⁻⁺

    # Ensure arrays are in the right memory space (CPU or GPU)
    ω̃  = (ω̃ isa Number) ? arr_type(fill(ω̃,n)) : arr_type(ω̃)
    #fᵗ  = (fᵗ isa Number) ? arr_type(fill(fᵗ,n)) : arr_type(fᵗ)
    
    ω̃̇ = arr_type(ω̃̇)
    ḟᵗ = arr_type(ḟᵗ)
    
    sz = size(ω̃̇)
    if ndims(ω̃̇) == 1
        tmpω̃̇ = reshape(ω̃̇, sz..., 1) .* ones(eltype(ω̃̇), 1,n)        
    elseif ndims(ω̃̇) == 2 && sz[2] == 1
        tmpω̃̇ = ω̃̇ .* ones(eltype(ω̃̇), 1,n)
    else
        tmpω̃̇ = ω̃̇
    end
    ω̃̇ = arr_type(zeros(7,n))
    ω̃̇[2:5,:] .= tmpω̃̇
    
    tmpḟᵗ = ḟᵗ
    ḟᵗ = arr_type(zeros(7,n))
    ḟᵗ[2:5] = tmpḟᵗ

    # Forward modified properties
    τ_mod = (1 .- fᵗ * ω̃) .* τAer
    ϖ_mod = (1 .- fᵗ) .* ω̃ ./ (1 .- fᵗ * ω̃)

    # Allocate linearized outputs
    τ̇_mod = arr_type(zeros(7,n)) #similar(τ̇Aer)
    ϖ̇_mod = arr_type(zeros(7,n)) #similar(ω̃̇)

    # Bug 18 fix: derivatives w.r.t. 7 aerosol sub-params [τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]
    τ̇_mod[1,:] .= (1 .- fᵗ * ω̃) .* τ̇Aer[1,:]
    ϖ̇_mod[1,:] .= 0.0
    # Vectorized form over iparam dimension
    # Dimensions: iparam × spectral
    #@show size(fᵗ * ω̃̇), size(ḟᵗ .* ω̃')  
    tmp = fᵗ * ω̃̇[2:5,:] .+ ḟᵗ[2:5] * ω̃'  # (iparam, :)
    #@show size(tmp), size(τ̇Aer), size(τAer)
    #@show size((1 .- fᵗ * ω̃)' .* τ̇Aer), size(tmp .* τAer')
    τ̇_mod[2:5,:] .= (1 .- fᵗ * ω̃)' .* τ̇Aer[2:5,:] .- tmp .* τAer'
    #@show size(ω̃̇ * (1 - fᵗ)), size(ḟᵗ * (ω̃ .* (1 .- ω̃))'), size((1 .- fᵗ * ω̃)'.^2)
    ϖ̇_mod[2:5,:] .= (ω̃̇[2:5] * (1 - fᵗ) .- ḟᵗ[2:5] .* (ω̃ .* (1 .- ω̃))') ./ (1 .- fᵗ * ω̃)'.^2

    τ̇_mod[6:7,:] .= (1 .- fᵗ * ω̃)' .* τ̇Aer[6:7,:]
    #@show size(ω̃̇ * (1 - fᵗ)), size(ḟᵗ * (ω̃ .* (1 .- ω̃))'), size((1 .- fᵗ * ω̃)'.^2)
    ϖ̇_mod[6:7,:] .= 0.0
    
    return CoreScatteringOpticalProperties(τ_mod, ϖ_mod, AerZ⁺⁺, AerZ⁻⁺), 
        CoreScatteringOpticalPropertiesLin(τ̇_mod, ϖ̇_mod, AerŻ⁺⁺, AerŻ⁻⁺)
end


# Extract scattering definitions and integrated absorptions for the source function!
"""
    extractEffectiveProps(lods, lods_lin)

Extract effective scattering properties from combined layer optical property arrays.

Computes the cumulative optical depth sums and determines scattering interfaces for
each layer. The scattering interface type (`ScatteringInterface_00`, `01`, `10`, `11`)
controls which adding method variant is dispatched during the interaction step.

# Returns
- `scattering_interfaces_all`: Array of `ScatteringInterface` types per layer.
- `τ_sum_all`: Cumulative optical depth `[nSpec × (Nz+1)]`.
- `τ̇_sum_all`: Cumulative τ derivative `[Nparams × nSpec × (Nz+1)]`.
"""
function extractEffectiveProps(
                lods::Array,#{CoreScatteringOpticalProperties{FT},1}
                lods_lin::Array) #where FT

    FT    = eltype(lods[1].τ)
    nSpec = length(lods[1].τ)
    nZ    = length(lods)
    nParams = size(lods_lin[1].τ̇)[1]
    # First the Scattering Interfaces:
    scattering_interface = ScatteringInterface_00()
    scattering_interfaces_all = []
    τ_sum_all = similar(lods[1].τ,(nSpec,nZ+1))
    τ_sum_all[:,1] .= 0
    τ̇_sum_all = similar(lods_lin[1].τ̇,(nParams,nSpec,nZ+1))
    τ̇_sum_all[:,:,1] .= 0
    #@show FT
    for iz =1:nZ
        # Need to check max entries in Z matrices here as well later!
        scatter = maximum(lods[iz].τ .* lods[iz].ϖ) > 2eps(FT)
        scattering_interface = get_scattering_interface(scattering_interface, scatter, iz)
        push!(scattering_interfaces_all, scattering_interface)
        @views τ_sum_all[:,iz+1] = τ_sum_all[:,iz] + lods[iz].τ 
        for ip=1:nParams
            τ̇_sum_all[ip,:,iz+1] = τ̇_sum_all[ip,:,iz] + lods_lin[iz].τ̇[ip,:] 
        end
    end
    return scattering_interfaces_all, τ_sum_all, τ̇_sum_all
end

"""
    expandOpticalProperties(in, in_lin, arr_type)

Expand optical properties and their derivatives to full spectral resolution by
replicating Z matrices along the spectral dimension if they are spectrally constant.

Ensures that `Z[nμ, nμ, nSpec]` and `Ż[Nparams, nμ, nμ, nSpec]` have matching
spectral dimensions with `τ` and `ϖ`.
"""
function expandOpticalProperties(in::CoreScatteringOpticalProperties, in_lin::CoreScatteringOpticalPropertiesLin,  arr_type)
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = in 
    @unpack τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺ = in_lin 
    @assert length(τ) == length(ϖ) "τ and ϖ sizes need to match"
    @assert length(τ̇) == length(ϖ̇) "τ̇ and ϖ̇ sizes need to match"

    #nParams = size(τ̇)[1]
    if size(Z⁺⁺,3) == 1
        Z⁺⁺ = _repeat(Z⁺⁺,1,1,length(τ))
        Z⁻⁺ = _repeat(Z⁻⁺,1,1,length(τ))
        Ż⁺⁺ = _repeat(Ż⁺⁺,1,1,1,length(τ))
        Ż⁻⁺ = _repeat(Ż⁻⁺,1,1,1,length(τ))
        return CoreScatteringOpticalProperties(arr_type(τ), arr_type(ϖ), arr_type(Z⁺⁺), arr_type(Z⁻⁺)), 
            CoreScatteringOpticalPropertiesLin(arr_type(τ̇), arr_type(ϖ̇), arr_type(Ż⁺⁺), arr_type(Ż⁻⁺))      
    else
        @assert size(Z⁺⁺,3) ==  length(τ) "Z and τ dimensions need to match "
        return CoreScatteringOpticalProperties(arr_type(τ), arr_type(ϖ), arr_type(Z⁺⁺), arr_type(Z⁻⁺)),
            CoreScatteringOpticalPropertiesLin(arr_type(τ̇), arr_type(ϖ̇), arr_type(Ż⁺⁺), arr_type(Ż⁻⁺))       
    end
end
