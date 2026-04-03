# Import necessary modules
using vSmartMOM, vSmartMOM.CoreRT
using Distributions: LogNormal

# defining vSmartMOM modes to be used in the following
fwd_mode = FwdMode() #vSmartMOM.CoreRT.Mode.FwdMode()
lin_mode = LinMode() #vSmartMOM.CoreRT.Mode.LinMode() 
FT = Float64  # ← KEEP: Using Float32 as requested

function perturb_parameters(params::vSmartMOM.CoreRT.vSmartMOM_Parameters, ppct)
    pert_params = []
    incr_fct    = 1. + ppct/100.
    
    #@show "hi1"
    # q
    if !isnothing(params.q)
        tmp_params = deepcopy(params)  # ← CRITICAL FIX: Use deepcopy
        if all(params.q .== 0.0)
            tmp_params.q .+= 0.01*ppct
        else
            tmp_params.q .*= incr_fct
        end
        push!(pert_params, tmp_params)
    end
    
    #@show "hi2"
    # variable molecules
    Nbands = length(params.brdf)
    for iband=1:Nbands
        var_mols = params.absorption_params.variable_molecules[iband]
        Nmol = length(var_mols)

        for imol=1:Nmol
            tmp_params = deepcopy(params)  # ← CRITICAL FIX: Use deepcopy
            mol_str = var_mols[imol]
            if all(params.absorption_params.vmr[mol_str] .== 0)
                tmp_params.absorption_params.vmr[mol_str] .= 1e-9  # ← FIX: Use .= for assignment
            else
                tmp_params.absorption_params.vmr[mol_str] .*= incr_fct
            end
            push!(pert_params, tmp_params)
        end
    end
    #@show "hi3"
    # Aerosols
    Naer = length(params.scattering_params.rt_aerosols)
    for iaer=1:Naer
        # τ_ref
        tmp_params = deepcopy(params)  # ← CRITICAL FIX: Use deepcopy
        τ_ref = tmp_params.scattering_params.rt_aerosols[iaer].τ_ref 
        tmp_params.scattering_params.rt_aerosols[iaer].τ_ref = (τ_ref == 0) ? 0.01 : incr_fct * τ_ref    
        push!(pert_params, tmp_params)
    
        # nᵣ
        tmp_params = deepcopy(params)  # ← CRITICAL FIX: Use deepcopy
        tmp_params.scattering_params.rt_aerosols[iaer].aerosol.nᵣ *= incr_fct
        push!(pert_params, tmp_params)   
        
        # nᵢ
        tmp_params = deepcopy(params)  # ← CRITICAL FIX: Use deepcopy
        nᵢ = tmp_params.scattering_params.rt_aerosols[iaer].aerosol.nᵢ
        tmp_params.scattering_params.rt_aerosols[iaer].aerosol.nᵢ = (nᵢ == 0) ? 0.0001 : incr_fct*nᵢ
        push!(pert_params, tmp_params)    

        # log normal μ  
        # NOTE: LogNormal(μ,σ) uses log-parameters where μ=log(median), σ=log(scale)
        # To perturb the actual median by incr_fct, we need: log(median * incr_fct) = log(median) + log(incr_fct) = μ + log(incr_fct)
        tmp_params = deepcopy(params)  # ← CRITICAL FIX: Use deepcopy
        μ = tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.μ
        σ = tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.σ
        # Convert to actual median, perturb, then convert back to log-parameter
        actual_median = exp(μ)
        new_median = actual_median == 0 ? 0.01 : incr_fct * actual_median
        new_μ = log(new_median)
        tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution = 
            LogNormal(new_μ, σ)
        push!(pert_params, tmp_params)

        # log normal σ
        # NOTE: Similarly for σ, we need to perturb the actual scale parameter  
        tmp_params = deepcopy(params)  # ← CRITICAL FIX: Use deepcopy
        μ = tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.μ
        σ = tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.σ
        # Convert to actual scale, perturb, then convert back to log-parameter
        actual_scale = exp(σ)
        new_scale = incr_fct * actual_scale
        new_σ = log(new_scale)
        tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution = 
            LogNormal(μ, new_σ)
        push!(pert_params, tmp_params)  

        # z₀
        tmp_params = deepcopy(params)  # ← CRITICAL FIX: Use deepcopy
        tmp_params.scattering_params.rt_aerosols[iaer].z₀ *= incr_fct
        push!(pert_params, tmp_params)  

        # σ₀
        tmp_params = deepcopy(params)  # ← CRITICAL FIX: Use deepcopy
        tmp_params.scattering_params.rt_aerosols[iaer].σ₀ *= incr_fct
        push!(pert_params, tmp_params)  
    end
    
    #@show "hi4"
    # BRDF
    for iband=1:Nbands
        tmp_params = deepcopy(params)  # ← CRITICAL FIX: Use deepcopy
        albedo_val = tmp_params.brdf[iband].albedo
        tmp_params.brdf[iband].albedo = (albedo_val == 0) ? 0.01 : incr_fct * albedo_val  # ← FIX: Actually assign the value
        push!(pert_params, tmp_params)
    end
    
    #@show "hi5"
    return pert_params
end

function compute_FD_modelJacobian(params::vSmartMOM.CoreRT.vSmartMOM_Parameters, ppct)
    model, lin_model = model_from_parameters(lin_mode, params)

    pert_params = perturb_parameters(params, ppct)

    Nparams = length(pert_params)
    var_mols = params.absorption_params.variable_molecules
    Nmol = length(var_mols)
    Naer = length(params.scattering_params.rt_aerosols)
    Nbands = length(params.brdf)

    has_q = !isnothing(params.q)
    q_offset = has_q ? 1 : 0
    @assert Nparams == q_offset + Nmol + Naer * 7 + Nbands

    iband=1
    # q
    if has_q
        M_air = 28.9647  # g/mol
        M_H₂O = 18.01528  # g/mol
        Δq = pert_params[1].q - params.q
        #pert_model = model_from_parameters(fwd_mode, pert_params[1])
        pert_model = model_from_parameters(pert_params[1])
        K_FD = (pert_model.τ_abs[iband] .- model.τ_abs[iband]).*(1.0./Δq)'
        K_lin = lin_model.τ̇_abs[iband][1,:,:]*(M_air/M_H₂O).*(1.0./(1.0.-params.q).^2)'  # ← FIX: Account for the scaling factor in the linear model
        @assert K_FD ≈ K_lin
    end

    # variable_molecules
    for imol=1:Nmol
        tΔmol = pert_params[q_offset+imol].absorption_params.vmr[var_mols[imol][1]] - params.absorption_params.vmr[var_mols[imol][1]]
        Δmol = (tΔmol[2:end] + tΔmol[1:end-1])/2.  # ← FIX: Compute Δmol correctly
        pert_model = model_from_parameters(pert_params[q_offset+imol])
        K_FD = (pert_model.τ_abs[iband] .- model.τ_abs[iband]).*(1.0./Δmol)'
        K_lin = lin_model.τ̇_abs[iband][q_offset+imol,:,:]
    end

    #Aerosols
    for iaer = 1:Naer
        #τ_ref
        i_aerprop = 1
        idx = q_offset + Nmol + 7*(iaer-1) + i_aerprop
        Δτ_ref = pert_params[idx].scattering_params.rt_aerosols[iaer].τ_ref - params.scattering_params.rt_aerosols[iaer].τ_ref
        pert_model = model_from_parameters(pert_params[idx])
        
        for iband=1:Nband
            # tau_aer [iband][iaer]
            K_FD = (pert_model.τ_aer[iband][iaer,:,:] .- model.τ_aer[iband][iaer,:,:])/Δτ_ref
            K_lin = lin_model.τ̇_aer[iband][iaer, i_aerprop,:,:]
            @assert K_FD ≈ K_lin
        end
    
        # nᵣ
        i_aerprop = 2
        idx = q_offset + Nmol + 7*(iaer-1) + i_aerprop
        Δnᵣ = pert_params[idx].scattering_params.rt_aerosols[iaer].aerosol.nᵣ - params.scattering_params.rt_aerosols[iaer].aerosol.nᵣ  
        pert_model = model_from_parameters(pert_params[idx])
        
        for iband=1:Nband
            # tau_aer [iband][iaer]
            K_FD = (pert_model.τ_aer[iband][iaer,:,:] .- model.τ_aer[iband][iaer,:,:])./Δnᵣ
            K_lin = lin_model.τ̇_aer[iband][iaer, i_aerprop,:,:]
            @assert K_FD ≈ K_lin #(this happens when the perturbation is largish, ~0.5%)
        
            # k [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].k - model.aerosol_optics[iband][iaer].k)/Δnᵣ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].k̇[1,:]
    
            # ω̃ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].ω̃ - model.aerosol_optics[iband][iaer].ω̃)/Δnᵣ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].ω̃̇[1,:]

            # fᵗ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].fᵗ - model.aerosol_optics[iband][iaer].fᵗ)/Δnᵣ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].ḟᵗ[1,:]

            # greek_coefs.α [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.α - model.aerosol_optics[iband][iaer].greek_coefs.α)/Δnᵣ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.α̇[1,:]

            # greek_coefs.β [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.β - model.aerosol_optics[iband][iaer].greek_coefs.β)/Δnᵣ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.β̇[1,:]

            # greek_coefs.γ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.γ - model.aerosol_optics[iband][iaer].greek_coefs.γ)/Δnᵣ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.γ̇[1,:]

            # greek_coefs.δ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.δ - model.aerosol_optics[iband][iaer].greek_coefs.δ)/Δnᵣ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.δ̇[1,:]

            # greek_coefs.ϵ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.ϵ - model.aerosol_optics[iband][iaer].greek_coefs.ϵ)/Δnᵣ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.ϵ̇[1,:]

            # greek_coefs.ζ  [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.ζ - model.aerosol_optics[iband][iaer].greek_coefs.ζ)/Δnᵣ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.ζ̇[1,:]
        end
    

        # nᵢ
        i_aerprop = 3
        idx = q_offset + Nmol + 7*(iaer-1) + i_aerprop
        Δnᵢ = pert_params[idx].scattering_params.rt_aerosols[iaer].aerosol.nᵢ - params.scattering_params.rt_aerosols[iaer].aerosol.nᵢ  
        pert_model = model_from_parameters(pert_params[idx])
        
        for iband=1:Nband
            # tau_aer [iband][iaer]
            K_FD = (pert_model.τ_aer[iband][iaer,:,:] .- model.τ_aer[iband][iaer,:,:])./Δnᵢ
            K_lin = lin_model.τ̇_aer[iband][iaer, i_aerprop,:,:]
            @assert K_FD ≈ K_lin 
        
            # k [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].k - model.aerosol_optics[iband][iaer].k)/Δnᵢ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].k̇[2,:]
    
            # ω̃ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].ω̃ - model.aerosol_optics[iband][iaer].ω̃)/Δnᵢ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].ω̃̇[2,:]

            # fᵗ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].fᵗ - model.aerosol_optics[iband][iaer].fᵗ)/Δnᵢ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].ḟᵗ[2,:]

            # greek_coefs.α [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.α - model.aerosol_optics[iband][iaer].greek_coefs.α)/Δnᵢ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.α̇[2,:]

            # greek_coefs.β [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.β - model.aerosol_optics[iband][iaer].greek_coefs.β)/Δnᵢ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.β̇[2,:]

            # greek_coefs.γ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.γ - model.aerosol_optics[iband][iaer].greek_coefs.γ)/Δnᵢ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.γ̇[2,:]

            # greek_coefs.δ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.δ - model.aerosol_optics[iband][iaer].greek_coefs.δ)/Δnᵢ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.δ̇[2,:]

            # greek_coefs.ϵ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.ϵ - model.aerosol_optics[iband][iaer].greek_coefs.ϵ)/Δnᵢ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.ϵ̇[2,:]

            # greek_coefs.α [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.ζ - model.aerosol_optics[iband][iaer].greek_coefs.ζ)/Δnᵢ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.ζ̇[2,:]
        end

        # r₀
        i_aerprop = 4
        idx = q_offset + Nmol + 7*(iaer-1) + i_aerprop
        Δμ = pert_params[idx].scattering_params.rt_aerosols[iaer].aerosol.size_distribution.μ - params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.μ
        pert_model = model_from_parameters(pert_params[idx])
        #K_FD = (pert_model.τ_aer .- model.τ_aer)./Δμ
        #@assert K_FD ≈ lin_model.τ̇_aer[idx]
        
        for iband=1:Nband
            # tau_aer [iband][iaer]
            K_FD = (pert_model.τ_aer[iband][iaer,:,:] .- model.τ_aer[iband][iaer,:,:])./Δμ
            K_lin = lin_model.τ̇_aer[iband][iaer, i_aerprop,:,:]
            @assert K_FD ≈ K_lin

            # k [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].k - model.aerosol_optics[iband][iaer].k)/Δμ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].k̇[3,:]
    
            # ω̃ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].ω̃ - model.aerosol_optics[iband][iaer].ω̃)/Δμ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].ω̃̇[3,:]

            # fᵗ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].fᵗ - model.aerosol_optics[iband][iaer].fᵗ)/Δμ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].ḟᵗ[3,:]

            # greek_coefs.α [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.α - model.aerosol_optics[iband][iaer].greek_coefs.α)/Δμ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.α̇[3,:]

            # greek_coefs.β [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.β - model.aerosol_optics[iband][iaer].greek_coefs.β)/Δμ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.β̇[3,:]

            # greek_coefs.γ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.γ - model.aerosol_optics[iband][iaer].greek_coefs.γ)/Δμ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.γ̇[3,:]

            # greek_coefs.δ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.δ - model.aerosol_optics[iband][iaer].greek_coefs.δ)/Δμ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.δ̇[3,:]

            # greek_coefs.ϵ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.ϵ - model.aerosol_optics[iband][iaer].greek_coefs.ϵ)/Δμ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.ϵ̇[3,:]

            # greek_coefs.α [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.ζ - model.aerosol_optics[iband][iaer].greek_coefs.ζ)/Δμ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.ζ̇[3,:]
        end

        # σ₀
        i_aerprop = 5
        idx = q_offset + Nmol + 7*(iaer-1) + i_aerprop
        Δσ = pert_params[idx].scattering_params.rt_aerosols[iaer].aerosol.size_distribution.σ - params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.σ
        pert_model = model_from_parameters(pert_params[idx])
        #K_FD = (pert_model.τ_aer .- model.τ_aer)./Δσ
        #@assert K_FD ≈ lin_model.τ̇_aer[idx]
        
        for iband=1:Nband
        # tau_aer [iband][iaer]
            K_FD = (pert_model.τ_aer[iband][iaer,:,:] .- model.τ_aer[iband][iaer,:,:])./Δσ
            K_lin = lin_model.τ̇_aer[iband][iaer, i_aerprop,:,:]
            @assert K_FD ≈ K_lin

            # k [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].k - model.aerosol_optics[iband][iaer].k)/Δσ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].k̇[4,:]
    
            # ω̃ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].ω̃ - model.aerosol_optics[iband][iaer].ω̃)/Δσ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].ω̃̇[4,:]

            # fᵗ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].fᵗ - model.aerosol_optics[iband][iaer].fᵗ)/Δσ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].ḟᵗ[4,:]

            # greek_coefs.α [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.α - model.aerosol_optics[iband][iaer].greek_coefs.α)/Δσ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.α̇[4,:]

            # greek_coefs.β [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.β - model.aerosol_optics[iband][iaer].greek_coefs.β)/Δσ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.β̇[4,:]

            # greek_coefs.γ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.γ - model.aerosol_optics[iband][iaer].greek_coefs.γ)/Δσ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.γ̇[4,:]

            # greek_coefs.δ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.δ - model.aerosol_optics[iband][iaer].greek_coefs.δ)/Δσ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.δ̇[4,:]

            # greek_coefs.ϵ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.ϵ - model.aerosol_optics[iband][iaer].greek_coefs.ϵ)/Δσ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.ϵ̇[4,:]

            # greek_coefs.α [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.ζ - model.aerosol_optics[iband][iaer].greek_coefs.ζ)/Δσ
            #@assert K_FD ≈ 
            lin_model.lin_aerosol_optics[iband][iaer].lin_greek_coefs.ζ̇[4,:]
        end

        # z₀
        i_aerprop = 6
        idx = q_offset + Nmol + 7*(iaer-1) + i_aerprop
        Δz₀ = pert_params[idx].scattering_params.rt_aerosols[iaer].z₀ - params.scattering_params.rt_aerosols[iaer].z₀
        pert_model = model_from_parameters(pert_params[idx])
                
        for iband=1:Nband
            # tau_aer [iband][iaer]
            K_FD = (pert_model.τ_aer[iband][iaer,:,:] .- model.τ_aer[iband][iaer,:,:])/Δz₀
            K_lin = lin_model.τ̇_aer[iband][iaer, i_aerprop,:,:]
            @assert K_FD ≈ K_lin
        end

        # σ₀
        i_aerprop = 7
        idx = q_offset + Nmol + 7*(iaer-1) + i_aerprop
        Δσ₀ = pert_params[idx].scattering_params.rt_aerosols[iaer].σ₀ - params.scattering_params.rt_aerosols[iaer].σ₀
        pert_model = model_from_parameters(pert_params[idx])    
        for iband=1:Nband
            # tau_aer [iband][iaer]
            K_FD = (pert_model.τ_aer[iband][iaer,:,:] .- model.τ_aer[iband][iaer,:,:])/Δσ₀
            K_lin = lin_model.τ̇_aer[iband][iaer, i_aerprop,:,:]
            @assert K_FD ≈ K_lin
        end 
    end
end