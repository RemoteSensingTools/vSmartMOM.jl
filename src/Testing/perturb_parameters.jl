# defining vSmartMOM modes to be used in the following
fwd_mode = FwdMode() #vSmartMOM.CoreRT.Mode.FwdMode()
lin_mode = LinMode() #vSmartMOM.CoreRT.Mode.LinMode() 

function perturb_parameters(params::vSmartMOM.CoreRT.vSmartMOM_Parameters{Float64}, ppct::FT)
    pert_params = []
    incr_fct    = 1 + ppct/100
@show "hi1"
    # q
    tmp_params = params
    if all(params.q .== 0)
        tmp_params.q .= 0
    else
        tmp_params.q .*= incr_fct
    end
    push!(pert_params, tmp_params)
@show "hi2"
    # variable molecules
    Nbands = length(params.brdf)
    for iband=1:Nbands
        var_mols = params.absorption_params.variable_molecules[iband]
        Nmol = length(var_mols)

        for imol=1:Nmol
            tmp_params = params
            mol_str = var_mols[imol]
            if all(params.absorption_params.vmr[mol_str] .== 0)
                tmp_params.absorption_params.vmr[mol_str] .== 1e-9
            else
                tmp_params.absorption_params.vmr[mol_str] .*= incr_fct
            end
            push!(pert_params, tmp_params)
        end
    end
@show "hi3"
    # Aerosols
    Naer = length(params.scattering_params.rt_aerosols)
    for iaer=1:Naer
        # τ_ref
        tmp_params = params
        tmp_params.scattering_params.rt_aerosols[iaer].τ_ref == 0 ? 0.01 : incr_fct * tmp_params.scattering_params.rt_aerosols[iaer].τ_ref    
        push!(pert_params, tmp_params)
    
        # nᵣ
        tmp_params = params
        tmp_params.scattering_params.rt_aerosols[iaer].aerosol.nᵣ *= incr_fct
        push!(pert_params, tmp_params)   
        
        # nᵢ
        tmp_params = params
        tmp_params.scattering_params.rt_aerosols[iaer].aerosol.nᵢ == 0 ? 0.0001 : incr_fct*tmp_params.scattering_params.rt_aerosols[iaer].aerosol.nᵢ 
        push!(pert_params, tmp_params)    

        # log normal μ
        tmp_params = params
        μ = tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.μ
        σ = tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.σ
        tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution = 
            LogNormal(μ == 0 ? 0.01 : incr_fct*μ, σ)
        push!(pert_params, tmp_params)

        # log normal σ
        tmp_params = params
        μ = tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.μ
        σ = tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.σ
        tmp_params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution = 
            LogNormal(μ, incr_fct*σ)
        push!(pert_params, tmp_params)  

        # z₀
        tmp_params = params
        tmp_params.scattering_params.rt_aerosols[iaer].z₀ *= incr_fct
        push!(pert_params, tmp_params)  

        # σ₀
        tmp_params = params
        tmp_params.scattering_params.rt_aerosols[iaer].σ₀ *= incr_fct
        push!(pert_params, tmp_params)  
    end
@show "hi4"
    # BRDF
    for iband=1:Nbands
        tmp_params = params
        tmp_params.brdf[iband].albedo == 0 ? 0.01 : incr_fct*tmp_params.brdf[iband].albedo
        push!(pert_params, tmp_params)
    end
@show "hi5"
    return pert_params
end

function compute_FD_modelJacobian(params::vSmartMOM_Parameters, ppct::FT)
    model, lin_model = model_from_parameters(lin_mode, params)

    pert_params = perturb_parameters(params, ppct)

    Nparams = length(pert_params)
    var_mols = params.absorption_params.variable_molecules
    Nmol = length(var_mols)
    Naer = length(params.scattering_params.rt_aerosols)
    Nbands = length(params.brdf)

    @assert Nparams == 1 + Nmol + Naer * 7 + Nbands
    
    # q
    Δq = pert_params[1].q - params.q
    pert_model = model_from_parameters(fwd_mode, pert_params[1])
    K_FD = (pert_model.τ_abs .- model.τ_abs)./Δq
    assert K_FD ≈ lin_model.τ̇_abs[1]

    # variable_molecules
    for imol=1:Nmol
        Δmol = pert_params[1+imol].absorption_params.vmr[var_mols[imol]] - params.absorption_params.vmr[var_mols[imol]]
        pert_model = model_from_parameters(fwd_mode, pert_params[1+imol])
        K_FD = (pert_model.τ_abs .- model.τ_abs)./Δmol
        assert K_FD ≈ lin_model.τ̇_abs[1+imol]
    end

    #Aerosols
    for iaer = 1:Naer
        #τ_ref
        idx = 1 + Nmol + 7*(iaer-1) + 1
        Δτ_ref = pert_params.scattering_params.rt_aerosols[iaer].τ_ref - params.scattering_params.rt_aerosols[iaer].τ_ref
        pert_model = model_from_parameters(fwd_mode, pert_params[idx])
        K_FD = (pert_model.τ_aer .- model.τ_aer)./Δτ_ref
        assert K_FD ≈ lin_model.τ̇_aer[idx]

        # nᵣ
        idx = 1 + Nmol + 7*(iaer-1) + 2
        Δnᵣ = pert_params.scattering_params.rt_aerosols[iaer].aerosol.nᵣ - params.scattering_params.rt_aerosols[iaer].aerosol.nᵢ
        pert_model = model_from_parameters(fwd_mode, pert_params[idx])
        K_FD = (pert_model.τ_aer .- model.τ_aer)./Δnᵣ
        assert K_FD ≈ lin_model.τ̇_aer[idx]
        
        for iband=1:Nband
            # k [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].k - model.aerosol_optics[iband][iaer].k)/Δnᵣ
            assert K_FD ≈ lin_model.lin_aerosol_optics[iband][iaer].k̇[1,:]
    
            # ω̃ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].ω̃ - model.aerosol_optics[iband][iaer].ω̃)/Δnᵣ
            assert K_FD ≈ lin_model.lin_aerosol_optics[iband][iaer].ω̃̇[1,:]

            # fᵗ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].fᵗ - model.aerosol_optics[iband][iaer].fᵗ)/Δnᵣ
            assert K_FD ≈ lin_model.lin_aerosol_optics[iband][iaer].ḟᵗ[1,:]

            # greek_coefs.α [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.α - model.aerosol_optics[iband][iaer].greek_coefs.α)/Δnᵣ
            assert K_FD ≈ lin_model.lin_aerosol_optics[iband][iaer].greek_coefs.α̇[1,:]

            # greek_coefs.β [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.β - model.aerosol_optics[iband][iaer].greek_coefs.β)/Δnᵣ
            assert K_FD ≈ lin_model.lin_aerosol_optics[iband][iaer].greek_coefs.β̇[1,:]

            # greek_coefs.γ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.γ - model.aerosol_optics[iband][iaer].greek_coefs.γ)/Δnᵣ
            assert K_FD ≈ lin_model.lin_aerosol_optics[iband][iaer].greek_coefs.γ̇[1,:]

            # greek_coefs.δ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.δ - model.aerosol_optics[iband][iaer].greek_coefs.δ)/Δnᵣ
            assert K_FD ≈ lin_model.lin_aerosol_optics[iband][iaer].greek_coefs.δ̇[1,:]

            # greek_coefs.ϵ [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.ϵ - model.aerosol_optics[iband][iaer].greek_coefs.ϵ)/Δnᵣ
            assert K_FD ≈ lin_model.lin_aerosol_optics[iband][iaer].greek_coefs.ϵ̇[1,:]

            # greek_coefs.α [iband][iaer]
            K_FD = (pert_model.aerosol_optics[iband][iaer].greek_coefs.ζ - model.aerosol_optics[iband][iaer].greek_coefs.ζ)/Δnᵣ
            assert K_FD ≈ lin_model.lin_aerosol_optics[iband][iaer].greek_coefs.ζ̇[1,:]
        end
    end
end