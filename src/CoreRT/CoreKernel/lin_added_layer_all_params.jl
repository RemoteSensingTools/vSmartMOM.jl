#=
 
This file contains RT interaction-related functions
 
=#

# No scattering in either the added layer or the composite layer
function lin_added_layer_all_params_helper!(SFI,
                                computed_layer_properties_lin, 
                                added_layer_lin::AddedLayerLin{FT}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    @unpack ap_ṙ⁺⁻, ap_ṙ⁻⁺, ap_ṫ⁻⁻, ap_ṫ⁺⁺ = added_layer_lin
    @unpack τ̇_λ, ϖ̇_λ, τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺, dτ̇, dτ̇_λ, expk_lin, τ̇_sum = computed_layer_properties_lin

    Nparams = size(computed_layer_properties_lin.dτ̇_λ)[1]
    for iparam=1:Nparams 
        # the following is placeholder code: check later for 
        # 1. use of dτ̇_λ/dϖ̇_λ vs. dτ̇/dϖ̇
        # 2. dimensions
        ap_ṫ⁺⁺[iparam,:] += added_layer_lin.ṫ⁺⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                            added_layer_lin.ṫ⁺⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                            added_layer_lin.ṫ⁺⁺[3,:,:,:].*dŻ⁺⁺[iparam] 
        ap_ṫ⁻⁻[iparam,:] += added_layer_lin.ṫ⁻⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                            added_layer_lin.ṫ⁻⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                            added_layer_lin.ṫ⁻⁻[3,:,:,:].*dŻ⁻⁻[iparam] 

        ap_ṙ⁻⁺[iparam,:] += added_layer_lin.ṙ⁻⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                            added_layer_lin.ṙ⁻⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                            added_layer_lin.ṙ⁻⁺[3,:,:,:].*dŻ⁻⁺[iparam]  
        ap_ṙ⁺⁻[iparam,:] += added_layer_lin.ṙ⁺⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                            added_layer_lin.ṙ⁺⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                            added_layer_lin.ṙ⁺⁻[3,:,:,:].*dŻ⁺⁻[iparam] 
        if SFI
            ap_J̇₀⁺[iparam,:] += added_layer_lin.J̇₀⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                                added_layer_lin.J̇₀⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                                added_layer_lin.J̇₀⁺[3,:,:,:].*dŻ⁺⁺[iparam] 
            ap_J̇₀⁻[iparam,:] += added_layer_lin.J̇₀⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                                added_layer_lin.J̇₀⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                                added_layer_lin.J̇₀⁻[3,:,:,:].*dŻ⁻⁺[iparam] 
        end
    end
end

"Compute interaction between composite and added layers"
function lin_added_layer_all_params!(SFI,
    computed_layer_properties_lin, 
    added_layer_lin::AddedLayerLin{FT}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    
    #@show A1[1,1,1], A2[1,1,1]
    lin_added_layer_all_params_helper!(SFI,
                    computed_layer_properties_lin, 
                    added_layer_lin)
    #A1 = Array(composite_layer.J₀⁻)
    #A2 = Array(composite_layer.J₀⁺)
    #@show A1[1,1,1], A2[1,1,1]
    synchronize_if_gpu()
end