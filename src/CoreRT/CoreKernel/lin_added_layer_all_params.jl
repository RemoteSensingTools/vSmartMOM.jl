#=
 
This file contains RT interaction-related functions
 
=#

# No scattering in either the added layer or the composite layer
function lin_added_layer_all_params_helper!(RS_type::noRS{FT}, 
                                pol_type, SFI, quad_points, 
                                computed_layer_properties_lin, 
                                added_layer_lin::AddedLayerLin{FT},
                                architecture) where {FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    @unpack ap_ṙ⁺⁻, ap_ṙ⁻⁺, ap_ṫ⁻⁻, ap_ṫ⁺⁺, ap_J̇₀⁺, ap_J̇₀⁻ = added_layer_lin
    @unpack τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺ = computed_layer_properties_lin
    @unpack D, n = pol_type
    @unpack qp_μ, μ₀, Nquad, iμ₀Nstart = quad_points
    @unpack F₀ = RS_type

    arr_type = array_type(architecture)

    nD=Int(size(Ż⁺⁺,2)/n)
    D_diag = repeat(arr_type(D), nD)             # full diagonal entries
    bigD = Diagonal(D_diag)                     # D-matrix
    
    nparams = size(computed_layer_properties_lin.τ̇)[1]
    #ap_ṙ⁺⁻ = zeros(Nparams, size(added_layer_lin.ṙ⁺⁻)[2], size(added_layer_lin.ṙ⁺⁻)[3], size(added_layer_lin.ṙ⁺⁻)[4])
    #ap_ṙ⁻⁺ = zeros(Nparams, size(added_layer_lin.ṙ⁺⁻)[2], size(added_layer_lin.ṙ⁺⁻)[3], size(added_layer_lin.ṙ⁺⁻)[4])
    #ap_ṫ⁺⁺ = zeros(Nparams, size(added_layer_lin.ṙ⁺⁻)[2], size(added_layer_lin.ṙ⁺⁻)[3], size(added_layer_lin.ṙ⁺⁻)[4])
    #ap_ṫ⁻⁻ = zeros(Nparams, size(added_layer_lin.ṙ⁺⁻)[2], size(added_layer_lin.ṙ⁺⁻)[3], size(added_layer_lin.ṙ⁺⁻)[4])
    #ap_J̇₀⁺ = zeros(Nparams, size(added_layer_lin.J̇₀⁺)[2], size(added_layer_lin.J̇₀⁺)[3], size(added_layer_lin.J̇₀⁺)[4])
    #ap_J̇₀⁻ = zeros(Nparams, size(added_layer_lin.J̇₀⁻)[2], size(added_layer_lin.J̇₀⁻)[3], size(added_layer_lin.J̇₀⁻)[4])   
    nspec = size(computed_layer_properties_lin.τ̇)[2]
    nbigD = size(bigD,1)
    #@show nD, n, nbigD
    i₀ = iμ₀Nstart:iμ₀Nstart+n-1
    #@show i₀
    Ż⁺⁺_I₀ = arr_type(zeros(nbigD, nspec))
    Ż⁻⁺_I₀ = arr_type(zeros(nbigD, nspec))
    Ż⁺⁺ = arr_type(Ż⁺⁺)
    Ż⁻⁺ = arr_type(Ż⁻⁺)
    #@show size(Ż⁺⁺), size(Ż⁺⁺_I₀), size(F₀)
    F₀ = arr_type(F₀)
    for iparam=1:nparams 
        # the following is placeholder code: check later for 
        # 1. use of dτ̇_λ/dϖ̇_λ vs. dτ̇/dϖ̇
        # 2. dimensions
        for ii = 1:nspec
            Ż⁺⁺_I₀[:,ii] = Ż⁺⁺[iparam,:,i₀,ii] * F₀[:,ii] #I₀[ii-i_start+1]
            Ż⁻⁺_I₀[:,ii] = Ż⁻⁺[iparam,:,i₀,ii] * F₀[:,ii] #I₀[ii-i_start+1] 
        end
        @views ap_ṫ⁺⁺[iparam,:,:,:] .= added_layer_lin.ṫ⁺⁺[1,:,:,:].*reshape(τ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṫ⁺⁺[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṫ⁺⁺[3,:,:,:].*Ż⁺⁺[iparam,:,:,:] 
        @views ap_ṫ⁻⁻[iparam,:,:,:] .= added_layer_lin.ṫ⁻⁻[1,:,:,:].*reshape(τ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṫ⁻⁻[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṫ⁻⁻[3,:,:,:].*(reshape(bigD,nbigD,nbigD,1).*Ż⁺⁺[iparam,:,:,:].*reshape(bigD,nbigD,nbigD,1)) #Ż⁻⁻[iparam,:,:,:] 

        @views ap_ṙ⁻⁺[iparam,:,:,:] .= added_layer_lin.ṙ⁻⁺[1,:,:,:].*reshape(τ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṙ⁻⁺[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṙ⁻⁺[3,:,:,:].*Ż⁻⁺[iparam,:,:,:]  
        @views ap_ṙ⁺⁻[iparam,:,:,:] .= added_layer_lin.ṙ⁺⁻[1,:,:,:].*reshape(τ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṙ⁺⁻[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṙ⁺⁻[3,:,:,:].*(reshape(bigD,nbigD,nbigD,1).*Ż⁻⁺[iparam,:,:,:].*reshape(bigD,nbigD,nbigD,1)) #Ż⁺⁻[iparam,:,:,:] 
        if SFI
            @views ap_J̇₀⁺[iparam,:,1,:] .= added_layer_lin.J̇₀⁺[1,:,1,:].*reshape(τ̇[iparam,:],1,nspec) + 
                                added_layer_lin.J̇₀⁺[2,:,1,:].*reshape(ϖ̇[iparam,:],1,nspec) + 
                                added_layer_lin.J̇₀⁺[3,:,1,:].*Ż⁺⁺_I₀
            @views ap_J̇₀⁻[iparam,:,1,:] .= added_layer_lin.J̇₀⁻[1,:,1,:].*reshape(τ̇[iparam,:],1,nspec) + 
                                added_layer_lin.J̇₀⁻[2,:,1,:].*reshape(ϖ̇[iparam,:],1,nspec) + 
                                added_layer_lin.J̇₀⁻[3,:,1,:].*Ż⁻⁺_I₀ 
        end
    end
end

"Compute interaction between composite and added layers"
function lin_added_layer_all_params!(RS_type::noRS{FT}, 
    pol_type, SFI, quad_points,
    computed_layer_properties_lin, 
    added_layer_lin::AddedLayerLin{FT}, architecture) where {FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    
    #@show A1[1,1,1], A2[1,1,1]
    lin_added_layer_all_params_helper!(RS_type, pol_type, 
                    SFI, quad_points,
                    computed_layer_properties_lin, 
                    added_layer_lin, architecture)
    #A1 = Array(composite_layer.J₀⁻)
    #A2 = Array(composite_layer.J₀⁺)
    #@show A1[1,1,1], A2[1,1,1]
    synchronize_if_gpu()
end