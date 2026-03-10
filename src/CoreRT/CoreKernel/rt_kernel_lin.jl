#=
 
This file implements rt_kernel!, which performs the core RT routines (elemental, doubling, interaction)
 
=#
#No Raman (default)
# Perform the Core RT routines (elemental, doubling, interaction)
#=
function rt_kernel!(RS_type::noRS, pol_type, SFI, 
        added_layer, added_layer_lin, 
        composite_layer, composite_layer_lin, 
        computed_layer_properties, computed_layer_properties_lin, 
        m, quad_points, I_static, architecture, qp_μN, iz) 

    @unpack τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface = computed_layer_properties
    # Downselect the following parameters as appropriate
    @unpack τ̇_λ, ϖ̇_λ, τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺, dτ̇, dτ̇_λ, expk_lin, τ̇_sum = computed_layer_properties_lin
    @unpack F₀ = RS_type

    Nparams = size(τ̇_λ, 2)
    #@show τ, ϖ, dτ_max, ndoubl
    # If there is scattering, perform the elemental and doubling steps
    if scatter  
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        τ_sum, dτ_λ, dτ, 
                                        ϖ_λ, ϖ, 
                                        Z⁺⁺, Z⁻⁺, 
                                        τ̇_sum, dτ̇_λ, dτ̇, 
                                        ϖ̇_λ, ϖ̇, 
                                        Ż⁺⁺, Ż⁻⁺, 
                                        F₀,
                                        m, ndoubl, 
                                        scatter, 
                                        quad_points,  
                                        added_layer,  
                                        added_layer_lin,  
                                        I_static, 
                                        architecture)
        #println("Elemental done...")
        @timeit "doubling"   doubling!(pol_type, SFI, expk, expk_lin, ndoubl, added_layer, added_layer_lin, I_static, architecture)
        #println("Doubling done...")
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.J₀⁻[:] .= 0;
        added_layer_lin.ṙ⁻⁺[:] .= 0;
        added_layer_lin.ṙ⁺⁻[:] .= 0;
        added_layer_lin.J̇₀⁻[:] .= 0;
        temp = Array(exp.(-τ_λ./qp_μN'))
        temp_lin = Array(exp.(-τ_λ./qp_μN') * (-1 ./ qp_μN))
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ_λ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
            
            # let ṫ, ṙ, snf J̇ in each layer be functions only of τ*, ϖ* and Z*, which in turn are composite functions of Nparams individual state parameters   
            added_layer_lin.ṫ⁺⁺[1,:,:,iλ] = Diagonal(temp_lin[iλ,:])
            added_layer_lin.ṫ⁻⁻[1,:,:,iλ] = Diagonal(temp_lin[iλ,:])
            
        end
    end
    #M1 = Array(added_layer.t⁺⁺)
    #M2 = Array(added_layer.r⁺⁻)
    #M3 = Array(added_layer.J₀⁻)
    #M4 = Array(added_layer.J₀⁺)
    #@show M1[1,1,1], M2[1,1,1], M3[1,1,1], M4[1,1,1]
    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        composite_layer.T⁺⁺[:], composite_layer.T⁻⁻[:] = (added_layer.t⁺⁺, added_layer.t⁻⁻)
        composite_layer.R⁻⁺[:], composite_layer.R⁺⁻[:] = (added_layer.r⁻⁺, added_layer.r⁺⁻)
        composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.J₀⁺, added_layer.J₀⁻ )
        
        # zero composite variables first
        composite_layer_lin.Ṫ⁺⁺[:] .= 0
        composite_layer_lin.Ṫ⁻⁻[:] .= 0
        composite_layer_lin.Ṙ⁻⁺[:] .= 0
        composite_layer_lin.Ṙ⁺⁻[:] .= 0
        composite_layer_lin.J̇₀⁺[:] .= 0
        composite_layer_lin.J̇₀⁻[:] .= 0

        for iparam = 1:Nparams
            # the following is placeholder code: check later for 
            # 1. use of dτ̇_λ/dϖ̇_λ vs. dτ̇/dϖ̇
            # 2. dimensions
            composite_layer_lin.Ṫ⁺⁺[iparam,:] += added_layer_lin.ṫ⁺⁺[:,:,:,1].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁺⁺[:,:,:,2].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁺⁺[:,:,:,3].*dŻ⁺⁺[iparam] 
            composite_layer_lin.Ṫ⁻⁻[iparam,:] += added_layer_lin.ṫ⁻⁻[:,:,:,1].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁻⁻[:,:,:,2].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁻⁻[:,:,:,3].*dŻ⁻⁻[iparam] 

            composite_layer_lin.Ṙ⁻⁺[iparam,:] += added_layer_lin.ṙ⁻⁺[:,:,:,1].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁻⁺[:,:,:,2].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁻⁺[:,:,:,3].*dŻ⁻⁺[iparam]  
            composite_layer_lin.Ṙ⁺⁻[iparam,:] += added_layer_lin.ṙ⁺⁻[:,:,:,1].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁺⁻[:,:,:,2].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁺⁻[:,:,:,3].*dŻ⁺⁻[iparam] 

            composite_layer_lin.J̇₀⁺[:,iparam] += added_layer_lin.J̇₀⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁺[3,:,:,:].*dŻ⁺⁺[iparam] 
            composite_layer_lin.J̇₀⁻[:,iparam] += added_layer_lin.J̇₀⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁻[3,:,:,:].*dŻ⁻⁺[iparam] 
        end
    # If this is not the TOA, perform the interaction step
    else
        @timeit "lin_added_layer_all_params" lin_added_layer_all_params!(SFI, 
                    computed_layer_properties_lin, 
                    added_layer_lin)
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, 
                            computed_layer_properties, computed_layer_properties_lin, 
                            composite_layer, composite_layer_lin, 
                            added_layer, added_layer_lin, 
                            I_static)
    end
end
=#

### New update: including towers/airborne sensors
function rt_kernel!(RS_type::noRS{FT}, 
                    pol_type, SFI, 
                    added_layer, added_layer_lin,
                    composite_layer, composite_layer_lin,
                    computed_layer_properties::CoreScatteringOpticalProperties, 
                    computed_layer_properties_lin::CoreScatteringOpticalPropertiesLin,
                    scattering_interface, 
                    τ_sum, τ̇_sum, 
                    m, quad_points, 
                    I_static, 
                    architecture, 
                    qp_μN, iz) where {FT}
    #@show array_type(architecture)
    @unpack qp_μ, μ₀, Nquad, iμ₀Nstart = quad_points
    @unpack F₀ = RS_type
    # Just unpack core optical properties from 
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
    @unpack τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺ = computed_layer_properties_lin
    @unpack D, n = pol_type

    arr_type = array_type(architecture)
    
    nD=Int(size(added_layer.t⁺⁺,1)/n)
    D_diag = repeat(arr_type(D), nD)             # full diagonal entries
    bigD = Diagonal(D_diag)                     # D-matrix
    # SUNITI, check? Also, better to write function here
    #@show τ, ϖ
    #@show maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ) #τ, ϖ
    dτ_max = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])

    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    # @show ndoubl
    scatter = scattering_interface==ScatteringInterface_10() || scattering_interface==ScatteringInterface_11() # edit later
    
    # Compute dτ vector
    dτ = τ ./ 2^ndoubl
    dτ̇ = τ̇ ./ 2^ndoubl

    expk = arr_type(exp.(-dτ /μ₀))
    expk_lin = arr_type(reshape(exp.(-dτ /μ₀)*(-1/μ₀), length(dτ), 1) .* dτ̇)
    Nparams = size(τ̇, 2)
#@show size(expk), size(expk_lin)
    #@show dτ, ndoubl
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show F₀
        #lin = LinMode()
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        dτ, arr_type(F₀),
                                        #τ̇_sum, dτ̇,
                                        computed_layer_properties,
                                        #computed_layer_properties_lin, 
                                        m, ndoubl, scatter, quad_points, 
                                        added_layer,  
                                        added_layer_lin,
                                        architecture)
        
        #@show 2  
        #=
        if m==0
            #m==0 ? 
            RayJ₀p = Array(added_layer.J₀⁺)
            RayJ₀m = Array(added_layer.J₀⁻)
            RayT   = Array(added_layer.t⁺⁺)
            RayR   = Array(added_layer.r⁻⁺)
            jldsave("/home/sanghavi/debugRay3.jld2"; RayJ₀p, RayJ₀m, RayT, RayR) 
        end                                
        =#
        
        #i_start  = pol_type.n*(quad_points.iμ₀-1) + 1 
        #AMF = FT(1/quad_points.μ₀) # AMF = 1/μ₀ 
        #i_end    = nStokes*iμ₀
        #println("Elemental done...")
        # Expanding derivatives to all parameters 
        @timeit "lin_added_layer_all_params" lin_added_layer_all_params!(
                    RS_type::noRS, pol_type,
                    SFI, quad_points, 
                    computed_layer_properties_lin, 
                    added_layer_lin, architecture)

        @timeit "doubling"   doubling!(pol_type, SFI, 
                                        expk, expk_lin, 
                                        arr_type(τ_sum), arr_type(τ̇_sum), 
                                        ndoubl, 
                                        #AMF,
                                        quad_points,
                                        added_layer, 
                                        added_layer_lin, 
                                        I_static, architecture)

        
                                        # Use the following between doubling and interaction steps to account for the fact that the added layer is not at the TOA. This is needed because the added layer is not at the TOA, so the derivatives of the added layer properties with respect to τ, ϖ and Z need to be scaled by exp(-τ_sum/μ₀) to account for the fact that the added layer is not at the TOA. This is done in the following lines of code.                    
        if SFI
            added_layer.J₀⁺[:, 1, :] .*= (exp.(-τ_sum[:]/μ₀))' #writing i_start:i_start to avoid scalar indexing errors with GPUArrays
            added_layer.J₀⁻[:, 1, :] .*= (exp.(-τ_sum[:]/μ₀))'

            for iparam=1:Nparams
                added_layer_lin.ap_J̇₀⁺[:, 1, :, iparam] = added_layer_lin.ap_J̇₀⁺[:, 1, :, iparam].*(exp.(-τ_sum[:]/μ₀))' +
                            added_layer.J₀⁺[:, 1, :] .* ((-1/μ₀) .* @view(τ̇_sum[:,iparam]))'
                added_layer_lin.ap_J̇₀⁻[:, 1, :, iparam] = added_layer_lin.ap_J̇₀⁻[:, 1, :, iparam].*(exp.(-τ_sum[:]/μ₀))' +
                            added_layer.J₀⁻[:, 1, :] .* ((-1/μ₀) .* @view(τ̇_sum[:,iparam]))'
            end
        ##    J̇₀⁺[2, :, 1, :] = J̇₀⁺[2, :, 1, :].*(exp.(-τ_sum[:]/μ₀))' 
        ##    J̇₀⁻[2, :, 1, :] = J̇₀⁻[2, :, 1, :].*(exp.(-τ_sum[:]/μ₀))' 
        ##    J̇₀⁺[3, :, 1, :] = J̇₀⁺[3, :, 1, :].*(exp.(-τ_sum[:]/μ₀))' 
        ##    J̇₀⁻[3, :, 1, :] = J̇₀⁻[3, :, 1, :].*(exp.(-τ_sum[:]/μ₀))' 
        end
        #@show 3        
                                        #=if m==0
            #m==0 ? 
            RayJ₀p = Array(added_layer.J₀⁺)
            RayJ₀m = Array(added_layer.J₀⁻)
            jldsave("/home/sanghavi/debugRay3.jld2"; RayJ₀p, RayJ₀m) 
        end=#
        
                                        #=if m==0
            #m==0 ? 
            RayJ₀ = Array(added_layer.J₀⁺)
            jldsave("/home/sanghavi/debugRay3.jld2"; RayJ₀) 
        end=#
        #println("Doubling done...")
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.J₀⁻[:] .= 0;
        temp = Array(exp.(-τ'./qp_μN))
        added_layer_lin.ap_ṙ⁻⁺[:] .= 0;
        added_layer_lin.ap_ṙ⁺⁻[:] .= 0;
        added_layer_lin.ap_J̇₀⁻[:] .= 0;
        temp_lin = Array(reshape(1, length(qp_μN), length(τ), exp.(-τ'./qp_μN)) .* reshape(1, length(qp_μN), 1, -1 ./ qp_μN) .* reshape(τ̇, 1, length(τ), Nparams)) 
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ_λ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[:,iλ]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[:,iλ]);
            
            for iparam=1:Nparams
                added_layer_lin.ṫ⁺⁺[:,:,iλ,iparam] = Diagonal(temp_lin[:, iλ, iparam])
                added_layer_lin.ṫ⁻⁻[:,:,iλ,iparam] = Diagonal(temp_lin[:, iλ, iparam])
            end
        end
    end
    #@show 4
    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        composite_layer.T⁺⁺ .= added_layer.t⁺⁺
        composite_layer.T⁻⁻ .= added_layer.t⁻⁻
        composite_layer.R⁻⁺ .= added_layer.r⁻⁺
        composite_layer.R⁺⁻ .= added_layer.r⁺⁻
        composite_layer.J₀⁺ .= added_layer.J₀⁺
        composite_layer.J₀⁻ .= added_layer.J₀⁻

        #composite_layer.T⁺⁺[:], composite_layer.T⁻⁻[:] = (added_layer.t⁺⁺, added_layer.t⁻⁻)
        #composite_layer.R⁻⁺[:], composite_layer.R⁺⁻[:] = (added_layer.r⁻⁺, added_layer.r⁺⁻)
        #composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.J₀⁺, added_layer.J₀⁻ )
    
        # zero composite variables first
        #=composite_layer_lin.Ṫ⁺⁺[:] .= 0
        composite_layer_lin.Ṫ⁻⁻[:] .= 0
        composite_layer_lin.Ṙ⁻⁺[:] .= 0
        composite_layer_lin.Ṙ⁺⁻[:] .= 0
        composite_layer_lin.J̇₀⁺[:] .= 0
        composite_layer_lin.J̇₀⁻[:] .= 0=#
        fill!(composite_layer_lin.Ṫ⁺⁺, 0)
        fill!(composite_layer_lin.Ṫ⁻⁻, 0)
        fill!(composite_layer_lin.Ṙ⁻⁺, 0)
        fill!(composite_layer_lin.Ṙ⁺⁻, 0)
        fill!(composite_layer_lin.J̇₀⁺, 0)
        fill!(composite_layer_lin.J̇₀⁻, 0)

        #nspec = size(added_layer.t⁺⁺,3)
        #nparams = size(τ̇)[1]
        #nbigD = size(bigD,1)
        #@show nD, n, nbigD
        #i₀ = iμ₀Nstart:iμ₀Nstart+n-1

        #Ż⁺⁺_I₀ = arr_type(zeros(nbigD, nspec))
        #Ż⁻⁺_I₀ = arr_type(zeros(nbigD, nspec))
        #Ż⁺⁺ = arr_type(Ż⁺⁺)
        #Ż⁻⁺ = arr_type(Ż⁻⁺)
        #F₀ = arr_type(F₀)
        for iparam = 1:Nparams
        
            @views composite_layer_lin.Ṫ⁺⁺[:,:,:,iparam] .= added_layer_lin.ap_ṫ⁺⁺[:,:,:,iparam]
            @views composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam] .= added_layer_lin.ap_ṫ⁻⁻[:,:,:,iparam]
            @views composite_layer_lin.Ṙ⁻⁺[:,:,:,iparam] .= added_layer_lin.ap_ṙ⁻⁺[:,:,:,iparam]
            @views composite_layer_lin.Ṙ⁺⁻[:,:,:,iparam] .= added_layer_lin.ap_ṙ⁺⁻[:,:,:,iparam]
            
            @views composite_layer_lin.J̇₀⁺[:,1,:,iparam] .= added_layer_lin.ap_J̇₀⁺[:,1,:,iparam]
            @views composite_layer_lin.J̇₀⁻[:,1,:,iparam] .= added_layer_lin.ap_J̇₀⁻[:,1,:,iparam]
            # the following is placeholder code: check later for 
            # 1. use of dτ̇_λ/dϖ̇_λ vs. dτ̇/dϖ̇
            # 2. dimensions
            #for ii = 1:nspec
            #    Ż⁺⁺_I₀[:,ii] = Ż⁺⁺[:,i₀,ii,iparam] * F₀[:,ii] #I₀[ii-i_start+1]
            #    Ż⁻⁺_I₀[:,ii] = Ż⁻⁺[:,i₀,ii,iparam] * F₀[:,ii] #I₀[ii-i_start+1] 
            #end
            #=
            @views composite_layer_lin.Ṫ⁺⁺[:,:,:,iparam] .= added_layer_lin.ṫ⁺⁺[:,:,:,1].*reshape(τ̇[:,iparam],1,1,nspec) .+ 
                                                added_layer_lin.ṫ⁺⁺[:,:,:,2].*reshape(ϖ̇[:,iparam],1,1,nspec) .+ 
                                                added_layer_lin.ṫ⁺⁺[:,:,:,3].*Ż⁺⁺[:,:,:,iparam] 
            @views composite_layer_lin.Ṫ⁻⁻[:,:,:,iparam] .= added_layer_lin.ṫ⁻⁻[:,:,:,1].*reshape(τ̇[:,iparam],1,1,nspec) .+ 
                                                added_layer_lin.ṫ⁻⁻[:,:,:,2].*reshape(ϖ̇[:,iparam],1,1,nspec) .+ 
                                                added_layer_lin.ṫ⁻⁻[:,:,:,3].*(reshape(bigD,nbigD,nbigD,1).*Ż⁺⁺[:,:,:,iparam].*reshape(bigD,nbigD,nbigD,1))

            @views composite_layer_lin.Ṙ⁻⁺[:,:,:,iparam] .= added_layer_lin.ṙ⁻⁺[:,:,:,1].*reshape(τ̇[:,iparam],1,1,nspec) .+ 
                                                added_layer_lin.ṙ⁻⁺[:,:,:,2].*reshape(ϖ̇[:,iparam],1,1,nspec) .+ 
                                                added_layer_lin.ṙ⁻⁺[:,:,:,3].*Ż⁻⁺[:,:,:,iparam]  
            @views composite_layer_lin.Ṙ⁺⁻[:,:,:,iparam] .= added_layer_lin.ṙ⁺⁻[:,:,:,1].*reshape(τ̇[:,iparam],1,1,nspec) .+ 
                                                added_layer_lin.ṙ⁺⁻[:,:,:,2].*reshape(ϖ̇[:,iparam],1,1,nspec) .+ 
                                                added_layer_lin.ṙ⁺⁻[:,:,:,3].*(reshape(bigD,nbigD,nbigD,1).*Ż⁻⁺[:,:,:,iparam].*reshape(bigD,nbigD,nbigD,1)) 

                                                #@show size(composite_layer_lin.J̇₀⁺), 
                                                #    size(added_layer_lin.J̇₀⁺), 
                                                #    size(dτ̇), size(ϖ̇), size(Ż⁺⁺)
            @views composite_layer_lin.J̇₀⁺[:,1,:,iparam] .= added_layer_lin.J̇₀⁺[:,1,:,1].*reshape(τ̇[:,iparam],1,nspec) .+ 
                                                added_layer_lin.J̇₀⁺[:,1,:,2].*reshape(ϖ̇[:,iparam],1,nspec) .+ 
                                                added_layer_lin.J̇₀⁺[:,1,:,3].*Ż⁺⁺_I₀
            @views composite_layer_lin.J̇₀⁻[:,1,:,iparam] .= added_layer_lin.J̇₀⁻[:,1,:,1].*reshape(τ̇[:,iparam],1,nspec) .+ 
                                                added_layer_lin.J̇₀⁻[:,1,:,2].*reshape(ϖ̇[:,iparam],1,nspec) .+ 
                                                added_layer_lin.J̇₀⁻[:,1,:,3].*Ż⁻⁺_I₀
            =#
        end
    # If this is not the TOA, perform the interaction step
    else
        
        #@timeit "lin_added_layer_all_params" lin_added_layer_all_params!(
        #            RS_type::noRS, pol_type,
        #            SFI, quad_points, 
        #            computed_layer_properties_lin, 
        #            added_layer_lin, architecture)

        @timeit "interaction" interaction!(scattering_interface, 
                    SFI, 
                    #computed_layer_properties, computed_layer_properties_lin, 
                    composite_layer, composite_layer_lin, 
                    added_layer, added_layer_lin, 
                    I_static)
        #=if iz==2
            M1 = Array(composite_layer.T⁺⁺);
            M2 = Array(composite_layer.R⁺⁻);
            M3 = Array(composite_layer.T⁻⁻);
            M4 = Array(composite_layer.R⁻⁺);
            M5 = Array(composite_layer.J₀⁻);
            M6 = Array(composite_layer.J₀⁺);
            #@show M1[1,1,1], M2[1,1,1], M3[1,1,1], M4[1,1,1], M5[1,1,1], M6[1,1,1]
        end=#
    end
    #=if m==0
        #m==0 ? 
        RayJ₀p = Array(composite_layer.J₀⁺)
        RayJ₀m = Array(composite_layer.J₀⁻)
        jldsave("/home/sanghavi/debugRay3.jld2"; RayJ₀p, RayJ₀m) 
    end=#
end


