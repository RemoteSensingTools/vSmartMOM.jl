#=
 
This file implements rt_kernel!, which performs the core RT routines (elemental, doubling, interaction)
 
=#
#No Raman (default)
# Perform the Core RT routines (elemental, doubling, interaction)
function rt_kernel!(RS_type::noRS, pol_type, SFI, 
        added_layer, added_layer_lin, 
        composite_layer, composite_layer_lin, 
        computed_layer_properties, computed_layer_properties_lin, 
        m, quad_points, I_static, architecture, qp_μN, iz) 

    @unpack τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface = computed_layer_properties
    # Downselect the following parameters as appropriate
    @unpack τ̇_λ, ϖ̇_λ, τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺, dτ̇, dτ̇_λ, expk_lin, τ̇_sum = computed_layer_properties_lin
    @unpack F₀ = RS_type

    Nparams = size(τ̇_λ)[1]
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
            composite_layer_lin.Ṫ⁺⁺[iparam,:] += added_layer_lin.ṫ⁺⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁺⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁺⁺[3,:,:,:].*dŻ⁺⁺[iparam] 
            composite_layer_lin.Ṫ⁻⁻[iparam,:] += added_layer_lin.ṫ⁻⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁻⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁻⁻[3,:,:,:].*dŻ⁻⁻[iparam] 

            composite_layer_lin.Ṙ⁻⁺[iparam,:] += added_layer_lin.ṙ⁻⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁻⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁻⁺[3,:,:,:].*dŻ⁻⁺[iparam]  
            composite_layer_lin.Ṙ⁺⁻[iparam,:] += added_layer_lin.ṙ⁺⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁺⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁺⁻[3,:,:,:].*dŻ⁺⁻[iparam] 

            composite_layer_lin.J̇₀⁺[iparam,:] += added_layer_lin.J̇₀⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁺[3,:,:,:].*dŻ⁺⁺[iparam] 
            composite_layer_lin.J̇₀⁻[iparam,:] += added_layer_lin.J̇₀⁻[1,:,:,:].*dτ̇_λ[iparam] + 
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
    @unpack qp_μ, μ₀ = quad_points
    @unpack F₀ = RS_type
    # Just unpack core optical properties from 
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
    @unpack τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺ = computed_layer_properties_lin
    
    # SUNITI, check? Also, better to write function here
    #@show τ, ϖ
    #@show maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ) #τ, ϖ
    dτ_max = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])

    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    # @show ndoubl
    scatter = true # edit later
    arr_type = array_type(architecture)
    # Compute dτ vector
    dτ = τ ./ 2^ndoubl
    dτ̇ = τ̇ ./ 2^ndoubl
    expk = arr_type(exp.(-dτ /μ₀))
    expk_lin = arr_type(exp.(-dτ /μ₀)*(-1/μ₀))
    #@show dτ, ndoubl
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show F₀
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        τ_sum, dτ, F₀,
                                        #τ̇_sum, dτ̇,
                                        computed_layer_properties,
                                        #computed_layer_properties_lin, 
                                        m, ndoubl, scatter, quad_points,  
                                        added_layer,  
                                        added_layer_lin,
                                        architecture)
        #@show "Done"  
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
        #println("Elemental done...")
        @timeit "doubling"   doubling!(pol_type, SFI, 
                                        expk, expk_lin, 
                                        ndoubl, 
                                        added_layer, added_layer_lin, 
                                        I_static, architecture)
                
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
        temp = Array(exp.(-τ./qp_μN'))
        added_layer_lin.ṙ⁻⁺[:] .= 0;
        added_layer_lin.ṙ⁺⁻[:] .= 0;
        added_layer_lin.J̇₀⁻[:] .= 0;
        temp_lin = Array(exp.(-τ./qp_μN') * (-1 ./ qp_μN))
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ_λ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);

            added_layer_lin.ṫ⁺⁺[1,:,:,iλ] = Diagonal(temp_lin[iλ,:])
            added_layer_lin.ṫ⁻⁻[1,:,:,iλ] = Diagonal(temp_lin[iλ,:])
            
        end
    end

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
            composite_layer_lin.Ṫ⁺⁺[iparam,:] += added_layer_lin.ṫ⁺⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁺⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁺⁺[3,:,:,:].*dŻ⁺⁺[iparam] 
            composite_layer_lin.Ṫ⁻⁻[iparam,:] += added_layer_lin.ṫ⁻⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁻⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁻⁻[3,:,:,:].*dŻ⁻⁻[iparam] 

            composite_layer_lin.Ṙ⁻⁺[iparam,:] += added_layer_lin.ṙ⁻⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁻⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁻⁺[3,:,:,:].*dŻ⁻⁺[iparam]  
            composite_layer_lin.Ṙ⁺⁻[iparam,:] += added_layer_lin.ṙ⁺⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁺⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁺⁻[3,:,:,:].*dŻ⁺⁻[iparam] 

            composite_layer_lin.J̇₀⁺[iparam,:] += added_layer_lin.J̇₀⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁺[3,:,:,:].*dŻ⁺⁺[iparam] 
            composite_layer_lin.J̇₀⁻[iparam,:] += added_layer_lin.J̇₀⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁻[3,:,:,:].*dŻ⁻⁺[iparam] 
        end
    # If this is not the TOA, perform the interaction step
    else
        @timeit "lin_added_layer_all_params" lin_added_layer_all_params!(SFI, 
                    computed_layer_properties_lin, 
                    added_layer_lin)
        @timeit "interaction" interaction!(RS_type, 
                    scattering_interface, 
                    SFI, 
                    computed_layer_properties, computed_layer_properties_lin, 
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


