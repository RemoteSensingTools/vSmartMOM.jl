#=
 
This file implements rt_kernel!, which performs the core RT routines (elemental, doubling, interaction)
 
=#
#No Raman (default)
# Perform the Core RT routines (elemental, doubling, interaction)
function rt_kernel!(RS_type::noRS, 
                    pol_type, SFI, 
                    added_layer, 
                    composite_layer, 
                    computed_layer_properties, 
                    m, quad_points, 
                    I_static, architecture, 
                    qp_μN, iz) 

    @unpack τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface = computed_layer_properties
    @unpack F₀ = RS_type
    #@show τ, ϖ, dτ_max, ndoubl
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        τ_sum, dτ_λ, dτ, 
                                        ϖ_λ, ϖ, 
                                        Z⁺⁺, Z⁻⁺, 
                                        F₀,
                                        m, ndoubl, 
                                        scatter, 
                                        quad_points,  
                                        added_layer,  
                                        I_static, 
                                        architecture)
        #println("Elemental done...")
        @timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
        #println("Doubling done...")
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.J₀⁻[:] .= 0;
        temp = Array(exp.(-τ_λ./qp_μN'))
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ_λ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
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
        
    # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static)
    end
end

#
#Rotational Raman Scattering (default)
# Perform the Core RT routines (elemental, doubling, interaction)
function rt_kernel!(RS_type::Union{RRS, VS_0to1, VS_1to0}, 
                    pol_type, SFI, 
                    added_layer, composite_layer, 
                    computed_layer_properties, 
                    m, quad_points, 
                    I_static, 
                    architecture, 
                    qp_μN, iz) 
    
    @unpack τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface = computed_layer_properties
    @unpack F₀ = RS_type
    @unpack Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀ = RS_type
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show τ, ϖ, RS_type.fscattRayl
        
        @timeit "elemental_inelastic" elemental_inelastic!(RS_type, 
                                                pol_type, SFI, 
                                                τ_sum, dτ_λ, ϖ_λ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                F₀,
                                                m, ndoubl, scatter, 
                                                quad_points,  added_layer,  
                                                I_static, architecture)
        #@show m, added_layer.ieJ₀⁺[1:3, 1, 1, 1], added_layer.ieJ₀⁺[1:3, 1, end, end]                                             
        #println("Elemental inelastic done...")                                        
        @timeit "elemental" elemental!(pol_type, SFI, 
                                    τ_sum, dτ_λ, dτ, 
                                    ϖ_λ, ϖ, 
                                    Z⁺⁺, Z⁻⁺, 
                                    F₀,
                                    m, ndoubl, scatter, 
                                    quad_points,  added_layer,  
                                    I_static, architecture)
        #@show added_layer.J₀⁺[1:3, 1, 1], added_layer.J₀⁺[1:3, 1, end]                                             
        #println("Elemental  done...")
        @timeit "doubling_inelastic" doubling_inelastic!(
            RS_type, pol_type, 
            SFI, expk, ndoubl, 
            added_layer, I_static, architecture)
        #println("Doubling done...")
        #@timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.J₀⁻[:] .= 0;
        added_layer.ier⁻⁺[:] .= 0;
        added_layer.ier⁺⁻[:] .= 0;
        added_layer.ieJ₀⁻[:] .= 0;
        added_layer.iet⁻⁻[:] .= 0;
        added_layer.iet⁺⁺[:] .= 0;
        added_layer.ieJ₀⁺[:] .= 0;
        temp = Array(exp.(-τ_λ./qp_μN'))
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ_λ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
        end
    end

    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        composite_layer.T⁺⁺[:], composite_layer.T⁻⁻[:] = (added_layer.t⁺⁺, added_layer.t⁻⁻)
        composite_layer.R⁻⁺[:], composite_layer.R⁺⁻[:] = (added_layer.r⁻⁺, added_layer.r⁺⁻)
        composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.J₀⁺, added_layer.J₀⁻ )
        composite_layer.ieT⁺⁺[:], composite_layer.ieT⁻⁻[:] = (added_layer.iet⁺⁺, added_layer.iet⁻⁻)
        composite_layer.ieR⁻⁺[:], composite_layer.ieR⁺⁻[:] = (added_layer.ier⁻⁺, added_layer.ier⁺⁻)
        composite_layer.ieJ₀⁺[:], composite_layer.ieJ₀⁻[:] = (added_layer.ieJ₀⁺, added_layer.ieJ₀⁻ )
    
    # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static)
    end
end


### New update:
function rt_kernel!(RS_type::noRS{FT}, 
                    pol_type, SFI, 
                    added_layer, 
                    composite_layer, 
                    computed_layer_properties::CoreScatteringOpticalProperties, 
                    scattering_interface, 
                    τ_sum, 
                    m, quad_points, 
                    I_static, 
                    architecture, 
                    qp_μN, iz) where {FT}
    #@show array_type(architecture)
    @unpack qp_μ, μ₀ = quad_points
    @unpack F₀ = RS_type
    # Just unpack core optical properties from 
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
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
    expk = arr_type(exp.(-dτ /μ₀))
    
    #@show dτ, ndoubl
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show F₀
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        τ_sum, dτ, F₀,
                                        computed_layer_properties, 
                                        m, ndoubl, scatter, quad_points,  
                                        added_layer,  architecture)
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
                                        expk, ndoubl, 
                                        added_layer, 
                                        I_static, architecture)
        
        #=if m==0
            #m==0 ? 
            RayJ₀p = Array(added_layer.J₀⁺)
            RayJ₀m = Array(added_layer.J₀⁻)
            RayT   = Array(added_layer.t⁺⁺)
            RayR   = Array(added_layer.r⁻⁺)
            jldsave("/home/sanghavi/debugRay3.jld2"; RayJ₀p, RayJ₀m, RayT, RayR) 
        end=#                                
                                        
        
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
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
        end
    end

    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        composite_layer.T⁺⁺[:], composite_layer.T⁻⁻[:] = (added_layer.t⁺⁺, added_layer.t⁻⁻)
        composite_layer.R⁻⁺[:], composite_layer.R⁺⁻[:] = (added_layer.r⁻⁺, added_layer.r⁺⁻)
        composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.J₀⁺, added_layer.J₀⁻ )
    # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static)
        if iz==2
            M1 = Array(composite_layer.T⁺⁺);
            M2 = Array(composite_layer.R⁺⁻);
            M3 = Array(composite_layer.T⁻⁻);
            M4 = Array(composite_layer.R⁻⁺);
            M5 = Array(composite_layer.J₀⁻);
            M6 = Array(composite_layer.J₀⁺);
            #@show M1[1,1,1], M2[1,1,1], M3[1,1,1], M4[1,1,1], M5[1,1,1], M6[1,1,1]
        end
    end
    #=if m==0
        #m==0 ? 
        RayJ₀p = Array(composite_layer.J₀⁺)
        RayJ₀m = Array(composite_layer.J₀⁻)
        jldsave("/home/sanghavi/debugRay3.jld2"; RayJ₀p, RayJ₀m) 
    end=#
end

function rt_kernel!(
            RS_type::Union{RRS{FT}, VS_0to1{FT}, VS_1to0{FT}, 
                RRS_plus{FT}, VS_0to1_plus{FT}, VS_1to0_plus{FT}}, 
            pol_type, SFI, 
            added_layer, 
            composite_layer, 
            computed_layer_properties::CoreScatteringOpticalProperties, 
            scattering_interface, 
            τ_sum, m, quad_points, 
            I_static, architecture, qp_μN, iz)  where {FT}
    @unpack qp_μ, μ₀ = quad_points
    @unpack F₀ = RS_type
    # Just unpack core optical properties from 
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
    # SUNITI, check? Also, better to write function here
    dτ_max = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])
    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    scatter = true # edit later
    arr_type = array_type(architecture)
    # Compute dτ vector
    dτ = τ ./ 2^ndoubl
    expk = arr_type(exp.(-dτ /μ₀))
    #@show dτ, ndoubl
    @unpack Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀ = RS_type
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        
        @timeit "elemental_inelastic" elemental_inelastic!(RS_type, 
                                                pol_type, SFI, 
                                                τ_sum, dτ, ϖ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                F₀,
                                                m, ndoubl, scatter, 
                                                quad_points, added_layer,  
                                                I_static, architecture)
    
        
        #println("Elemental inelastic done...")                                        
        @timeit "elemental" elemental!(pol_type, SFI, τ_sum, dτ, F₀,
            computed_layer_properties, m, ndoubl, 
            scatter, quad_points, added_layer, architecture)

        #println("Elemental  done...")
       
        @timeit "doubling_inelastic" doubling_inelastic!(
            RS_type, pol_type, 
            SFI, expk, ndoubl, 
            added_layer, I_static, architecture)

            #=if m==0
                #m==0 ? 
                CabJ₀p = Array(added_layer.J₀⁺)
                CabJ₀m = Array(added_layer.J₀⁻)
                CabT   = Array(added_layer.t⁺⁺)
                CabR   = Array(added_layer.r⁻⁺)
                jldsave("/home/sanghavi/debugCab4.jld2"; CabJ₀p, CabJ₀m, CabT, CabR) 

                ieJ₀p = Array(added_layer.ieJ₀⁺)
                ieJ₀m = Array(added_layer.ieJ₀⁻)
                ieT   = Array(added_layer.iet⁺⁺)
                ieR   = Array(added_layer.ier⁻⁺)
                jldsave("/home/sanghavi/debugRRS4.jld2"; ieJ₀p, ieJ₀m, ieT, ieR) # : nothing
            end=#        
            
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.J₀⁻[:] .= 0;
        added_layer.ier⁻⁺[:] .= 0;
        added_layer.ier⁺⁻[:] .= 0;
        added_layer.ieJ₀⁻[:] .= 0;
        added_layer.iet⁻⁻[:] .= 0;
        added_layer.iet⁺⁺[:] .= 0;
        added_layer.ieJ₀⁺[:] .= 0;
        temp = Array(exp.(-τ_λ./qp_μN'))
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ_λ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
        end
    end

    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        composite_layer.T⁺⁺[:], composite_layer.T⁻⁻[:] = (added_layer.t⁺⁺, added_layer.t⁻⁻)
        composite_layer.R⁻⁺[:], composite_layer.R⁺⁻[:] = (added_layer.r⁻⁺, added_layer.r⁺⁻)
        composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.J₀⁺, added_layer.J₀⁻ )
        composite_layer.ieT⁺⁺[:], composite_layer.ieT⁻⁻[:] = (added_layer.iet⁺⁺, added_layer.iet⁻⁻)
        composite_layer.ieR⁻⁺[:], composite_layer.ieR⁺⁻[:] = (added_layer.ier⁻⁺, added_layer.ier⁺⁻)
        composite_layer.ieJ₀⁺[:], composite_layer.ieJ₀⁻[:] = (added_layer.ieJ₀⁺, added_layer.ieJ₀⁻ )
    
    # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, 
            SFI, composite_layer, added_layer, I_static)
    end

end
#=
function rt_kernel!(
            RS_type::Union{RRS_plus{FT}, VS_0to1_plus{FT}, VS_1to0_plus{FT}}, 
            pol_type, SFI, 
            added_layer, 
            composite_layer, 
            computed_layer_properties::CoreScatteringOpticalProperties, 
            scattering_interface, 
            τ_sum, m, quad_points, 
            I_static, architecture, qp_μN, iz)  where {FT}
    @unpack qp_μ, μ₀ = quad_points
    # Just unpack core optical properties from 
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
    # SUNITI, check? Also, better to write function here
    dτ_max = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])
    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    scatter = true # edit later
    arr_type = array_type(architecture)
    # Compute dτ vector
    dτ = τ ./ 2^ndoubl
    expk = arr_type(exp.(-dτ /μ₀))

    @unpack Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀ = RS_type
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show τ, ϖ, RS_type.fscattRayl
        
        @timeit "elemental_inelastic" elemental_inelastic!(RS_type, 
                                                pol_type, SFI, 
                                                τ_sum, dτ, ϖ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                m, ndoubl, scatter, 
                                                quad_points, added_layer,  
                                                I_static, architecture)
        #println("Elemental inelastic done...")                                        
        @timeit "elemental" elemental!(pol_type, SFI, τ_sum, dτ, 
            computed_layer_properties, m, ndoubl, 
            scatter, quad_points, added_layer, architecture)
        #println("Elemental  done...")
        @timeit "doubling_inelastic" doubling_inelastic!(RS_type, pol_type, 
                    SFI, expk, ndoubl, added_layer, I_static, architecture)
        #println("Doubling done...")
        #@timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.J₀⁻[:] .= 0;
        added_layer.ier⁻⁺[:] .= 0;
        added_layer.ier⁺⁻[:] .= 0;
        added_layer.ieJ₀⁻[:] .= 0;
        added_layer.iet⁻⁻[:] .= 0;
        added_layer.iet⁺⁺[:] .= 0;
        added_layer.ieJ₀⁺[:] .= 0;
        temp = Array(exp.(-τ./qp_μN'))
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
        end
    end

    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        composite_layer.T⁺⁺[:], composite_layer.T⁻⁻[:] = (added_layer.t⁺⁺, added_layer.t⁻⁻)
        composite_layer.R⁻⁺[:], composite_layer.R⁺⁻[:] = (added_layer.r⁻⁺, added_layer.r⁺⁻)
        composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.J₀⁺, added_layer.J₀⁻ )
        composite_layer.ieT⁺⁺[:], composite_layer.ieT⁻⁻[:] = (added_layer.iet⁺⁺, added_layer.iet⁻⁻)
        composite_layer.ieR⁻⁺[:], composite_layer.ieR⁺⁻[:] = (added_layer.ier⁻⁺, added_layer.ier⁺⁻)
        composite_layer.ieJ₀⁺[:], composite_layer.ieJ₀⁻[:] = (added_layer.ieJ₀⁺, added_layer.ieJ₀⁻ )
    
    # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, 
            SFI, composite_layer, added_layer, I_static)
    end
end
=#