#=
 
This file implements rt_kernel!, which performs the core RT routines (elemental, doubling, interaction)
 
=#
#No Raman (default)
# Perform the Core RT routines (elemental, doubling, interaction)

### 
function rt_kernel!(RS_type::noRS{FT}, 
                    pol_type, SFI, 
                    added_layer, lin_added_layer,
                    composite_layer, lin_composite_layer,
                    computed_layer_properties::M, 
                    lin_computed_layer_properties, 
                    scattering_interface, 
                    τ_sum, 
                    lin_τ_sum,
                    m, quad_points, 
                    I_static, 
                    architecture, 
                    qp_μN, iz) where {FT,M}
    #@show array_type(architecture)
    
    @unpack qp_μ, μ₀ = quad_points
    # Just unpack core optical properties from 
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
    @unpack lin_τ, lin_ϖ, lin_Z⁺⁺, lin_Z⁻⁺ = lin_computed_layer_properties
    
    nparams = size(lin_τ,1)   
    # @show ndoubl
    scatter = true # edit later
    
    dτ, lin_dτ, ndoubl, expk, dexpk_fctr = init_layer(
                        computed_layer_properties, 
                        lin_computed_layer_properties, 
                        quad_points, pol_type, architecture)

    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show typeof(computed_layer_properties)
        @timeit "elemental" elemental!(pol_type, SFI, 
                                τ_sum, #lin_τ_sum, 
                                dτ, #lin_dτ,
                                computed_layer_properties, 
                                #lin_computed_layer_properties,
                                m, ndoubl, scatter, quad_points,  
                                added_layer,  
                                lin_added_layer,
                                architecture)
        #println("Elemental done...")
        @timeit "doubling"   doubling!(pol_type, SFI, 
                                expk, #lin_expk, 
                                dexpk_fctr,
                                ndoubl, 
                                added_layer, lin_added_layer,
                                I_static, architecture)
        #@show added_layer.r⁻⁺[1:2,1,1], added_layer.r⁺⁻[1:2,1,1],added_layer.t⁺⁺[1:2,1,1], added_layer.t⁻⁻[1:2,1,1] 
        #@show dτ, ndoubl, expk
    #println("Doubling done...")
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.j₀⁻[:] .= 0;
        lin_added_layer.dr⁻⁺[:] .= 0;
        lin_added_layer.dr⁺⁻[:] .= 0;
        lin_added_layer.dj₀⁻[:] .= 0;
        temp = Array(exp.(-τ_λ./qp_μN'))
        dtemp = Array(-exp.(-τ_λ./qp_μN')./qp_μN')
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ_λ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
            lin_added_layer.dt⁺⁺[1,:,:,iλ] = Diagonal(dtemp[iλ,:]);
            lin_added_layer.dt⁻⁻[1,:,:,iλ] = Diagonal(dtemp[iλ,:]);
        end
    end

    # @assert !any(isnan.(added_layer.t⁺⁺))

    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        
        composite_layer.T⁺⁺[:], composite_layer.T⁻⁻[:] = (added_layer.t⁺⁺, added_layer.t⁻⁻)
        composite_layer.R⁻⁺[:], composite_layer.R⁺⁻[:] = (added_layer.r⁻⁺, added_layer.r⁺⁻)
        composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.j₀⁺, added_layer.j₀⁻)
        
        for ctr=1:nparams
            lin_composite_layer.dT⁺⁺[ctr,:,:,:] = lin_τ[ctr,:].*lin_added_layer.dt⁺⁺[1,:,:,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dt⁺⁺[2,:,:,:] + 
                                                lin_Z⁺⁺[ctr,:,:,:].*lin_added_layer.dt⁺⁺[3,:,:,:] +
    
            lin_composite_layer.dT⁻⁻[ctr,:,:,:] = lin_τ[ctr,:].*lin_added_layer.dt⁻⁻[1,:,:,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dt⁻⁻[2,:,:,:] + 
                                                lin_Z⁺⁺[ctr,:,:,:].*lin_added_layer.dt⁻⁻[3,:,:,:] +

            lin_composite_layer.dR⁻⁺[ctr,:,:,:] = lin_τ[ctr,:].*lin_added_layer.dr⁻⁺[1,:,:,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dr⁻⁺[2,:,:,:] +
                                                lin_Z⁻⁺[ctr,:,:,:].*lin_added_layer.dr⁻⁺[3,:,:,:]
            
            lin_composite_layer.dR⁺⁻[ctr,:,:,:] = lin_τ[ctr,:].*lin_added_layer.dr⁺⁻[1,:,:,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dr⁺⁻[2,:,:,:] +
                                                lin_Z⁻⁺[ctr,:,:,:].*lin_added_layer.dr⁺⁻[3,:,:,:]

            lin_composite_layer.dJ₀⁺[ctr,:,1,:] = lin_τ[ctr,:].*lin_added_layer.dj₀⁺[1,:,1,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dj₀⁺[2,:,1,:] +
                                                lin_Z⁺⁺[ctr,:,:,:].*lin_added_layer.dj₀⁺[3,:,1,:] +
                                                lin_τ_sum[ctr,:,:,:].*lin_added_layer.dj₀⁺[4,:,1,:]
    
            lin_composite_layer.dJ₀⁻[ctr,:,1,:] = lin_τ[ctr,:].*lin_added_layer.dj₀⁻[1,:,1,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dj₀⁻[2,:,1,:] +
                                                lin_Z⁻⁺[ctr,:,:,:].*lin_added_layer.dj₀⁻[3,:,1,:] +
                                                lin_τ_sum[ctr,:,:,:].*lin_added_layer.dj₀⁻[4,:,1,:]
        end
        # If this is not the TOA, perform the interaction step
    else        
        for ctr=1:nparams
            lin_added_layer.dxt⁺⁺[ctr,:,:,:] = lin_τ[ctr,:].*lin_added_layer.dt⁺⁺[1,:,:,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dt⁺⁺[2,:,:,:] + 
                                                lin_Z⁺⁺[ctr,:,:,:].*lin_added_layer.dt⁺⁺[3,:,:,:] +
    
            lin_added_layer.dxt⁻⁻[ctr,:,:,:] = lin_τ[ctr,:].*lin_added_layer.dt⁻⁻[1,:,:,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dt⁻⁻[2,:,:,:] + 
                                                lin_Z⁺⁺[ctr,:,:,:].*lin_added_layer.dt⁻⁻[3,:,:,:] +

            lin_added_layer.dxr⁻⁺[ctr,:,:,:] = lin_τ[ctr,:].*lin_added_layer.dr⁻⁺[1,:,:,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dr⁻⁺[2,:,:,:] +
                                                lin_Z⁻⁺[ctr,:,:,:].*lin_added_layer.dr⁻⁺[3,:,:,:]
            
            lin_added_layer.dxr⁺⁻[ctr,:,:,:] = lin_τ[ctr,:].*lin_added_layer.dr⁺⁻[1,:,:,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dr⁺⁻[2,:,:,:] +
                                                lin_Z⁻⁺[ctr,:,:,:].*lin_added_layer.dr⁺⁻[3,:,:,:]

            lin_added_layer.dxj₀⁺[ctr,:,1,:] = lin_τ[ctr,:].*lin_added_layer.dj₀⁺[1,:,1,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dj₀⁺[2,:,1,:] +
                                                lin_Z⁺⁺[ctr,:,:,:].*lin_added_layer.dj₀⁺[3,:,1,:] +
                                                lin_τ_sum[ctr,:,:,:].*lin_added_layer.dj₀⁺[4,:,1,:]
    
            lin_added_layer.dxj₀⁻[ctr,:,1,:] = lin_τ[ctr,:].*lin_added_layer.dj₀⁻[1,:,1,:] +
                                                lin_ϖ[ctr]*lin_added_layer.dj₀⁻[2,:,1,:] +
                                                lin_Z⁻⁺[ctr,:,:,:].*lin_added_layer.dj₀⁻[3,:,1,:] +
                                                lin_τ_sum[ctr,:,:,:].*lin_added_layer.dj₀⁻[4,:,1,:]
        end
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, lin_composite_layer, added_layer, lin_added_layer, I_static)
    end
end

#=
function get_dtau_ndoubl(computed_layer_properties::CoreScatteringOpticalProperties,
        lin_computed_layer_properties::linCoreScatteringOpticalProperties, 
        quad_points::QuadPoints{FT}) where {FT}
    @unpack qp_μ  = quad_points
    @unpack τ, ϖ  = computed_layer_properties
    @unpack lin_τ    = lin_computed_layer_properties
    dτ_max = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])
    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    # Compute dτ vector
    dτ = τ ./ 2^ndoubl
    lin_dτ = lin_τ / 2^ndoubl

    return dτ, lin_dτ, ndoubl
end

function get_dtau_ndoubl(computed_layer_properties::CoreDirectionalScatteringOpticalProperties, 
        lin_computed_layer_properties::linCoreDirectionalScatteringOpticalProperties, 
        quad_points::QuadPoints{FT}) where {FT}
    @unpack qp_μ,iμ₀  = quad_points
    @unpack τ, ϖ, G  = computed_layer_properties
    @unpack lin_τ    = lin_computed_layer_properties
    gfct = Array(G)[iμ₀]
    dτ_max = minimum([maximum(gfct * τ .* ϖ), FT(0.001) * minimum(qp_μ)])
    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    # Compute dτ vector
    dτ = τ ./ 2^ndoubl
    lin_dτ = lin_τ / 2^ndoubl

    return dτ, lin_dτ, ndoubl
end

function init_layer(computed_layer_properties::CoreDirectionalScatteringOpticalProperties, 
                    lin_computed_layer_properties::linCoreDirectionalScatteringOpticalProperties,
                    quad_points, pol_type, architecture)
    arr_type = array_type(architecture) 
    @unpack μ₀, iμ₀ = quad_points
    @unpack G = computed_layer_properties
    dτ, lin_dτ, ndoubl = get_dtau_ndoubl(
                            computed_layer_properties, 
                            lin_computed_layer_properties, 
                            quad_points)
    gfct = Array(G)[iμ₀]
    expk = exp.(-dτ*gfct/μ₀)
    #lin_expk = -expk.*lin_dτ*gfct/μ₀
    dexpk_fctr = -gfct/μ₀

    return dτ, lin_dτ, ndoubl, arr_type(expk), FT(dexpk_fctr) #arr_type(lin_expk)
end

function init_layer(computed_layer_properties::CoreScatteringOpticalProperties, 
                    lin_computed_layer_properties::linCoreScatteringOpticalProperties, 
                    quad_points, pol_type, architecture)
    arr_type = array_type(architecture)
    @unpack μ₀ = quad_points
    dτ, lin_dτ, ndoubl = get_dtau_ndoubl(
                            computed_layer_properties, 
                            lin_computed_layer_properties, 
                            quad_points)
    expk = exp.(-dτ/μ₀)
    #lin_expk = -expk.*lin_dτ/μ₀
    dexpk_fctr = -1/μ₀

    return dτ, lin_dτ, ndoubl, arr_type(expk), FT(dexpk_fctr)#arr_type(lin_expk)
end

=#
function rt_kernel!(RS_type::Union{RRS{FT}, VS_0to1{FT}, VS_1to0{FT}}, 
        pol_type, SFI, 
        added_layer, lin_added_layer,
        composite_layer, lin_composite_layer,
        computed_layer_properties::CoreScatteringOpticalProperties, 
        lin_computed_layer_properties::linCoreScatteringOpticalProperties, 
        scattering_interface, 
        τ_sum,
        lin_τ_sum,
        m, quad_points, 
        I_static, architecture, qp_μN, iz)  where {FT}
    @unpack qp_μ, μ₀ = quad_points
    # Just unpack core optical properties from 
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
    @unpack lin_τ, lin_ϖ, lin_Z⁺⁺, lin_Z⁻⁺ = lin_computed_layer_properties
    
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
                                                τ_sum, lin_τ_sum, 
                                                dτ, ϖ, 
                                                lin_dτ, lin_ϖ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                m, ndoubl, scatter, 
                                                quad_points,  
                                                added_layer, lin_added_layer,
                                                I_static, architecture)
        #println("Elemental inelastic done...")                                        
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        τ_sum, lin_τ_sum,   
                                        dτ, lin_dτ,
                                        computed_layer_properties, 
                                        lin_computed_layer_properties,
                                        m, ndoubl, scatter, quad_points,  
                                        added_layer, lin_added_layer,  
                                        architecture)
        #println("Elemental  done...")
        @timeit "doubling_inelastic" doubling_inelastic!(RS_type, 
                                        pol_type, SFI, 
                                        expk, lin_expk,
                                        ndoubl, 
                                        added_layer, lin_added_layer,
                                        I_static, architecture)
        #println("Doubling done...")
        #@timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.j₀⁻[:] .= 0;
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
        composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.j₀⁺, added_layer.j₀⁻ )
        composite_layer.ieT⁺⁺[:], composite_layer.ieT⁻⁻[:] = (added_layer.iet⁺⁺, added_layer.iet⁻⁻)
        composite_layer.ieR⁻⁺[:], composite_layer.ieR⁺⁻[:] = (added_layer.ier⁻⁺, added_layer.ier⁺⁻)
        composite_layer.ieJ₀⁺[:], composite_layer.ieJ₀⁻[:] = (added_layer.ieJ₀⁺, added_layer.ieJ₀⁻ )
    
    # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, 
                            composite_layer, lin_composite_layer, 
                            added_layer, lin_added_layer,
                            I_static)
    end
end



function rt_kernel!(
            RS_type::Union{RRS_plus{FT}, VS_0to1_plus{FT}, VS_1to0_plus{FT}}, 
            pol_type, SFI, 
            added_layer, 
            composite_layer, 
            computed_layer_properties::CoreScatteringOpticalProperties, 
            lin_computed_layer_properties::linCoreScatteringOpticalProperties, 
            scattering_interface, 
            τ_sum,
            lin_τ_sum,
            m, quad_points, 
            I_static, architecture, qp_μN, iz)  where {FT}
    @unpack qp_μ, μ₀ = quad_points
    # Just unpack core optical properties from 
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
    @unpack lin_τ, lin_ϖ, lin_Z⁺⁺, lin_Z⁻⁺ = lin_computed_layer_properties
    
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
                                                τ_sum, lin_τ_sum, 
                                                dτ, ϖ, 
                                                lin_dτ, lin_ϖ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                m, ndoubl, scatter, 
                                                quad_points,  
                                                added_layer, lin_added_layer,  
                                                I_static, architecture)
        #println("Elemental inelastic done...")                                        
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        τ_sum, lin_τ_sum, 
                                        dτ, lin_dτ, 
                                        computed_layer_properties, 
                                        lin_computed_layer_properties,
                                        m, ndoubl, scatter, 
                                        quad_points,  
                                        added_layer,  lin_added_layer,
                                        architecture)
        #println("Elemental  done...")
        @timeit "doubling_inelastic" doubling_inelastic!(RS_type, pol_type, 
                                        SFI, 
                                        expk, lin_expk, 
                                        ndoubl, 
                                        added_layer, lin_added_layer, 
                                        I_static, architecture)
        #println("Doubling done...")
        #@timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.j₀⁻[:] .= 0;
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
        composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.j₀⁺, added_layer.j₀⁻ )
        composite_layer.ieT⁺⁺[:], composite_layer.ieT⁻⁻[:] = (added_layer.iet⁺⁺, added_layer.iet⁻⁻)
        composite_layer.ieR⁻⁺[:], composite_layer.ieR⁺⁻[:] = (added_layer.ier⁻⁺, added_layer.ier⁺⁻)
        composite_layer.ieJ₀⁺[:], composite_layer.ieJ₀⁻[:] = (added_layer.ieJ₀⁺, added_layer.ieJ₀⁻ )
    
    # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, 
                composite_layer, lin_composite_layer,
                added_layer, lin_added_layer, I_static)
    end
end


