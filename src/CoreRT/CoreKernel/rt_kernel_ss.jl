#=
 
This file implements rt_kernel_ss!, which performs the core RT routines (elemental, doubling, interaction) to compute ONLY the single scattering compoent of the total signal
 
=#
#No Raman (default)
# Perform the Core RT routines (elemental, doubling, interaction)
function rt_kernel_ss!(RS_type::noRS, pol_type, SFI, added_layer, composite_layer, computed_layer_properties, m, quad_points, I_static, architecture, qp_μN, iz) 

    (; τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface) = computed_layer_properties
    (; F₀) = RS_type

    ndoubl_ss = Int(0)
    #@show τ, ϖ, dτ_max, ndoubl
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        τ_sum, τ_λ, τ, 
                                        ϖ_λ, ϖ, 
                                        Z⁺⁺, Z⁻⁺, 
                                        F₀,
                                        m, ndoubl_ss, 
                                        scatter, 
                                        quad_points,  
                                        added_layer,  
                                        I_static, 
                                        architecture)
        #println("Elemental done...")
        #@timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
        #println("Doubling done...")
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.j₀⁻[:] .= 0;
        temp = collect(exp.(-τ_λ./qp_μN'))
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ_λ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
        end
    end
    #M1 = Array(added_layer.t⁺⁺)
    #M2 = Array(added_layer.r⁺⁻)
    #M3 = Array(added_layer.j₀⁻)
    #M4 = Array(added_layer.j₀⁺)
    #@show M1[1,1,1], M2[1,1,1], M3[1,1,1], M4[1,1,1]
    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        composite_layer.T⁺⁺[:], composite_layer.T⁻⁻[:] = (added_layer.t⁺⁺, added_layer.t⁻⁻)
        composite_layer.R⁻⁺[:], composite_layer.R⁺⁻[:] = (added_layer.r⁻⁺, added_layer.r⁺⁻)
        composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.j₀⁺, added_layer.j₀⁻ )
        
    # If this is not the TOA, perform the interaction step
    else
        interaction_ss!(SFI,
            composite_layer, 
            added_layer, 
            τ_sum,
            τ_λ,
            quad_points,
            architecture)
        #@timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static)
    end
end

#
#Rotational Raman Scattering (default)
# Perform the Core RT routines (elemental, doubling, interaction)
function rt_kernel_ss!(RS_type::Union{RRS, VS_0to1, VS_1to0}, pol_type, SFI, added_layer, composite_layer, computed_layer_properties, m, quad_points, I_static, architecture, qp_μN, iz) 
    
    (; τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface) = computed_layer_properties
    (; F₀) = RS_type
    (; Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀) = RS_type
    # If there is scattering, perform the elemental and doubling steps
    ndoubl_ss = Int(0)
    if scatter
        #@show τ, ϖ, RS_type.fscattRayl
        @timeit "elemental_inelastic" elemental_inelastic!(RS_type, 
                                                pol_type, SFI, 
                                                τ_sum, τ_λ, ϖ_λ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                F₀,
                                                m, ndoubl_ss, scatter, 
                                                quad_points,  added_layer,  
                                                I_static, architecture)
        #println("Elemental inelastic done...")                                        
        @timeit "elemental" elemental!(pol_type, SFI, 
                                    τ_sum, τ_λ, τ, 
                                    ϖ_λ, ϖ, 
                                    Z⁺⁺, Z⁻⁺, 
                                    F₀,
                                    m, ndoubl_ss, scatter, 
                                    quad_points,  added_layer,  
                                    I_static, architecture)
        #println("Elemental  done...")
        #@timeit "doubling_inelastic" doubling_inelastic!(RS_type, pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
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
        temp = collect(exp.(-τ_λ./qp_μN'))
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
        interaction_ss!(SFI,
            composite_layer, 
            added_layer, 
            τ_sum,
            τ_λ,
            quad_points,
            architecture)
        
        interaction_inelastic_ss!(RS_type,
            SFI,
            composite_layer, 
            added_layer, 
            τ_sum,
            τ_λ,
            quad_points,
            architecture)
        #@timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static)
    end
end

### New update:
function rt_kernel_ss!(RS_type::noRS{FT},
                    pol_type, SFI,
                    added_layer,
                    composite_layer,
                    computed_layer_properties::CoreScatteringOpticalProperties,
                    scattering_interface,
                    τ_sum,
                    m, quad_points,
                    I_static,
                    architecture,
                    qp_μN, iz;
                    dτ_max_threshold::Union{Nothing,Real} = nothing,
                    dτ_min_floor::Union{Nothing,Real} = nothing) where {FT}
    (; μ₀) = quad_points
    (; F₀) = RS_type
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = computed_layer_properties
    # Centralised dτ/ndoubl — same helper as rt_kernel! to keep the
    # zero-weight-stream filter and absolute floor consistent across the
    # full RT, single-scatter, and multi-sensor paths.
    dτ, ndoubl = get_dtau_ndoubl(computed_layer_properties, quad_points;
                                 dτ_max_threshold = dτ_max_threshold,
                                 dτ_min_floor = dτ_min_floor)
    scatter = true # edit later
    arr_type = array_type(architecture)
    expk = arr_type(exp.(-dτ /μ₀))

    ndoubl_ss = Int(0)
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        τ_sum, τ, F₀,
                                        computed_layer_properties, 
                                        m, ndoubl_ss, scatter, quad_points,  
                                        added_layer,  architecture)
        #println("Elemental done...")
        #@timeit "doubling"   doubling!(pol_type, SFI, 
        #                                expk, ndoubl, 
        #                                added_layer, 
        #                                I_static, architecture)
        #println("Doubling done...")
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.j₀⁻[:] .= 0;
        temp = collect(exp.(-τ./qp_μN'))
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
    # If this is not the TOA, perform the interaction step
    else
        #@timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static)
        interaction_ss!(SFI,
            composite_layer, 
            added_layer, 
            τ_sum,
            τ,
            quad_points,
            architecture)
        
        if iz==2
            M1 = collect(composite_layer.T⁺⁺);
            M2 = collect(composite_layer.R⁺⁻);
            M3 = collect(composite_layer.T⁻⁻);
            M4 = collect(composite_layer.R⁻⁺);
            M5 = collect(composite_layer.J₀⁻);
            M6 = collect(composite_layer.J₀⁺);
            #@show M1[1,1,1], M2[1,1,1], M3[1,1,1], M4[1,1,1], M5[1,1,1], M6[1,1,1]
        end
    end
end

function rt_kernel_ss!(
            RS_type::Union{RRS{FT}, VS_0to1{FT}, VS_1to0{FT},
                VS_0to1_plus{FT}, VS_1to0_plus{FT}},
            pol_type, SFI,
            added_layer,
            composite_layer,
            computed_layer_properties::CoreScatteringOpticalProperties,
            scattering_interface,
            τ_sum, m, quad_points,
            I_static, architecture, qp_μN, iz;
            dτ_max_threshold::Union{Nothing,Real} = nothing,
            dτ_min_floor::Union{Nothing,Real} = nothing)  where {FT}
    (; μ₀) = quad_points
    (; F₀) = RS_type
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = computed_layer_properties
    # Centralised dτ/ndoubl helper (see noRS rt_kernel_ss! variant).
    dτ, ndoubl = get_dtau_ndoubl(computed_layer_properties, quad_points;
                                 dτ_max_threshold = dτ_max_threshold,
                                 dτ_min_floor = dτ_min_floor)
    scatter = true # edit later
    arr_type = array_type(architecture)
    expk = arr_type(exp.(-dτ /μ₀))

    (; Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀) = RS_type
    ndoubl_ss = Int(0)
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show τ, ϖ, RS_type.fscattRayl
        
        @timeit "elemental_inelastic" elemental_inelastic!(RS_type, 
                                                pol_type, SFI, 
                                                τ_sum, τ, ϖ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                F₀,
                                                m, ndoubl_ss, scatter, 
                                                quad_points, added_layer,  
                                                I_static, architecture)
        #println("Elemental inelastic done...")                                        
        @timeit "elemental" elemental!(pol_type, SFI, τ_sum, τ, F₀,
            computed_layer_properties, m, ndoubl_ss, 
            scatter, quad_points, added_layer, architecture)
        #println("Elemental  done...")
        #@timeit "doubling_inelastic" doubling_inelastic!(RS_type, pol_type, 
        #            SFI, expk, ndoubl, added_layer, I_static, architecture)
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
        temp = collect(exp.(-τ./qp_μN'))
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
        #@timeit "interaction" interaction!(RS_type, scattering_interface, 
        #    SFI, composite_layer, added_layer, I_static)
        interaction_ss!(SFI,
            composite_layer, 
            added_layer, 
            τ_sum,
            τ,
            quad_points,
            architecture)
        
        interaction_inelastic_ss!(RS_type,
            SFI,
            composite_layer, 
            added_layer, 
            τ_sum,
            τ,
            quad_points,
            architecture)
    end
end
#=
function rt_kernel!(
            RS_type::Union{VS_0to1_plus{FT}, VS_1to0_plus{FT}},
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
        added_layer.j₀⁻[:] .= 0;
        added_layer.ier⁻⁺[:] .= 0;
        added_layer.ier⁺⁻[:] .= 0;
        added_layer.ieJ₀⁻[:] .= 0;
        added_layer.iet⁻⁻[:] .= 0;
        added_layer.iet⁺⁺[:] .= 0;
        added_layer.ieJ₀⁺[:] .= 0;
        temp = collect(exp.(-τ./qp_μN'))
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
        @timeit "interaction" interaction!(RS_type, scattering_interface, 
            SFI, composite_layer, added_layer, I_static)
    end
end
=#
