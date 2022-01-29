#=
 
This file implements rt_kernel!, which performs the core RT routines (elemental, doubling, interaction)
 
=#
#No Raman (default)
# Perform the Core RT routines (elemental, doubling, interaction)
function rt_kernel!(RS_type::noRS, pol_type, SFI, added_layer, composite_layer, computed_layer_properties, m, quad_points, I_static, architecture, qp_μN, iz) 

    @unpack τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface = computed_layer_properties

    # If there is scattering, perform the elemental and doubling steps
    if scatter
        
        @timeit "elemental" elemental!(pol_type, SFI, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, quad_points,  added_layer,  I_static, architecture)
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
function rt_kernel!(RS_type::Union{RRS, VS_0to1, VS_1to0}, pol_type, SFI, added_layer, composite_layer, computed_layer_properties, m, quad_points, I_static, architecture, qp_μN, iz) 
    
    @unpack τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface = computed_layer_properties
    @unpack Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀ = RS_type
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show τ, ϖ, RS_type.fscattRayl
        @timeit "elemental_inelastic" elemental_inelastic!(RS_type, 
                                                pol_type, SFI, 
                                                τ_sum, dτ_λ, ϖ_λ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                m, ndoubl, scatter, 
                                                quad_points,  added_layer,  
                                                I_static, architecture)
        #println("Elemental inelastic done...")                                        
        @timeit "elemental" elemental!(pol_type, SFI, 
                                    τ_sum, dτ_λ, dτ, 
                                    ϖ_λ, ϖ, 
                                    Z⁺⁺, Z⁻⁺, 
                                    m, ndoubl, scatter, 
                                    quad_points,  added_layer,  
                                    I_static, architecture)
        #println("Elemental  done...")
        @timeit "doubling_inelastic" doubling_inelastic!(RS_type, pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
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