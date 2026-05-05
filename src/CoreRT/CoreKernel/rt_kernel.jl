#=
 
This file implements rt_kernel!, which performs the core RT routines (elemental, doubling, interaction)
 
=#

"""
    rt_kernel!

Core radiative-transfer kernel that processes a single atmospheric layer.

For each layer (indexed by `iz` from top-of-atmosphere downward):

1. **Elemental** – builds the thin single-scattering layer matrices (r, t, j).
2. **Doubling** – doubles the elemental layer `ndoubl` times to obtain the
   full homogeneous layer.
3. **Interaction** – combines the newly doubled (added) layer with the
   accumulated composite layer from above.

At the first layer (`iz == 1`) the added-layer matrices are simply copied
into the composite layer.  For non-scattering layers, Beer-law transmission
is assigned directly.

Multiple methods are dispatched on the Raman-scattering type (`noRS`, `RRS`,
`VS_0to1`, etc.) and on the optical-property container type
(`CoreScatteringOpticalProperties` vs. pre-unpacked fields).
"""

"""
    _set_transmission_noscat!(t⁺⁺, t⁻⁻, τ_vals, qp_μN)

Set transmission matrices for non-scattering layers using Beer's law.
Constructs batch diagonal matrices: `t[j,j,iλ] = exp(-τ[iλ]/μ[j])`.
Works on both CPU and GPU arrays without device transfers.
"""
@inline function _set_transmission_noscat!(t⁺⁺, t⁻⁻, τ_vals, qp_μN)
    fill!(t⁺⁺, 0)
    fill!(t⁻⁻, 0)
    temp = exp.(-τ_vals ./ qp_μN')
    nμ = size(temp, 2)
    for j in 1:nμ
        @views t⁺⁺[j,j,:] .= temp[:, j]
        @views t⁻⁻[j,j,:] .= temp[:, j]
    end
end
#No Raman (default)
# Perform the Core RT routines (elemental, doubling, interaction)
function rt_kernel!(RS_type::noRS, 
            pol_type, SFI, 
            added_layer, 
            composite_layer, 
            computed_layer_properties, 
            m, 
            quad_points, 
            I_static, 
            architecture, 
            qp_μN, iz) 

    (; τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface) = computed_layer_properties
    (; F₀) = RS_type
    # If there is scattering, perform the elemental and doubling steps
    if scatter

        @timeit "elemental" elemental!(pol_type, SFI, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, F₀, m, ndoubl, scatter, quad_points,  added_layer,  I_static, architecture)
        #println("Elemental done...")
        @timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
        #println("Doubling done...")
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        zero_added_noscat!(added_layer, τ_λ, qp_μN)
    end
    #M1 = Array(added_layer.t⁺⁺)
    #M2 = Array(added_layer.r⁺⁻)
    #M3 = Array(added_layer.j₀⁻)
    #M4 = Array(added_layer.j₀⁺)
    #@show M1[1,1,1], M2[1,1,1], M3[1,1,1], M4[1,1,1]
    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        copy_added_to_composite!(composite_layer, added_layer)
        
    # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static)
    end
end

"""
    rt_kernel_canopy!(RS_type::noRS, pol_type, SFI, added_layer, composite_layer, ...)

RT kernel for canopy layers: uses `elemental_canopy!` instead of `elemental!` for
directional scattering (G factor). Otherwise identical to `rt_kernel!`.
"""
function rt_kernel_canopy!(RS_type::noRS, pol_type, SFI, added_layer, composite_layer, computed_layer_properties, m, quad_points, I_static, architecture, qp_μN, iz) 

    (; τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface) = computed_layer_properties
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        
        @timeit "elemental_canopy" elemental_canopy!(pol_type, SFI, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, quad_points,  added_layer,  I_static, architecture)
        #println("Elemental done...")
        @timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
        #println("Doubling done...")
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        zero_added_noscat!(added_layer, τ_λ, qp_μN)
    end
    #M1 = Array(added_layer.t⁺⁺)
    #M2 = Array(added_layer.r⁺⁻)
    #M3 = Array(added_layer.j₀⁻)
    #M4 = Array(added_layer.j₀⁺)
    #@show M1[1,1,1], M2[1,1,1], M3[1,1,1], M4[1,1,1]
    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        copy_added_to_composite!(composite_layer, added_layer)
        
    # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static)
    end
end

#
#Rotational Raman Scattering (default)
# Perform the Core RT routines (elemental, doubling, interaction)
function rt_kernel!(RS_type::Union{RRS, VS_0to1, VS_1to0}, pol_type, SFI, added_layer, composite_layer, computed_layer_properties, m, quad_points, I_static, architecture, qp_μN, iz) 
    
    (; τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface) = computed_layer_properties
    (; Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀) = RS_type
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show τ, ϖ, RS_type.fscattRayl
        @timeit "elemental_inelastic" elemental_inelastic!(RS_type, 
                                                pol_type, SFI, 
                                                τ_sum, dτ_λ, ϖ_λ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                RS_type.F₀,
                                                m, ndoubl, scatter, 
                                                quad_points,  added_layer,  
                                                I_static, architecture)
        #println("Elemental inelastic done...")                                        
        @timeit "elemental" elemental!(pol_type, SFI,
                                    τ_sum, dτ_λ, dτ,
                                    ϖ_λ, ϖ,
                                    Z⁺⁺, Z⁻⁺,
                                    RS_type.F₀,
                                    m, ndoubl, scatter,
                                    quad_points,  added_layer,
                                    I_static, architecture)
        #println("Elemental  done...")
        @timeit "doubling_inelastic" doubling_inelastic!(RS_type, pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
        #println("Doubling done...")
        #@timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        zero_added_noscat_ie!(added_layer, τ_λ, qp_μN)
    end

    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        copy_added_to_composite_ie!(composite_layer, added_layer)
    
    # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static)
    end
end

### 
function rt_kernel!(RS_type::noRS{FT},
                    pol_type, SFI,
                    added_layer,
                    composite_layer,
                    computed_layer_properties::M,
                    scattering_interface,
                    τ_sum,
                    m, quad_points,
                    I_static,
                    architecture,
                    qp_μN, iz;
                    workspace=nothing) where {FT,M}
    #@show array_type(architecture)
    
    (; qp_μ, μ₀) = quad_points
    (; F₀) = RS_type
    # Just unpack core optical properties from
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = computed_layer_properties

    scatter = maximum(τ .* ϖ) > 2eps(FT)

    # If there is scattering, perform the elemental and doubling steps
    if scatter
        dτ, ndoubl, expk = init_layer(computed_layer_properties, quad_points, pol_type, architecture)
        #@show typeof(computed_layer_properties)
        @timeit "elemental" elemental!(pol_type, SFI,
                                τ_sum, dτ, F₀,
                                computed_layer_properties,
                                m, ndoubl, scatter, quad_points,
                                added_layer,  architecture)
        #println("Elemental done...")
        @timeit "doubling"   doubling!(pol_type, SFI, 
                                expk, ndoubl, 
                                added_layer,
                                I_static, architecture)
        #@show added_layer.r⁻⁺[1:2,1,1], added_layer.r⁺⁻[1:2,1,1],added_layer.t⁺⁺[1:2,1,1], added_layer.t⁻⁻[1:2,1,1] 
        #@show dτ, ndoubl, expk
    #println("Doubling done...")
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        zero_added_noscat!(added_layer, τ, qp_μN)
    end

    # @assert !any(isnan.(added_layer.t⁺⁺))

    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        copy_added_to_composite!(composite_layer, added_layer)
        # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static)
    end
end


"""
    get_dtau_ndoubl(computed_layer_properties, quad_points)

Compute elemental optical depth `dτ` and doubling count `ndoubl` for the adding-doubling method.

The elemental layer optical depth is ``d\\tau = \\tau / 2^{n_d}``, where `ndoubl` is chosen so that
``d\\tau \\cdot \\varpi`` is sufficiently small for accurate single-scattering (typically < 0.001).
Uses `doubling_number` to determine `ndoubl` from the maximum allowed `dτ_max`.

# Returns
- `dτ`: Elemental optical depth vector `[nSpec]`.
- `ndoubl`: Number of doubling iterations.
"""
@inline function get_dtau_ndoubl(computed_layer_properties::CoreScatteringOpticalProperties, quad_points::QuadPoints{FT}) where {FT}
    (; qp_μ) = quad_points
    (; τ, ϖ) = computed_layer_properties
    dτ_max = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])
    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    # Compute dτ vector
    dτ = τ ./ 2^ndoubl
    return dτ, ndoubl
end

"""
    get_dtau_ndoubl(computed_layer_properties::CoreDirectionalScatteringOpticalProperties, quad_points)

Variant for directional scattering: uses `G[iμ₀]` to scale optical depth for the solar beam
direction when computing `dτ_max` for the doubling criterion.
"""
@inline function get_dtau_ndoubl(computed_layer_properties::CoreDirectionalScatteringOpticalProperties, quad_points::QuadPoints{FT}) where {FT}
    (; qp_μ, iμ₀) = quad_points
    (; τ, ϖ, G) = computed_layer_properties
    gfct = Array(G)[iμ₀]
    dτ_max = minimum([maximum(gfct * τ .* ϖ), FT(0.001) * minimum(qp_μ)])
    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    # Compute dτ vector
    dτ = τ ./ 2^ndoubl
    return dτ, ndoubl
end

"""
    init_layer(computed_layer_properties, quad_points, pol_type, architecture)

Initialize layer quantities for the doubling step: elemental `dτ`, doubling count `ndoubl`,
and beam attenuation factor ``e^{-d\\tau/\\mu_0}`` (or ``e^{-d\\tau \\cdot G/\\mu_0}`` for directional scattering).

# Returns
- `dτ`: Elemental optical depth.
- `ndoubl`: Number of doubling iterations.
- `expk`: Beam attenuation factor ``e^{-d\\tau/\\mu_0}`` `[nSpec]`, on the correct architecture.
"""
@inline function init_layer(computed_layer_properties::CoreDirectionalScatteringOpticalProperties, quad_points, pol_type, architecture)
    arr_type = array_type(architecture) 
    (; μ₀, iμ₀) = quad_points
    (; G) = computed_layer_properties
    dτ, ndoubl = get_dtau_ndoubl(computed_layer_properties, quad_points)
    gfct = Array(G)[iμ₀]
    expk = exp.(-dτ*gfct/μ₀)
    return dτ, ndoubl, arr_type(expk)
end

@inline function init_layer(computed_layer_properties::CoreScatteringOpticalProperties, quad_points, pol_type, architecture)
    arr_type = array_type(architecture)
    (; μ₀) = quad_points
    dτ, ndoubl = get_dtau_ndoubl(computed_layer_properties, quad_points)
    expk = exp.(-dτ/μ₀)
    return dτ, ndoubl, arr_type(expk)
end


function rt_kernel!(RS_type::Union{RRS{FT}, VS_0to1{FT}, VS_1to0{FT}}, pol_type, SFI, added_layer, composite_layer, computed_layer_properties::CoreScatteringOpticalProperties, scattering_interface, τ_sum, m, quad_points, I_static, architecture, qp_μN, iz;
                    workspace::Union{InteractionWorkspace, Nothing}=nothing)  where {FT}
    (; qp_μ, μ₀) = quad_points
    # Just unpack core optical properties from
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = computed_layer_properties
    dτ_max = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])
    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    scatter = true # edit later
    arr_type = array_type(architecture)
    dτ = τ ./ 2^ndoubl
    expk = arr_type(exp.(-dτ /μ₀))

    (; Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀) = RS_type
    if scatter
        @timeit "elemental_inelastic" elemental_inelastic!(RS_type,
                                                pol_type, SFI,
                                                τ_sum, dτ, ϖ,
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀,
                                                RS_type.F₀,
                                                m, ndoubl, scatter,
                                                quad_points,  added_layer,
                                                I_static, architecture)
        @timeit "elemental" elemental!(pol_type, SFI, τ_sum, dτ, RS_type.F₀, computed_layer_properties, m, ndoubl, scatter, quad_points,  added_layer,  architecture)
        @timeit "doubling_inelastic" doubling_inelastic!(RS_type, pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    else
        zero_added_noscat_ie!(added_layer, τ_λ, qp_μN)
    end

    if (iz == 1)
        copy_added_to_composite_ie!(composite_layer, added_layer)
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static;
                                            workspace=workspace)
    end
end



function rt_kernel!(
            RS_type::Union{RRS_plus{FT}, VS_0to1_plus{FT}, VS_1to0_plus{FT}},
            pol_type, SFI,
            added_layer,
            composite_layer,
            computed_layer_properties::CoreScatteringOpticalProperties,
            scattering_interface,
            τ_sum, m, quad_points,
            I_static, architecture, qp_μN, iz;
            workspace::Union{InteractionWorkspace, Nothing}=nothing)  where {FT}
    (; qp_μ, μ₀) = quad_points
    # Just unpack core optical properties from 
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = computed_layer_properties
    # SUNITI, check? Also, better to write function here
    dτ_max = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])
    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    scatter = true # edit later
    arr_type = array_type(architecture)
    # Compute dτ vector
    dτ = τ ./ 2^ndoubl
    expk = arr_type(exp.(-dτ /μ₀))

    (; Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀) = RS_type
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        #@show τ, ϖ, RS_type.fscattRayl
        
        @timeit "elemental_inelastic" elemental_inelastic!(RS_type, 
                                                pol_type, SFI, 
                                                τ_sum, dτ, ϖ, 
                                                Z⁺⁺_λ₁λ₀, Z⁻⁺_λ₁λ₀, 
                                                RS_type.F₀,
                                                m, ndoubl, scatter, 
                                                quad_points,  added_layer,
                                                I_static, architecture)
        #println("Elemental inelastic done...")
        @timeit "elemental" elemental!(pol_type, SFI, τ_sum, dτ, RS_type.F₀, computed_layer_properties, m, ndoubl, scatter, quad_points,  added_layer,  architecture)
        #println("Elemental  done...")
        @timeit "doubling_inelastic" doubling_inelastic!(RS_type, pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
        #println("Doubling done...")
        #@timeit "doubling"   doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        zero_added_noscat_ie!(added_layer, τ, qp_μN)
    end

    # @assert !any(isnan.(added_layer.t⁺⁺))

    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        copy_added_to_composite_ie!(composite_layer, added_layer)

    # If this is not the TOA, perform the interaction step
    else
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, composite_layer, added_layer, I_static;
                                            workspace=workspace)
    end
end

