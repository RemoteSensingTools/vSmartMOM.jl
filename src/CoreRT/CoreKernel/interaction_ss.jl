#=
 
This file contains RT interaction-related functions
 
=#
@kernel function get_interaction_ss!(τ_sum, τ_λ, qp_μN,
                j₀⁺, j₀⁻, J₀⁺, J₀⁻)
    i, _, n = @index(Global, NTuple)
    J₀⁺[i,1,n] = J₀⁺[i,1,n] * exp(-τ_λ[n]/qp_μN[i]) + j₀⁺[i,1,n]
    J₀⁻[i,1,n] = J₀⁻[i,1,n] + j₀⁻[i,1,n] * exp(-τ_sum[n]/qp_μN[i])
end

function interaction_ss!(SFI::Bool,
            composite_layer::Union{CompositeLayer{FT},CompositeLayerRS{FT}}, 
            added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}}, 
            τ_sum::AbstractArray,
            τ_λ::AbstractArray{FT,1},
            quad_points::QuadPoints{FT2},
            architecture) where {FT<:Real, FT2}
    
    #@unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer #these are aliases to the respective struct elements  
    @unpack J₀⁺, J₀⁻ = composite_layer #these are aliases to the respective struct elements 
    @unpack qp_μN = quad_points
    arr_type = array_type(architecture)
    device = devi(architecture)
    qp_μN = arr_type(qp_μN)
    τ_sum = arr_type(τ_sum)
    τ_λ = arr_type(τ_λ)
    J₀⁺ = arr_type(J₀⁺)
    J₀⁻ = arr_type(J₀⁻)

    kernel! = get_interaction_ss!(device)
    event = kernel!(τ_sum, τ_λ, qp_μN, 
                    arr_type(added_layer.j₀⁺), 
                    arr_type(added_layer.j₀⁻),
                    J₀⁺, J₀⁻, ndrange=size(J₀⁻))
    #wait(device, event)
    synchronize_if_gpu()
end

function interaction_inelastic_ss!(RS_type::RRS,
    SFI::Bool,
    composite_layer::Union{CompositeLayer{FT},CompositeLayerRS{FT}}, 
    added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}}, 
    τ_sum::AbstractArray,
    τ_λ::AbstractArray{FT,1},
    quad_points::QuadPoints{FT2},
    architecture) where {FT<:Real, FT2}

    #@unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer #these are aliases to the respective struct elements  
    @unpack i_λ₁λ₀ = RS_type
    @unpack ieJ₀⁺, ieJ₀⁻ = composite_layer #these are aliases to the respective struct elements 
    @unpack qp_μN = quad_points

    atype = array_type(architecture)
    device = devi(architecture)
    qp_μN = atype(qp_μN)
    τ_sum = atype(τ_sum)
    τ_λ = atype(τ_λ)
    ieJ₀⁺ = atype(ieJ₀⁺)
    ieJ₀⁻ = atype(ieJ₀⁻)
    aa = getKernelDimSFI(RS_type, ieJ₀⁻)
    kernel! = get_interaction_ss_RRS!(device)
    event = kernel!(τ_sum, τ_λ, qp_μN, atype(i_λ₁λ₀),
                atype(added_layer.ieJ₀⁺), atype(added_layer.ieJ₀⁻),
                ieJ₀⁺, ieJ₀⁻,
                ndrange=getKernelDimSFI(RS_type, ieJ₀⁻))
    #wait(device, event)
    synchronize_if_gpu()
end

@kernel function get_interaction_ss_RRS!(τ_sum, τ_λ, 
                    qp_μN,
                    i_λ₁λ₀,
                    iej₀⁺, iej₀⁻, ieJ₀⁺, ieJ₀⁻)
    i, _, n₁, Δn = @index(Global, NTuple)
    n₀  = n₁ + i_λ₁λ₀[Δn]
    nMax = length(τ_λ) 
    if (1 ≤ n₀ ≤ nMax) 
        ieJ₀⁺[i,1,n₁,Δn] = ieJ₀⁺[i,1,n₁,Δn] * exp(-τ_λ[n₁]/qp_μN[i]) + iej₀⁺[i,1,n₁,Δn]
        ieJ₀⁻[i,1,n₁,Δn] = ieJ₀⁻[i,1,n₁,Δn] + iej₀⁻[i,1,n₁,Δn] * exp(-τ_sum[n₁]/qp_μN[i])
    end
end

function interaction_inelastic_ss!(
    RS_type::Union{VS_0to1_plus, VS_1to0_plus},
    SFI::Bool,
    composite_layer::Union{CompositeLayer{FT},CompositeLayerRS{FT}}, 
    added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}}, 
    τ_sum::AbstractArray,
    τ_λ::AbstractArray{FT,1},
    quad_points::QuadPoints{FT2},
    architecture) where {FT<:Real, FT2}

@unpack i_λ₁λ₀_all = RS_type
@unpack ieJ₀⁺, ieJ₀⁻ = composite_layer #these are aliases to the respective struct elements 
@unpack qp_μN = quad_points

atype = array_type(architecture)
device = devi(architecture)
qp_μN = atype(qp_μN)
τ_sum = atype(τ_sum)
τ_λ = atype(τ_λ)
ieJ₀⁺ = atype(ieJ₀⁺)
ieJ₀⁻ = atype(ieJ₀⁻)

kernel! = get_interaction_ss_VS!(device)
event = kernel!(τ_sum, τ_λ, qp_μN, atype(i_λ₁λ₀_all),
            atype(added_layer.ieJ₀⁺), atype(added_layer.ieJ₀⁻),
            ieJ₀⁺, ieJ₀⁻,
            ndrange = getKernelDimSFI(RS_type,ieJ₀⁻,RS_type.i_λ₁λ₀_all))
#wait(device, event)
synchronize_if_gpu()
end

@kernel function get_interaction_ss_VS!(
                    τ_sum, τ_λ, 
                    qp_μN,
                    i_λ₁λ₀_all,
                    iej₀⁺, iej₀⁻, ieJ₀⁺, ieJ₀⁻)
    i, Δn = @index(Global, NTuple)
    n₁ =  i_λ₁λ₀_all[Δn]
     
    if (n₁ > 0) 
        ieJ₀⁺[i,1,n₁,1] = ieJ₀⁺[i,1,n₁,1] * exp(-τ_λ[n₁]/qp_μN[i]) + iej₀⁺[i,1,n₁,1]
        ieJ₀⁻[i,1,n₁,1] = ieJ₀⁻[i,1,n₁,1] + iej₀⁻[i,1,n₁,1] * exp(-τ_sum[n₁]/qp_μN[i])
    end
end