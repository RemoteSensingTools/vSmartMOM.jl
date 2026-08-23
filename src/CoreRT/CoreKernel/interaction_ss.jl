#=

This file contains RT interaction-related functions

=#
"""
    get_interaction_ss!(τ_sum, τ_λ, qp_μN, j₀⁺, j₀⁻, J₀⁺, J₀⁻)

KernelAbstractions single-scattering source interaction kernel. Each workitem
owns one stream/spectral source element, attenuates the existing composite
downwelling source through the added layer, and adds the new upwelling source
with attenuation through the optical depth above the added layer.
"""
@kernel function get_interaction_ss!(@Const(τ_sum), @Const(τ_λ), @Const(qp_μN),
                @Const(j₀⁺), @Const(j₀⁻), J₀⁺, J₀⁻)
    i, _, n = @index(Global, NTuple)
    J₀⁺[i,1,n] = J₀⁺[i,1,n] * exp(-τ_λ[n]/qp_μN[i]) + j₀⁺[i,1,n]
    J₀⁻[i,1,n] = J₀⁻[i,1,n] + j₀⁻[i,1,n] * exp(-τ_sum[n]/qp_μN[i])
end

function interaction_ss!(SFI::Bool,
            composite_layer::Union{CompositeLayer{FT},CompositeLayerRS{FT}},
            added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}},
            τ_sum::AbstractArray,
            τ_λ::AbstractVector,
            quad_points::QuadPoints{FT2},
            architecture) where {FT<:Real, FT2}

    (; qp_μN) = quad_points
    arr_type = array_type(architecture)
    device = devi(architecture)
    # Optical construction can legitimately promote τ to a wider type than
    # the operator workspace (for example Float64 Mie optics in a Float32 RT
    # model).  The interaction kernel writes into FT composite arrays, so cast
    # all attenuation inputs to that authoritative workspace type here.
    qp_μN = arr_type(FT.(qp_μN))
    τ_sum = arr_type(FT.(τ_sum))
    τ_λ = arr_type(FT.(τ_λ))

    kernel! = get_interaction_ss!(device)
    event = kernel!(τ_sum, τ_λ, qp_μN,
                    arr_type(added_layer.j₀⁺),
                    arr_type(added_layer.j₀⁻),
                    composite_layer.J₀⁺, composite_layer.J₀⁻,
                    ndrange=size(composite_layer.J₀⁻))
    synchronize_if_gpu()
end

function interaction_inelastic_ss!(RS_type::RRS,
    SFI::Bool,
    composite_layer::Union{CompositeLayer{FT},CompositeLayerRS{FT}},
    added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}},
    τ_sum::AbstractArray,
    τ_λ::AbstractVector,
    quad_points::QuadPoints{FT2},
    architecture) where {FT<:Real, FT2}

    (; i_λ₁λ₀) = RS_type
    (; qp_μN) = quad_points

    atype = array_type(architecture)
    device = devi(architecture)
    qp_μN = atype(FT.(qp_μN))
    τ_sum = atype(FT.(τ_sum))
    τ_λ = atype(FT.(τ_λ))
    kernel! = get_interaction_ss_RRS!(device)
    event = kernel!(τ_sum, τ_λ, qp_μN, atype(i_λ₁λ₀),
                atype(added_layer.ieJ₀⁺), atype(added_layer.ieJ₀⁻),
                composite_layer.ieJ₀⁺, composite_layer.ieJ₀⁻,
                ndrange=getKernelDimSFI(RS_type, composite_layer.ieJ₀⁻))
    synchronize_if_gpu()
end

"""
    get_interaction_ss_RRS!(τ_sum, τ_λ, qp_μN, i_λ₁λ₀, iej₀⁺, iej₀⁻, ieJ₀⁺, ieJ₀⁻)

KernelAbstractions single-scattering interaction kernel for rotational Raman
spectral coupling. Each workitem maps an output wavelength `n₁` and Raman
offset `Δn` to its incident wavelength, skips out-of-band couplings, then
updates the inelastic source vectors with the same attenuation convention as
the elastic single-scattering interaction.
"""
@kernel function get_interaction_ss_RRS!(@Const(τ_sum), @Const(τ_λ),
                    @Const(qp_μN),
                    @Const(i_λ₁λ₀),
                    @Const(iej₀⁺), @Const(iej₀⁻), ieJ₀⁺, ieJ₀⁻)
    i, _, n₁, Δn = @index(Global, NTuple)
    n₀  = n₁ + i_λ₁λ₀[Δn]
    nMax = length(τ_λ)
    if (1 <= n₀) & (n₀ <= nMax)
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
    τ_λ::AbstractVector,
    quad_points::QuadPoints{FT2},
    architecture) where {FT<:Real, FT2}

(; i_λ₁λ₀_all) = RS_type
(; qp_μN) = quad_points

atype = array_type(architecture)
device = devi(architecture)
qp_μN = atype(FT.(qp_μN))
τ_sum = atype(FT.(τ_sum))
τ_λ = atype(FT.(τ_λ))

kernel! = get_interaction_ss_VS!(device)
event = kernel!(τ_sum, τ_λ, qp_μN, atype(i_λ₁λ₀_all),
            atype(added_layer.ieJ₀⁺), atype(added_layer.ieJ₀⁻),
            composite_layer.ieJ₀⁺, composite_layer.ieJ₀⁻,
            ndrange = getKernelDimSFI(RS_type, composite_layer.ieJ₀⁻, RS_type.i_λ₁λ₀_all))
synchronize_if_gpu()
end

"""
    get_interaction_ss_VS!(τ_sum, τ_λ, qp_μN, i_λ₁λ₀_all, iej₀⁺, iej₀⁻, ieJ₀⁺, ieJ₀⁻)

KernelAbstractions single-scattering interaction kernel for vibrational Raman
couplings represented by a sparse wavelength map. Each workitem owns one
active Raman offset, skips inactive map entries, and updates the inelastic
source vectors at the mapped wavelength.
"""
@kernel function get_interaction_ss_VS!(
                    @Const(τ_sum), @Const(τ_λ),
                    @Const(qp_μN),
                    @Const(i_λ₁λ₀_all),
                    @Const(iej₀⁺), @Const(iej₀⁻), ieJ₀⁺, ieJ₀⁻)
    i, Δn = @index(Global, NTuple)
    n₁ =  i_λ₁λ₀_all[Δn]

    if (n₁ > 0)
        ieJ₀⁺[i,1,n₁,1] = ieJ₀⁺[i,1,n₁,1] * exp(-τ_λ[n₁]/qp_μN[i]) + iej₀⁺[i,1,n₁,1]
        ieJ₀⁻[i,1,n₁,1] = ieJ₀⁻[i,1,n₁,1] + iej₀⁻[i,1,n₁,1] * exp(-τ_sum[n₁]/qp_μN[i])
    end
end
