function interlayer_flux_helper!(RS_type::noRS, 
    I_static::AbstractArray{FT2},
    itopR⁺⁻::AbstractArray{FT}, ibotR⁻⁺::AbstractArray{FT},
    itopJ₀⁺::AbstractArray{FT}, ibotJ₀⁻::AbstractArray{FT},
    otdwJ::AbstractArray{FT}, otuwJ::AbstractArray{FT}) where {FT<:Real,FT2}

    tmpR = similar(itopR⁺⁻)
    # elastic
    @timeit "interlayer_flux inv1" batch_inv!(
        tmpR, I_static .- itopR⁺⁻ ⊠ ibotR⁻⁺) #Suniti
     
    otdwJ .= tmpR ⊠ (itopJ₀⁺ .+ itopR⁺⁻ ⊠ ibotJ₀⁻)

    @timeit "interlayer_flux inv2" batch_inv!(
        tmpR, I_static .- ibotR⁻⁺ ⊠ itopR⁺⁻) #Suniti
    
    otuwJ .= tmpR ⊠ (ibotJ₀⁻ .+ ibotR⁻⁺ ⊠ itopJ₀⁺)
end

"""
    interlayer_flux_helper!(::noRS, I_static,
        topR⁺⁻, botR⁻⁺, topJ₀⁺, botJ₀⁻, dwJ, uwJ,
        topṘ⁺⁻, botṘ⁻⁺, topJ̇₀⁺, botJ̇₀⁻, dwJ̇, uwJ̇)

Compute the elastic interlayer downwelling/upwelling fields and their
Jacobians. The fourth dimension of every dotted array is the physical
parameter dimension.

For the downwelling field, with `A = topR⁺⁻`, `B = botR⁻⁺`, and
`G = (I - A·B)⁻¹`, the tangent of the inverse is evaluated analytically as

```math
\\dot G = G\\,(\\dot A\\,B + A\\,\\dot B)\\,G.
```

The mirror-direction solve uses the corresponding `(I - B·A)⁻¹` identity.
Both forward inverses are therefore computed once and reused for every
parameter column; no AD is applied through `batch_inv!`.
"""
function interlayer_flux_helper!(RS_type::noRS,
    I_static::AbstractArray{FT2},
    itopR⁺⁻::AbstractArray{FT,3}, ibotR⁻⁺::AbstractArray{FT,3},
    itopJ₀⁺::AbstractArray{FT,3}, ibotJ₀⁻::AbstractArray{FT,3},
    otdwJ::AbstractArray{FT,3}, otuwJ::AbstractArray{FT,3},
    itopṘ⁺⁻::AbstractArray{FT,4}, ibotṘ⁻⁺::AbstractArray{FT,4},
    itopJ̇₀⁺::AbstractArray{FT,4}, ibotJ̇₀⁻::AbstractArray{FT,4},
    otdwJ̇::AbstractArray{FT,4}, otuwJ̇::AbstractArray{FT,4}) where {FT<:Real,FT2}

    Nparams = size(itopṘ⁺⁻, 4)
    size(ibotṘ⁻⁺, 4) == Nparams || throw(DimensionMismatch(
        "top and bottom reflectance Jacobians must have the same parameter dimension"))
    size(itopJ̇₀⁺, 4) == Nparams || throw(DimensionMismatch(
        "top source and reflectance Jacobians must have the same parameter dimension"))
    size(ibotJ̇₀⁻, 4) == Nparams || throw(DimensionMismatch(
        "bottom source and reflectance Jacobians must have the same parameter dimension"))
    size(otdwJ̇, 4) == Nparams || throw(DimensionMismatch(
        "downwelling output Jacobian has the wrong parameter dimension"))
    size(otuwJ̇, 4) == Nparams || throw(DimensionMismatch(
        "upwelling output Jacobian has the wrong parameter dimension"))

    # Allocate one inverse-tangent and one source-tangent workspace, then reuse
    # them for every parameter and for both propagation directions. All arrays
    # remain on the architecture selected by the caller.
    G = similar(itopR⁺⁻)
    Ġ = similar(itopR⁺⁻)
    rhs = similar(itopJ₀⁺)
    rhṡ = similar(itopJ₀⁺)

    # Downwelling: D = (I - A B)⁻¹ (a + A b).
    batch_inv!(
        G, I_static .- itopR⁺⁻ ⊠ ibotR⁻⁺)
    rhs .= itopJ₀⁺ .+ itopR⁺⁻ ⊠ ibotJ₀⁻
    otdwJ .= G ⊠ rhs

    @inbounds for iparam in 1:Nparams
        @views begin
            Ġ .= G ⊠ (itopṘ⁺⁻[:, :, :, iparam] ⊠ ibotR⁻⁺ .+
                         itopR⁺⁻ ⊠ ibotṘ⁻⁺[:, :, :, iparam]) ⊠ G
            rhṡ .= itopJ̇₀⁺[:, :, :, iparam] .+
                    itopṘ⁺⁻[:, :, :, iparam] ⊠ ibotJ₀⁻ .+
                    itopR⁺⁻ ⊠ ibotJ̇₀⁻[:, :, :, iparam]
            otdwJ̇[:, :, :, iparam] .= Ġ ⊠ rhs .+ G ⊠ rhṡ
        end
    end

    # Upwelling: U = (I - B A)⁻¹ (b + B a).
    batch_inv!(
        G, I_static .- ibotR⁻⁺ ⊠ itopR⁺⁻)
    rhs .= ibotJ₀⁻ .+ ibotR⁻⁺ ⊠ itopJ₀⁺
    otuwJ .= G ⊠ rhs

    @inbounds for iparam in 1:Nparams
        @views begin
            Ġ .= G ⊠ (ibotṘ⁻⁺[:, :, :, iparam] ⊠ itopR⁺⁻ .+
                         ibotR⁻⁺ ⊠ itopṘ⁺⁻[:, :, :, iparam]) ⊠ G
            rhṡ .= ibotJ̇₀⁻[:, :, :, iparam] .+
                    ibotṘ⁻⁺[:, :, :, iparam] ⊠ itopJ₀⁺ .+
                    ibotR⁻⁺ ⊠ itopJ̇₀⁺[:, :, :, iparam]
            otuwJ̇[:, :, :, iparam] .= Ġ ⊠ rhs .+ G ⊠ rhṡ
        end
    end

    return nothing
end


function interlayer_flux_helper!(RS_type::RRS, 
        I_static::AbstractArray{FT2},
        itopR⁺⁻::AbstractArray{FT}, ibotR⁻⁺::AbstractArray{FT},
        itopJ₀⁺::AbstractArray{FT}, ibotJ₀⁻::AbstractArray{FT},
        otdwJ::AbstractArray{FT}, otuwJ::AbstractArray{FT},
        itopieR⁺⁻::AbstractArray{FT}, ibotieR⁻⁺::AbstractArray{FT},
        itopieJ₀⁺::AbstractArray{FT}, ibotieJ₀⁻::AbstractArray{FT},
        otdwieJ::AbstractArray{FT}, otuwieJ::AbstractArray{FT}) where {FT<:Real,FT2}
    (; i_λ₁λ₀) = RS_type
    tmpR = similar(itopR⁺⁻)
    # elastic
    #@show size(itopR⁺⁻)
    @timeit "interlayer_flux inv1" batch_inv!(
        tmpR, I_static .- itopR⁺⁻ ⊠ ibotR⁻⁺) #Suniti
    otdwJ[:] = tmpR ⊠ (itopJ₀⁺ .+ itopR⁺⁻ ⊠ ibotJ₀⁻)
    # inelastic
    #RRS
    (; i_λ₁λ₀) = RS_type
    for Δn=1:size(itopieJ₀⁺,4)
    #for n₁ = 1:size(itopieJ₀⁺,3)
        #for Δn=1:size(itopieJ₀⁺,4) #eachindex itopieJ₀⁺[1,1,1,:]
            #n₀  = n₁ + i_λ₁λ₀[Δn]
            n₀, n₁ = get_n₀_n₁(itopieR⁺⁻,i_λ₁λ₀[Δn])    
            #for t =1:size(topieJ₀⁺,5)
            @inbounds @views otdwieJ[:,:,n₁,Δn] = tmpR[:,:,n₁] ⊠ (
                        itopieJ₀⁺[:,:,n₁,Δn] + 
                        itopieR⁺⁻[:,:,n₁,Δn] ⊠ ibotJ₀⁻[:,:,n₀] +
                        itopR⁺⁻[:,:,n₁] ⊠ ibotieJ₀⁻[:,:,n₁,Δn] +
                        (itopR⁺⁻[:,:,n₁] ⊠ ibotieR⁻⁺[:,:,n₁,Δn] + 
                        itopieR⁺⁻[:,:,n₁,Δn] ⊠ ibotR⁻⁺[:,:,n₀]) ⊠ 
                        tmpR[:,:,n₀] ⊠
                        (itopJ₀⁺[:,:,n₀] + itopR⁺⁻[:,:,n₀] ⊠ ibotJ₀⁻[:,:,n₀]))
        #end
    end
    
    # elastic
    @timeit "interlayer_flux inv2" batch_inv!(
        tmpR, I_static .- ibotR⁻⁺ ⊠ itopR⁺⁻) #Suniti
    otuwJ[:] = tmpR ⊠ (ibotJ₀⁻ .+ ibotR⁻⁺ ⊠ itopJ₀⁺)

    # inelastic
    for Δn=1:size(itopieJ₀⁺,4)
    #for n₁ = 1:size(itopieJ₀⁺,3)
        #for Δn=1:size(itopieJ₀⁺,4) #eachindex itopieJ₀⁺[1,1,1,:]
            n₀, n₁ = get_n₀_n₁(itopieR⁺⁻,i_λ₁λ₀[Δn])     
            #for t =1:size(topieJ₀⁺,5)
            @inbounds @views otuwieJ[:,:,n₁,Δn] = tmpR[:,:,n₁] ⊠ (
                        ibotieJ₀⁻[:,:,n₁,Δn] + 
                        ibotieR⁻⁺[:,:,n₁,Δn] ⊠ itopJ₀⁺[:,:,n₀] +
                        ibotR⁻⁺[:,:,n₁] ⊠ itopieJ₀⁺[:,:,n₁,Δn] +
                        (ibotR⁻⁺[:,:,n₁] ⊠ itopieR⁺⁻[:,:,n₁,Δn] + 
                        ibotieR⁻⁺[:,:,n₁,Δn] ⊠ itopR⁺⁻[:,:,n₀]) ⊠ 
                        tmpR[:,:,n₀] ⊠
                        (ibotJ₀⁻[:,:,n₀] + ibotR⁻⁺[:,:,n₀] ⊠ itopJ₀⁺[:,:,n₀]))
        #end
    end
end

function interlayer_flux_helper!(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
        I_static::AbstractArray{FT2},
        itopR⁺⁻::AbstractArray{FT}, ibotR⁻⁺::AbstractArray{FT},
        itopJ₀⁺::AbstractArray{FT}, ibotJ₀⁻::AbstractArray{FT},
        otdwJ::AbstractArray{FT}, otuwJ::AbstractArray{FT},
        itopieR⁺⁻::AbstractArray{FT}, ibotieR⁻⁺::AbstractArray{FT},
        itopieJ₀⁺::AbstractArray{FT}, ibotieJ₀⁻::AbstractArray{FT},
        otdwieJ::AbstractArray{FT}, otuwieJ::AbstractArray{FT}) where {FT<:Real,FT2}
    
    (; i_λ₁λ₀_all) = RS_type
    
    tmpR = similar(itopR⁺⁻)
    # elastic
    @timeit "interlayer_flux inv1" batch_inv!(
        tmpR, I_static .- itopR⁺⁻ ⊠ ibotR⁻⁺) #Suniti
    otdwJ[:] = tmpR ⊠ (itopJ₀⁺ .+ itopR⁺⁻ ⊠ ibotJ₀⁻)
    # inelastic
    for Δn = 1:length(i_λ₁λ₀_all) # in eachindex ieJ₁⁺[1,1,:,1]
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        if (n₁>0)
            otdwieJ[:,:,n₁,n₀] = tmpR[:,:,n₁] * (
                                    itopieJ₀⁺[:,1,n₁,n₀] + 
                                    itopieR⁺⁻[:,:,n₁,n₀] * ibotJ₀⁻[:,1,n₀] +
                                    itopR⁺⁻[:,:,n₁] * ibotieJ₀⁻[:,1,n₁,n₀] +
                                    (itopR⁺⁻[:,:,n₁] * ibotieR⁻⁺[:,:,n₁,n₀] + 
                                    itopieR⁺⁻[:,:,n₁,n₀] * ibotR⁻⁺[:,:,n₀]) * 
                                    tmpR[:,:,n₀] *
                                    (itopJ₀⁺[:,1,n₀] + itopR⁺⁻[:,:,n₀] * ibotJ₀⁻[:,1,n₀]))
        end
    end

    @timeit "interlayer_flux inv2" batch_inv!(
        tmpR, I_static .- ibotR⁻⁺ ⊠ itopR⁺⁻) #Suniti
    otuwJ[:] = tmpR ⊠ (ibotJ₀⁻ .+ ibotR⁻⁺ ⊠ itopJ₀⁺)

    # inelastic
    for Δn = 1:length(i_λ₁λ₀_all) # in eachindex ieJ₁⁺[1,1,:,1]
        n₁ = i_λ₁λ₀_all[Δn]
        n₀ = 1
        if (n₁>0)
            otuwieJ[:,:,n₁,n₀] = tmpR[:,:,n₁] * (
                                    ibotieJ₀⁻[:,1,n₁,n₀] + 
                                    ibotieR⁻⁺[:,:,n₁,n₀] * itopJ₀⁺[:,1,n₀] +
                                    ibotR⁻⁺[:,:,n₁] * itopieJ₀⁺[:,1,n₁,n₀] +
                                    (ibotR⁻⁺[:,:,n₁] * itopieR⁺⁻[:,:,n₁,n₀] + 
                                    ibotieR⁻⁺[:,:,n₁,n₀] * itopR⁺⁻[:,:,n₀]) * 
                                    tmpR[:,:,n₀] *
                                    (ibotJ₀⁻[:,1,n₀] .+ ibotR⁻⁺[:,:,n₀] * itopJ₀⁺[:,1,n₀]))
        end
    end
end

"Compute interaction between composite and added layers with inelastic scattering"
function compute_interlayer_flux!(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus}, 
                        #scattering_interface::AbstractScatteringInterface, SFI,
                        #composite_layer::Union{CompositeLayer,CompositeLayerRS}, 
                        #added_layer::Union{AddedLayer,AddedLayerRS},
                        I_static::AbstractArray{FT2},
                        itopR⁺⁻::AbstractArray{FT}, ibotR⁻⁺::AbstractArray{FT}, 
                        itopJ₀⁺::AbstractArray{FT}, ibotJ₀⁻::AbstractArray{FT},
                        otdwJ::AbstractArray{FT}, otuwJ::AbstractArray{FT},
                        itopieR⁺⁻::AbstractArray{FT}, ibotieR⁻⁺::AbstractArray{FT}, 
                        itopieJ₀⁺::AbstractArray{FT}, ibotieJ₀⁻::AbstractArray{FT},
                        otdwieJ::AbstractArray{FT}, 
                        otuwieJ::AbstractArray{FT},
                        arr_type) where {FT<:Real,FT2}
    topR⁺⁻ = arr_type(itopR⁺⁻) 
    botR⁻⁺ = arr_type(ibotR⁻⁺)
    topJ₀⁺ = arr_type(itopJ₀⁺) 
    botJ₀⁻ = arr_type(ibotJ₀⁻)
    dwJ    = arr_type(otdwJ)
    uwJ    = arr_type(otuwJ)
    topieR⁺⁻ = arr_type(itopieR⁺⁻) 
    botieR⁻⁺ = arr_type(ibotieR⁻⁺)
    topieJ₀⁺ = arr_type(itopieJ₀⁺) 
    botieJ₀⁻ = arr_type(ibotieJ₀⁻)
    dwieJ    = arr_type(otdwieJ)
    uwieJ    = arr_type(otuwieJ)
    interlayer_flux_helper!(RS_type, I_static,
        topR⁺⁻, botR⁻⁺,
        topJ₀⁺, botJ₀⁻,
        dwJ, uwJ,
        topieR⁺⁻, botieR⁻⁺,
        topieJ₀⁺, botieJ₀⁻,
        dwieJ, uwieJ)

    itopR⁺⁻.= collect(topR⁺⁻) 
    ibotR⁻⁺ .= collect(botR⁻⁺)
    itopJ₀⁺ .= collect(topJ₀⁺) 
    ibotJ₀⁻ .= collect(botJ₀⁻)
    otdwJ    .= collect(dwJ)
    otuwJ    .= collect(uwJ)
    itopieR⁺⁻ .= collect(topieR⁺⁻) 
    ibotieR⁻⁺ .= collect(botieR⁻⁺)
    itopieJ₀⁺ .= collect(topieJ₀⁺) 
    ibotieJ₀⁻ .= collect(botieJ₀⁻)
    otdwieJ    .= collect(dwieJ)
    otuwieJ    .= collect(uwieJ)
    
    #scattering_interface, SFI, composite_layer, added_layer, I_static)
    synchronize_if_gpu()
    
end

"Compute interaction between composite and added layers with inelastic scattering"
function compute_interlayer_flux!(RS_type::noRS, 
                        #scattering_interface::AbstractScatteringInterface, SFI,
                        #composite_layer::Union{CompositeLayer,CompositeLayerRS}, 
                        #added_layer::Union{AddedLayer,AddedLayerRS},
                        I_static::AbstractArray{FT2},
                        itopR⁺⁻::AbstractArray{FT}, ibotR⁻⁺::AbstractArray{FT}, 
                        itopJ₀⁺::AbstractArray{FT}, ibotJ₀⁻::AbstractArray{FT},
                        otdwJ::AbstractArray{FT}, otuwJ::AbstractArray{FT}, 
                        arr_type) where {FT<:Real,FT2}

    interlayer_flux_helper!(RS_type, I_static,
        itopR⁺⁻, ibotR⁻⁺,
        itopJ₀⁺, ibotJ₀⁻,
        otdwJ, otuwJ)
    
    synchronize_if_gpu()
    
end

"""Compute elastic interlayer fields and their analytic parameter Jacobians."""
function compute_interlayer_flux!(RS_type::noRS,
                        I_static::AbstractArray{FT2},
                        itopR⁺⁻::AbstractArray{FT,3}, ibotR⁻⁺::AbstractArray{FT,3},
                        itopJ₀⁺::AbstractArray{FT,3}, ibotJ₀⁻::AbstractArray{FT,3},
                        otdwJ::AbstractArray{FT,3}, otuwJ::AbstractArray{FT,3},
                        itopṘ⁺⁻::AbstractArray{FT,4}, ibotṘ⁻⁺::AbstractArray{FT,4},
                        itopJ̇₀⁺::AbstractArray{FT,4}, ibotJ̇₀⁻::AbstractArray{FT,4},
                        otdwJ̇::AbstractArray{FT,4}, otuwJ̇::AbstractArray{FT,4},
                        arr_type) where {FT<:Real,FT2}

    interlayer_flux_helper!(RS_type, I_static,
        itopR⁺⁻, ibotR⁻⁺,
        itopJ₀⁺, ibotJ₀⁻,
        otdwJ, otuwJ,
        itopṘ⁺⁻, ibotṘ⁻⁺,
        itopJ̇₀⁺, ibotJ̇₀⁻,
        otdwJ̇, otuwJ̇)

    synchronize_if_gpu()
    return nothing
end
