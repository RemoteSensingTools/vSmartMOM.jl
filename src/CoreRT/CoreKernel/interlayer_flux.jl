function interlayer_flux_helper!(RS_type::noRS, 
    I_static::AbstractArray{FT2},
    itopR⁺⁻::AbstractArray{FT}, ibotR⁻⁺::AbstractArray{FT},
    itopJ₀⁺::AbstractArray{FT}, ibotJ₀⁻::AbstractArray{FT},
    otdwJ::AbstractArray{FT}, otuwJ::AbstractArray{FT}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    tmpR = similar(itopR⁺⁻)
    # elastic
    @timeit "interlayer_flux inv1" batch_inv!(
        tmpR, I_static .- itopR⁺⁻ ⊠ ibotR⁻⁺) #Suniti
     
    otdwJ .= tmpR ⊠ (itopJ₀⁺ .+ itopR⁺⁻ ⊠ ibotJ₀⁻)

    @timeit "interlayer_flux inv2" batch_inv!(
        tmpR, I_static .- ibotR⁻⁺ ⊠ itopR⁺⁻) #Suniti
    
    otuwJ .= tmpR ⊠ (ibotJ₀⁻ .+ ibotR⁻⁺ ⊠ itopJ₀⁺)
end


function interlayer_flux_helper!(RS_type::RRS, 
        I_static::AbstractArray{FT2},
        itopR⁺⁻::AbstractArray{FT}, ibotR⁻⁺::AbstractArray{FT},
        itopJ₀⁺::AbstractArray{FT}, ibotJ₀⁻::AbstractArray{FT},
        otdwJ::AbstractArray{FT}, otuwJ::AbstractArray{FT},
        itopieR⁺⁻::AbstractArray{FT}, ibotieR⁻⁺::AbstractArray{FT},
        itopieJ₀⁺::AbstractArray{FT}, ibotieJ₀⁻::AbstractArray{FT},
        otdwieJ::AbstractArray{FT}, otuwieJ::AbstractArray{FT}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    @unpack i_λ₁λ₀ = RS_type
    tmpR = similar(itopR⁺⁻)
    # elastic
    #@show size(itopR⁺⁻)
    @timeit "interlayer_flux inv1" batch_inv!(
        tmpR, I_static .- itopR⁺⁻ ⊠ ibotR⁻⁺) #Suniti
    otdwJ[:] = tmpR ⊠ (itopJ₀⁺ .+ itopR⁺⁻ ⊠ ibotJ₀⁻)
    # inelastic
    #RRS
    @unpack i_λ₁λ₀ = RS_type
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
        otdwieJ::AbstractArray{FT}, otuwieJ::AbstractArray{FT}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    
    @unpack i_λ₁λ₀_all = RS_type
    
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
                        arr_type) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
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

    itopR⁺⁻.= Array(topR⁺⁻) 
    ibotR⁻⁺ .= Array(botR⁻⁺)
    itopJ₀⁺ .= Array(topJ₀⁺) 
    ibotJ₀⁻ .= Array(botJ₀⁻)
    otdwJ    .= Array(dwJ)
    otuwJ    .= Array(uwJ)
    itopieR⁺⁻ .= Array(topieR⁺⁻) 
    ibotieR⁻⁺ .= Array(botieR⁻⁺)
    itopieJ₀⁺ .= Array(topieJ₀⁺) 
    ibotieJ₀⁻ .= Array(botieJ₀⁻)
    otdwieJ    .= Array(dwieJ)
    otuwieJ    .= Array(uwieJ)
    
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
                        arr_type) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    interlayer_flux_helper!(RS_type, I_static,
        itopR⁺⁻, ibotR⁻⁺,
        itopJ₀⁺, ibotJ₀⁻,
        otdwJ, otuwJ)
    
    synchronize_if_gpu()
    
end