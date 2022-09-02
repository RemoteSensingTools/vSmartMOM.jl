#=
 
This file contains RT interaction-related functions
 
=#
# Scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer which is added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_hdrf_helper!(SFI,
    composite_layer::CompositeLayer{FT}, 
    added_layer::AddedLayer{FT}, 
    m, pol_type, NquadN, wt_μN,
    hdr_J₀⁻, bhr_J₀⁻, bhr_J₀⁺) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack r⁻⁺, j₀⁻ = added_layer     #these are aliases to the respective struct elements  
    @unpack J₀⁺      = composite_layer #these are aliases to the respective struct elements 

    hdr_J₀⁻ .= r⁻⁺ ⊠ J₀⁺ .+ j₀⁻

    if m==0
        for i = 1:pol_type.n
            #bhr_J₀⁻[i,:] .= 0
            #bhr_J₀⁺[i,:] .= 0
            #for j=i:pol_type.n:NquadN
            j=i:pol_type.n:NquadN
            bhr_J₀⁻[i,:] .= Array(sum(hdr_J₀⁻[j,1,:].*wt_μN[j], dims=1)')
            bhr_J₀⁺[i,:] .= Array(sum(J₀⁺[j,1,:].*wt_μN[j], dims=1)')
            #end
        end
    end
end

"Compute interaction between composite and added layers"
function interaction_hdrf!(SFI,
    composite_layer::CompositeLayer{FT}, 
    added_layer::AddedLayer{FT},
    m, pol_type, NquadN, wt_μN,
    hdr_J₀⁻::AbstractArray{FT2}, 
    bhr_J₀⁻::AbstractArray{FT2}, 
    bhr_J₀⁺::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    interaction_hdrf_helper!(SFI, 
                composite_layer, added_layer, 
                m, pol_type, NquadN, wt_μN,
                hdr_J₀⁻, bhr_J₀⁻, bhr_J₀⁺)
    
    synchronize_if_gpu()
end