#=
 
This file contains RT interaction-related functions
 
=#
# Scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer which is added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_hdrf!(SFI,
    composite_layer::CompositeLayer{FT}, 
    added_layer::AddedLayer{FT}, 
    m, pol_type, quad_points,
    hdr_J₀⁻, bhr_J₀⁻, bhr_J₀⁺) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    @unpack r⁻⁺, j₀⁻, j₀⁺ = added_layer     #these are aliases to the respective struct elements  
    @unpack J₀⁺, J₀⁻      = composite_layer #these are aliases to the respective struct elements 
    @unpack Nquad, wt_μN, iμ₀, qp_μN = quad_points
    NquadN =  Nquad * pol_type.n
    hdr_J₀⁻ .= r⁻⁺ ⊠ J₀⁺ .+ j₀⁻
    # @show hdr_J₀⁻./ J₀⁺
    
    qp = Array(qp_μN)
    if m==0

        for i = 1:pol_type.n
            #bhr_J₀⁻[i,:] .= 0
            #bhr_J₀⁺[i,:] .= 0
            #for j=i:pol_type.n:NquadN
            j=i:pol_type.n:NquadN
            bhr_J₀⁻[i,:] .= Array(sum(hdr_J₀⁻[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)')
            bhr_J₀⁺[i,:] .= Array(sum(J₀⁺[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)' .+ j₀⁺[iμ₀,1,:] .* qp[iμ₀]) #TODO: Use Radau quadrature and include insolation in the quadrature sum
            #@show j₀⁺[iμ₀,1,1:3].* qp[iμ₀], J₀⁺[iμ₀,1,1:3] .* qp[iμ₀], bhr_J₀⁺[i,1:3], bhr_J₀⁻[i,1:3]
            #@show Array(sum(J₀⁻[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)')[1], Array(sum(hdr_J₀⁻[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)')[1],Array(sum(j₀⁻[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)')[1]
            #end
        end
    end
end

function interaction_hdrf_canopy!(SFI,
    #composite_layer::CompositeLayer{FT}, 
    dwJ, uwJ, solJ₀,
    m, pol_type, quad_points,
    hdr_J₀⁻, bhr_J₀⁻, bhr_J₀⁺) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

    #@unpack topJ₀⁺, botJ₀⁻      = composite_layer #these are aliases to the respective struct elements 
    @unpack Nquad, wt_μN, iμ₀, qp_μN = quad_points
    NquadN =  Nquad * pol_type.n
    
    wt = Array(wt_μN)
    qp = Array(qp_μN)
    hdr_J₀⁻ .= uwJ#r⁻⁺ ⊠ J₀⁺ .+ j₀⁻
    # @show hdr_J₀⁻./ J₀⁺
    
    #@show size(bhr_J₀⁻), size(bhr_J₀⁺)
    #@show size(uwJ), size(dwJ)
    #qp = Array(qp_μN)
    if m==0
        
        @show dwJ[:,1,end]
        @show "========="
        @show uwJ[:,1,end]
        #showp
        for i = 1:pol_type.n
            #bhr_J₀⁻[i,:] .= 0
            #bhr_J₀⁺[i,:] .= 0
            #for j=i:pol_type.n:NquadN
            j=i:pol_type.n:NquadN
            #@show typeof(uwJ[j,1,:]),  typeof(wt_μN[j])
            bhr_J₀⁻[i,:] .= (sum(uwJ[j,1,:].*wt[j].*qp[j], dims=1)')
            #@show size(bhr_J₀⁺[i,:]), size(solJ₀[i,:])
            bhr_J₀⁺[i,:] .= (sum(dwJ[j,1,:].*wt[j].*qp[j], dims=1)') .+ (solJ₀[i,:] .* qp[iμ₀]) #TODO: Use Radau quadrature and include insolation in the quadrature sum
            #@show j₀⁺[iμ₀,1,1:3].* qp[iμ₀], J₀⁺[iμ₀,1,1:3] .* qp[iμ₀], bhr_J₀⁺[i,1:3], bhr_J₀⁻[i,1:3]
            #@show Array(sum(J₀⁻[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)')[1], Array(sum(hdr_J₀⁻[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)')[1],Array(sum(j₀⁻[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)')[1]
            #end
        end
        @show bhr_J₀⁻./bhr_J₀⁺
        @show solJ₀[1,:]
    end

end

"Compute interaction between composite and added layers"
#==function interaction_hdrf!(SFI,
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
==#