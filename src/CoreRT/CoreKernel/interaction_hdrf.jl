# Scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer which is added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function interaction_hdrf!(SFI,
    composite_layer::AbstractLayer, 
    added_layer::AbstractLayer, 
    m, pol_type, quad_points,
    hdr_J₀⁻, bhr_J₀⁻, bhr_J₀⁺)

    r⁻⁺ = added_layer.r⁻⁺
    j₀⁻ = added_layer.j₀⁻
    j₀⁺ = added_layer.j₀⁺
    @unpack J₀⁺, J₀⁻ = composite_layer 
    @unpack Nquad, wt_μN, iμ₀, iμ₀Nstart, qp_μN = quad_points
    NquadN =  Nquad * pol_type.n
    hdr_J₀⁻ .= r⁻⁺ ⊠ J₀⁺ .+ j₀⁻
    # @show hdr_J₀⁻./ J₀⁺
    #@show hdr_J₀⁻[1,1,:]
    #@show J₀⁺[1,1,:]
    #@show iμ₀
    #@show j₀⁺[iμ₀Nstart,1,:]
    qp = collect(qp_μN)
    if m==0

        for i = 1:pol_type.n
            #bhr_J₀⁻[i,:] .= 0
            #bhr_J₀⁺[i,:] .= 0
            #for j=i:pol_type.n:NquadN
            j=i:pol_type.n:NquadN
            bhr_J₀⁻[i,:] .= collect(sum(hdr_J₀⁻[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)')
            bhr_J₀⁺[i,:] .= collect(sum(J₀⁺[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)' .+ j₀⁺[iμ₀Nstart,1,:] .* qp[iμ₀Nstart]) 
            #if i==1
            #    @show bhr_J₀⁻[i,:], bhr_J₀⁺[i,:]
            #end
            #TODO: Use Radau quadrature and include insolation in the quadrature sum
            
            #@show j₀⁺[iμ₀,1,1:3].* qp[iμ₀], J₀⁺[iμ₀,1,1:3] .* qp[iμ₀], bhr_J₀⁺[i,1:3], bhr_J₀⁻[i,1:3]
            #@show collect(sum(J₀⁻[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)')[1], collect(sum(hdr_J₀⁻[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)')[1],collect(sum(j₀⁻[j,1,:].*wt_μN[j].*qp_μN[j], dims=1)')[1]
            #end
        end
    end
end

function interaction_hdrf_canopy!(SFI,
    #composite_layer::CompositeLayer{FT}, 
    dwJ, uwJ, solJ₀,
    m, pol_type, quad_points,
    hdr_J₀⁻, bhr_J₀⁻, bhr_J₀⁺)

    #@unpack topJ₀⁺, botJ₀⁻      = composite_layer #these are aliases to the respective struct elements 
    @unpack Nquad, wt_μN, iμ₀, qp_μN = quad_points
    NquadN =  Nquad * pol_type.n
    
    wt = collect(wt_μN)
    qp = collect(qp_μN)
    _uwJ = collect(uwJ)
    _dwJ = collect(dwJ)
    #@show typeof(hdr_J₀⁻), typeof(uwJ)
    hdr_J₀⁻ .= collect(uwJ)#r⁻⁺ ⊠ J₀⁺ .+ j₀⁻
    # @show hdr_J₀⁻./ J₀⁺
    
    #@show size(bhr_J₀⁻), size(bhr_J₀⁺)
    #@show size(uwJ), size(dwJ)
    #qp = collect(qp_μN)
    if m==0
        
        #@show dwJ[:,1,end]
        #@show "========="
        #@show uwJ[:,1,end]
        #showp
        for i = 1:pol_type.n
            #bhr_J₀⁻[i,:] .= 0
            #bhr_J₀⁺[i,:] .= 0
            #for j=i:pol_type.n:NquadN
            j=i:pol_type.n:NquadN
            #@show typeof(uwJ[j,1,:]),  typeof(wt_μN[j])
            #collect(sum(collect(uwJ[j,1,:]).*wt[j].*qp[j], dims=1)')
            bhr_J₀⁻[i,:] .= (sum(_uwJ[j,1,:].*wt[j].*qp[j], dims=1)')
            #@show size(bhr_J₀⁺[i,:]), size(solJ₀[i,:])
            bhr_J₀⁺[i,:] .= (sum(_dwJ[j,1,:].*wt[j].*qp[j], dims=1)') .+ (solJ₀[i,:] .* qp[iμ₀]) #TODO: Use Radau quadrature and include insolation in the quadrature sum
            direct = collect(solJ₀[i,:] .* qp[iμ₀])
            diffuse = (sum(_dwJ[j,1,:].*wt[j].*qp[j], dims=1)')
        end
    end

end

"Compute interaction between composite and added layers"