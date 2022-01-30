#=
 
This file contains RT doubling-related functions
 
=#

"""
    $(FUNCTIONNAME)(pol_type, SFI, expk, ndoubl::Int, added_layer::AddedLayer, I_static::AbstractArray{FT}, 
                    architecture) where {FT}

Compute homogenous layer matrices from its elemental layer using Doubling 
"""
function doubling_helper!(RS_type::RRS,
    pol_type, 
    SFI, 
    expk, 
    ndoubl::Int, 
    added_layer::Union{AddedLayer,AddedLayerRS},
    I_static::AbstractArray{FT}, 
    architecture) where {FT}

    # Unpack the added layer
    @unpack i_λ₁λ₀ = RS_type 
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @unpack ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻ = added_layer
    # Device architecture
    dev = devi(architecture)

    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing
    nQuad, _, nSpec = size(r⁺⁻)
    nRaman = length(i_λ₁λ₀);
    # Geometric progression of reflections (1-RR)⁻¹
    gp_refl      = similar(t⁺⁺)
    tt⁺⁺_gp_refl = similar(t⁺⁺)

    if SFI
        # Dummy for source 
        J₁⁺ = similar(J₀⁺)
        # Dummy for J
        J₁⁻ = similar(J₀⁻)

        # Dummy for source 
        ieJ₁⁺ = similar(ieJ₀⁺)
        # Dummy for J
        ieJ₁⁻ = similar(ieJ₀⁻)
    end

    # Loop over number of doublings
    for n = 1:ndoubl

        # T⁺⁺(λ)[I - R⁺⁻(λ)R⁻⁺(λ)]⁻¹, for doubling R⁺⁻,R⁻⁺ and T⁺⁺,T⁻⁻ is identical
        batch_inv!(gp_refl, I_static .- r⁻⁺ ⊠ r⁻⁺)
        tt⁺⁺_gp_refl[:] = t⁺⁺ ⊠ gp_refl
        if SFI
            # J⁺₂₁(λ) = J⁺₁₀(λ).exp(-τ(λ)/μ₀)
            J₁⁺[:,1,:] = J₀⁺[:,1,:] .* expk'
            ieJ₁⁺[:,1,:,:] = ieJ₀⁺[:,1,:,:] .* expk'

            # J⁻₁₂(λ)  = J⁻₀₁(λ).exp(-τ(λ)/μ₀)
            J₁⁻[:,1,:] = J₀⁻[:,1,:] .* expk'
            ieJ₁⁻[:,1,:,:] = ieJ₀⁻[:,1,:,:] .* expk'
            #@show size(ieJ₁⁺)
            tmp1 = gp_refl ⊠  (J₀⁺ + r⁻⁺ ⊠ J₁⁻)
            tmp2 = gp_refl ⊠  (J₁⁻ + r⁻⁺ ⊠ J₀⁺)
            for Δn = 1:nRaman
                n₀, n₁ = get_n₀_n₁(ieJ₁⁺,i_λ₁λ₀[Δn])
                #@show length(n₁), length(n₀), length(n₁_), length(n₀_)
                @inbounds @views ieJ₀⁺[:,:,n₁,Δn] = ieJ₁⁺[:,:,n₁,Δn] + 
                                (tt⁺⁺_gp_refl[:,:,n₁] ⊠ 
                                (ieJ₀⁺[:,:,n₁,Δn] + r⁻⁺[:,:,n₁] ⊠ ieJ₁⁻[:,:,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] ⊠ J₁⁻[:,:,n₀] + 
                                (r⁻⁺[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀]) ⊠ 
                                tmp1[:,:,n₀])) + 
                                iet⁺⁺[:,:,n₁,Δn] ⊠ tmp1[:,:,n₀]
                @inbounds @views ieJ₀⁻[:,:,n₁,Δn] = ieJ₀⁻[:,:,n₁,Δn] + 
                                (tt⁺⁺_gp_refl[:,:,n₁] ⊠ 
                                (ieJ₁⁻[:,:,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] ⊠ J₀⁺[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ieJ₀⁺[:,:,n₁,Δn] + 
                                (ier⁻⁺[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) ⊠ 
                                tmp2[:,:,n₀])) +
                                iet⁻⁻[:,:,n₁,Δn] ⊠ tmp2[:,:,n₀]
            end
            #=
            for n₁ = 1:size(ieJ₁⁺,3) # in eachindex ieJ₁⁺[1,1,:,1]
                for Δn = 1:size(ieJ₁⁺,4) #in eachindex ieJ₁⁺[1,1,1,:]

                    n₀  = n₁ + i_λ₁λ₀[Δn]
                    if 1 ≤ n₀ ≤ size(ieJ₁⁺,3)
                        # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.16 in Raman paper draft)
                        @inbounds ieJ₀⁺[:,1,n₁,Δn] = ieJ₁⁺[:,1,n₁,Δn] + 
                                (tt⁺⁺_gp_refl[:,:,n₁] * 
                                (ieJ₀⁺[:,1,n₁,Δn] + r⁻⁺[:,:,n₁] * ieJ₁⁻[:,1,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] * J₁⁻[:,1,n₀] + 
                                (r⁻⁺[:,:,n₁] * ier⁻⁺[:,:,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] * r⁻⁺[:,:,n₀]) * 
                                gp_refl[:,:,n₀] * (J₀⁺[:,1,n₀] + r⁻⁺[:,:,n₀] * J₁⁻[:,1,n₀]))) + 
                                iet⁺⁺[:,:,n₁,Δn] * gp_refl[:,:,n₀] * 
                                (J₀⁺[:,1,n₀] + r⁻⁺[:,:,n₀] * J₁⁻[:,1,n₀])

                        # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.17 in Raman paper draft)
                        @inbounds ieJ₀⁻[:,1,n₁,Δn] = ieJ₀⁻[:,1,n₁,Δn] + 
                                (tt⁺⁺_gp_refl[:,:,n₁] * 
                                (ieJ₁⁻[:,1,n₁,Δn] + ier⁻⁺[:,:,n₁,Δn] * J₀⁺[:,1,n₀] + r⁻⁺[:,:,n₁] * ieJ₀⁺[:,1,n₁,Δn] + 
                                (ier⁻⁺[:,:,n₁,Δn] * r⁻⁺[:,:,n₀] + r⁻⁺[:,:,n₁] * ier⁻⁺[:,:,n₁,Δn]) * 
                                gp_refl[:,:,n₀] * (J₁⁻[:,1,n₀] + r⁻⁺[:,:,n₀] * J₀⁺[:,1,n₀]))) +
                                iet⁻⁻[:,:,n₁,Δn] * gp_refl[:,:,n₀] * (J₁⁻[:,1,n₀] + 
                                r⁻⁺[:,:,n₀]*J₀⁺[:,1,n₀])
                    end
                end
            end 
            =#
        
            # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.8 in Raman paper draft)
            J₀⁻[:] = J₀⁻ + (tt⁺⁺_gp_refl ⊠ (J₁⁻ + r⁻⁺ ⊠ J₀⁺)) 

            # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.8 in Raman paper draft)
            J₀⁺[:] = J₁⁺ + (tt⁺⁺_gp_refl ⊠ (J₀⁺ + r⁻⁺ ⊠ J₁⁻))

            expk[:] = expk.^2
        end  
        #println("Doubling part 1 done")
        for Δn = 1:nRaman
                n₀, n₁ = get_n₀_n₁(ieJ₁⁺,i_λ₁λ₀[Δn])
                #@show n₁, n₀
                #@show length(n₀)
                @inbounds @views iet⁺⁺[:,:,n₁,Δn] = t⁺⁺[:,:,n₁] ⊠ gp_refl[:,:,n₁] ⊠ 
                        (iet⁺⁺[:,:,n₁,Δn] + 
                        (ier⁻⁺[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) ⊠ 
                        gp_refl[:,:,n₀] ⊠ t⁺⁺[:,:,n₀]) + 
                        iet⁺⁺[:,:,n₁,Δn] ⊠ gp_refl[:,:,n₀] ⊠  t⁺⁺[:,:,n₀]
                @inbounds @views ier⁻⁺[:,:,n₁,Δn] = ier⁻⁺[:,:,n₁,Δn] + 
                        t⁺⁺[:,:,n₁] ⊠ gp_refl[:,:,n₁] ⊠ r⁻⁺[:,:,n₁] ⊠  
                        (iet⁺⁺[:,:,n₁,Δn] + 
                        (ier⁻⁺[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + r⁻⁺[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) ⊠ 
                        gp_refl[:,:,n₀] ⊠ t⁺⁺[:,:,n₀]) + 
                        (iet⁺⁺[:,:,n₁,Δn] ⊠ gp_refl[:,:,n₀] ⊠ r⁻⁺[:,:,n₀] + 
                        t⁺⁺[:,:,n₁] ⊠ gp_refl[:,:,n₁] ⊠ ier⁻⁺[:,:,n₁,Δn]) ⊠ t⁺⁺[:,:,n₀]
        end
        #=for n₁ = 1:size(ieJ₁⁺,3) #in eachindex ieJ₁⁺[1,1,:,1]
            for Δn = 1:size(ieJ₁⁺,4) #in eachindex ieJ₁⁺[1,1,1,:]
                n₀  = n₁ + i_λ₁λ₀[Δn]
                # (see Eqs.12 in Raman paper draft)
                if 1 ≤ n₀ ≤ size(ieJ₁⁺,3)
                    iet⁺⁺[:,:,n₁,Δn] = t⁺⁺[:,:,n₁] * gp_refl[:,:,n₁] * 
                        (iet⁺⁺[:,:,n₁,Δn] + 
                        (ier⁻⁺[:,:,n₁,Δn] * r⁻⁺[:,:,n₀] + r⁻⁺[:,:,n₁] * ier⁻⁺[:,:,n₁,Δn]) * 
                        gp_refl[:,:,n₀] * t⁺⁺[:,:,n₀]) + 
                        iet⁺⁺[:,:,n₁,Δn] * gp_refl[:,:,n₀] * t⁺⁺[:,:,n₀]

                    # (see Eqs.14 in Raman paper draft)
                    ier⁻⁺[:,:,n₁,Δn] = ier⁻⁺[:,:,n₁,Δn] + 
                        t⁺⁺[:,:,n₁] * gp_refl[:,:,n₁] * r⁻⁺[:,:,n₁] * 
                        (iet⁺⁺[:,:,n₁,Δn] + 
                        (ier⁻⁺[:,:,n₁,Δn] * r⁻⁺[:,:,n₀] + r⁻⁺[:,:,n₁] * ier⁻⁺[:,:,n₁,Δn]) * 
                        gp_refl[:,:,n₀] * t⁺⁺[:,:,n₀]) + 
                        (iet⁺⁺[:,:,n₁,Δn] * gp_refl[:,:,n₀] * r⁻⁺[:,:,n₀] + 
                        t⁺⁺[:,:,n₁] * gp_refl[:,:,n₁] * ier⁻⁺[:,:,n₁,Δn]) * t⁺⁺[:,:,n₀]
                end
            end
        end
        =#
        # R⁻⁺₂₀(λ) = R⁻⁺₁₀(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹R⁻⁺₂₁(λ)T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        r⁻⁺[:]  = r⁻⁺ + (tt⁺⁺_gp_refl ⊠ r⁻⁺ ⊠ t⁺⁺)

        # T⁺⁺₂₀(λ) = T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        t⁺⁺[:]  = tt⁺⁺_gp_refl ⊠ t⁺⁺
    end

    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    # For SFI, after doubling, revert D(DJ₀⁻)->J₀⁻

    synchronize_if_gpu()

    apply_D_matrix!(pol_type.n, added_layer.r⁻⁺, added_layer.t⁺⁺, added_layer.r⁺⁻, added_layer.t⁻⁻)
    apply_D_matrix_IE!(RS_type, pol_type.n, added_layer.ier⁻⁺, added_layer.iet⁺⁺, added_layer.ier⁺⁻, added_layer.iet⁻⁻)
    SFI && apply_D_matrix_SFI!(pol_type.n, added_layer.J₀⁻)
    SFI && apply_D_matrix_SFI_IE!(RS_type, pol_type.n, added_layer.ieJ₀⁻)

    return nothing 
end


function doubling_helper!(RS_type::Union{VS_0to1, VS_1to0},
                        pol_type, 
                        SFI, 
                        expk, 
                        ndoubl::Int, 
                        added_layer::AddedLayer,
                        I_static::AbstractArray{FT}, 
                        architecture) where {FT}
    # Unpack the added layer
    @unpack i_λ₁λ₀ = RS_type 
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻ = added_layer
    @unpack ier⁺⁻, ier⁻⁺, iet⁻⁻, iet⁺⁺, ieJ₀⁺, ieJ₀⁻ = added_layer
    # Device architecture
    dev = devi(architecture)

    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing
    
    # Geometric progression of reflections (1-RR)⁻¹
    gp_refl      = similar(t⁺⁺)
    tt⁺⁺_gp_refl = similar(t⁺⁺)
 
    if SFI
        # Dummy for source 
        J₁⁺ = similar(J₀⁺)
        # Dummy for J
        J₁⁻ = similar(J₀⁻)

        # Dummy for source 
        ieJ₁⁺ = similar(ieJ₀⁺)
        # Dummy for J
        ieJ₁⁻ = similar(ieJ₀⁻)
    end

    # Loop over number of doublings
    for n = 1:ndoubl
        
        # T⁺⁺(λ)[I - R⁺⁻(λ)R⁻⁺(λ)]⁻¹, for doubling R⁺⁻,R⁻⁺ and T⁺⁺,T⁻⁻ is identical
        batch_inv!(gp_refl, I_static .- r⁻⁺ ⊠ r⁻⁺)
        tt⁺⁺_gp_refl[:] = t⁺⁺ ⊠ gp_refl

        if SFI

            # J⁺₂₁(λ) = J⁺₁₀(λ).exp(-τ(λ)/μ₀)
            J₁⁺[:,1,:] = J₀⁺[:,1,:] .* expk'
            ieJ₁⁺[:,1,:] = ieJ₀⁺[:,1,:] .* expk'

            # J⁻₁₂(λ)  = J⁻₀₁(λ).exp(-τ(λ)/μ₀)
            J₁⁻[:,1,:] = J₀⁻[:,1,:] .* expk'
            ieJ₁⁻[:,1,:] = ieJ₀⁻[:,1,:] .* expk'

            #for n₁ in eachindex ieJ₁⁺[1,1,:,1]
            for Δn in eachindex ieJ₁⁺[1,1,1,:]
                n₁  = i_λ₁λ₀[Δn]
                # TODO: replace all the n₀ indices with ref_r, ref_t, and ref_J
                # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.16 in Raman paper draft)
                ieJ₀⁺[:,1,n₁,1] = ieJ₁⁺[:,1,n₁,1] + 
                        (tt⁺⁺_gp_refl[:,:,n₁] * (ieJ₀⁺[:,1,n₁,1] + 
                        r⁻⁺[:,:,n₁] * ieJ₁⁻[:,1,n₁,1] + 
                        ier⁻⁺[:,:,n₁,1] * J₁⁻[:,1,n₀] + 
                        (r⁻⁺[:,:,n₁] * ier⁻⁺[:,:,n₁,1] + 
                        ier⁻⁺[:,:,n₁,1] * r⁻⁺[:,:,n₀]) * gp_refl[:,:,n₀] * (J₀⁺[:,1,n₀] + 
                        r⁻⁺[:,:,n₀] * J₁⁻[:,1,n₀]))) + 
                        iet⁺⁺[:,:,n₁,1] * gp_refl[:,:,n₀] * 
                        (J₀⁺[:,1,n₀] + r⁻⁺[:,:,n₀] * J₁⁻[:,1,n₀])
        
                # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.17 in Raman paper draft)
                ieJ₀⁻[:,1,n₁,1] = ieJ₀⁻[:,1,n₁,1] + 
                        (tt⁺⁺_gp_refl[:,:,n₁] * (ieJ₁⁻[:,1,n₁,1] + 
                        ier⁻⁺[:,:,n₁,1] * J₀⁺[:,1,n₀] + 
                        r⁻⁺[:,:,n₁] * ieJ₀⁺[:,1,n₁,1] + 
                        (ier⁻⁺[:,:,n₁,1] * r⁻⁺[:,:,n₀] + 
                        r⁻⁺[:,:,n₁] * ier⁻⁺[:,:,n₁,1]) * gp_refl[:,:,n₀] * (J₁⁻[:,1,n₀] + 
                        r⁻⁺[:,:,n₀] * J₀⁺[:,1,n₀]))) +
                        iet⁻⁻[:,:,n₁,1] * gp_refl[:,:,n₀] * (J₁⁻[:,1,n₀] + 
                        r⁻⁺[:,:,n₀]*J₀⁺[:,1,n₀])
            end 
            #end

            
            # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.8 in Raman paper draft)
            J₀⁻[:] = J₀⁻ + (tt⁺⁺_gp_refl ⊠ (J₁⁻ + r⁻⁺ ⊠ J₀⁺)) 

            # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.8 in Raman paper draft)
            J₀⁺[:] = J₁⁺ + (tt⁺⁺_gp_refl ⊠ (J₀⁺ + r⁻⁺ ⊠ J₁⁻))
             
            expk[:] = expk.^2
        end  

        #for n₁ in eachindex ieJ₁⁺[1,1,:,1]
        for Δn in eachindex ieJ₁⁺[1,1,1,:]
            n₁ = i_λ₁λ₀[Δn]
            # (see Eqs.12 in Raman paper draft)
            iet⁺⁺[:,:,n₁,1] = t⁺⁺[:,:,n₁] * gp_refl[:,:,n₁] * 
                    (iet⁺⁺[:,:,n₁,1] + 
                    (ier⁻⁺[:,:,n₁,1] * r⁻⁺[:,:,n₀] + 
                    r⁻⁺[:,:,n₁] * ier⁻⁺[:,:,n₁,1]) * 
                    gp_refl[:,:,n₀] * t⁺⁺[:,:,n₀]) + 
                    iet⁺⁺[:,:,n₁,1] * gp_refl[:,:,n₀] * t⁺⁺[:,:,n₀]

            # (see Eqs.14 in Raman paper draft)
            ier⁻⁺[:,:,n₁,1] = ier⁻⁺[:,:,n₁,1] + 
                    t⁺⁺[:,:,n₁] * gp_refl[:,:,n₁] * r⁻⁺[:,:,n₁] * 
                    (iet⁺⁺[:,:,n₁,1] + 
                    (ier⁻⁺[:,:,n₁,1] * r⁻⁺[:,:,n₀] + 
                    r⁻⁺[:,:,n₁] * ier⁻⁺[:,:,n₁,1]) * 
                    gp_refl[:,:,n₀] * t⁺⁺[:,:,n₀]) + 
                    (iet⁺⁺[:,:,n₁,1] * gp_refl[:,:,n₀] * r⁻⁺[:,:,n₀] + 
                    t⁺⁺[:,:,n₁] * gp_refl[:,:,n₁] * ier⁻⁺[:,:,n₁,1]) * t⁺⁺[:,:,n₀]
        end
    
        # R⁻⁺₂₀(λ) = R⁻⁺₁₀(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹R⁻⁺₂₁(λ)T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        r⁻⁺[:]  = r⁻⁺ + (tt⁺⁺_gp_refl ⊠ r⁻⁺ ⊠ t⁺⁺)

        # T⁺⁺₂₀(λ) = T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        t⁺⁺[:]  = tt⁺⁺_gp_refl ⊠ t⁺⁺
    end

    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    # For SFI, after doubling, revert D(DJ₀⁻)->J₀⁻

    synchronize_if_gpu()

    apply_D_matrix!(pol_type.n, 
        added_layer.r⁻⁺, added_layer.t⁺⁺, added_layer.r⁺⁻, added_layer.t⁻⁻)
    apply_D_matrix_IE!(RS_type, pol_type.n, 
        added_layer.ier⁻⁺, added_layer.iet⁺⁺, added_layer.ier⁺⁻, added_layer.iet⁻⁻)
    SFI && apply_D_matrix_SFI!(pol_type.n, added_layer.J₀⁻)
    SFI && apply_D_matrix_SFI_IE!(RS_type, pol_type.n, added_layer.ieJ₀⁻)
    
    return nothing 
end

function doubling_inelastic!(RS_type, 
                    pol_type, SFI, 
                    expk, ndoubl::Int, 
                    added_layer::Union{AddedLayer,AddedLayerRS},#{FT},
                    I_static::AbstractArray{FT}, 
                    architecture) where {FT}

    doubling_helper!(RS_type, 
                pol_type, SFI, 
                expk, ndoubl, 
                added_layer, 
                I_static, 
                architecture)

    synchronize_if_gpu()
end

#@kernel function apply_D!(n_stokes::Int,  r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
#    iμ, jμ, n = @index(Global, NTuple)
#    i = mod(iμ, n_stokes)
#    j = mod(jμ, n_stokes)
#
#    if (i > 2)
#        r⁻⁺[iμ,jμ,n] = - r⁻⁺[iμ, jμ, n]
#        r⁻⁺[iμ,jμ,n] = - r⁻⁺[iμ, jμ, n]
#    end
#    
#    if ((i <= 2) & (j <= 2)) | ((i > 2) & (j > 2))
#        r⁺⁻[iμ,jμ,n] = r⁻⁺[iμ,jμ,n]
#        t⁻⁻[iμ,jμ,n] = t⁺⁺[iμ,jμ,n]
#    else
#        r⁺⁻[iμ,jμ,n] = - r⁻⁺[iμ,jμ,n]
#        t⁻⁻[iμ,jμ,n] = - t⁺⁺[iμ,jμ,n]
#    end
#end

@kernel function apply_D_IE!(RS_type::RRS, n_stokes::Int,  
                        ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻)
    iμ, jμ, n, Δn  = @index(Global, NTuple)
    @unpack i_λ₁λ₀ = RS_type 
    n₀  = n + i_λ₁λ₀[Δn]
    i = mod(iμ, n_stokes)
    j = mod(jμ, n_stokes)

    if (i > 2)
        ier⁻⁺[iμ,jμ,n,n₀] = - ier⁻⁺[iμ, jμ, n, n₀]
    end
    
    if ((i <= 2) & (j <= 2)) | ((i > 2) & (j > 2))
        ier⁺⁻[iμ,jμ,n,n₀] = ier⁻⁺[iμ,jμ,n,n₀]
        iet⁻⁻[iμ,jμ,n,n₀] = iet⁺⁺[iμ,jμ,n,n₀]
    else
        ier⁺⁻[iμ,jμ,n,n₀] = - ier⁻⁺[iμ,jμ,n,n₀]
        iet⁻⁻[iμ,jμ,n,n₀] = - iet⁺⁺[iμ,jμ,n,n₀]
    end

end

@kernel function apply_D_IE!(RS_type::Union{VS_0to1, VS_1to0}, n_stokes::Int,  
                        ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻)
    iμ, jμ, Δn  = @index(Global, NTuple)
    @unpack i_λ₁λ₀ = RS_type 
    n  = i_λ₁λ₀[Δn]
    i = mod(iμ, n_stokes)
    j = mod(jμ, n_stokes)

    if (i > 2)
        ier⁻⁺[iμ,jμ,n,1] = - ier⁻⁺[iμ, jμ, n, 1]
    end
    
    if ((i <= 2) & (j <= 2)) | ((i > 2) & (j > 2))
        ier⁺⁻[iμ,jμ,n,1] = ier⁻⁺[iμ,jμ,n,1]
        iet⁻⁻[iμ,jμ,n,1] = iet⁺⁺[iμ,jμ,n,1]
    else
        ier⁺⁻[iμ,jμ,n,1] = - ier⁻⁺[iμ,jμ,n,1]
        iet⁻⁻[iμ,jμ,n,1] = - iet⁺⁺[iμ,jμ,n,1]
    end

end

#@kernel function apply_D_SFI!(n_stokes::Int, J₀⁻)
#    iμ, _, n = @index(Global, NTuple)
#    i = mod(iμ, n_stokes)
#
#    if (i > 2)
#        J₀⁻[iμ, 1, n] = - J₀⁻[iμ, 1, n] 
#    end
#end

@kernel function apply_D_SFI_IE!(RS_type::RRS, n_stokes::Int, ieJ₀⁻)
    iμ, n, Δn = @index(Global, NTuple)
    @unpack i_λ₁λ₀ = RS_type
    n₀ = n + i_λ₁λ₀[Δn] 
    i = mod(iμ, n_stokes)

    if (i > 2)
        ieJ₀⁻[iμ, 1, n, n₀] = - ieJ₀⁻[iμ, 1, n, Δn] 
    end
end

@kernel function apply_D_SFI_IE!(RS_type::Union{VS_0to1, VS_1to0}, 
                                n_stokes::Int, ieJ₀⁻)
    iμ, _, Δn = @index(Global, NTuple)
    @unpack i_λ₁λ₀ = RS_type
    n = i_λ₁λ₀[Δn] 
    i = mod(iμ, n_stokes)

    if (i > 2)
        ieJ₀⁻[iμ, 1, n, 1] = - ieJ₀⁻[iμ, 1, n, 1] 
    end
end

#Suniti: is it possible to  use the same kernel for the 3D elastic and 4D inelastic terms or do we need to call two different kernels separately? 
#function apply_D_matrix!(n_stokes::Int, r⁻⁺::CuArray{FT,3}, t⁺⁺::CuArray{FT,3}, r⁺⁻::CuArray{FT,3}, t⁻⁻::CuArray{FT,3}) where {FT}
#    
#    if n_stokes == 1
#        r⁺⁻[:] = r⁻⁺
#        t⁻⁻[:] = t⁺⁺    
#        
#        return nothing
#    else 
#        device = devi(architecture(r⁻⁺))
#        applyD_kernel! = apply_D!(device)
#        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺)); #Suniti: is it possible to  use the same kernel for the 3D elastic and 4D inelastic terms or do we need to call two different kernels separately? 
#        wait(device, event);
#        synchronize_if_gpu();
#        return nothing
#    end
#end

function apply_D_matrix_IE!(RS_type, n_stokes::Int, ier⁻⁺::Array{FT,4}, iet⁺⁺::Array{FT,4}, ier⁺⁻::Array{FT,4}, iet⁻⁻::Array{FT,4}) where {FT}
    if n_stokes == 1
        ier⁺⁻[:] = ier⁻⁺
        iet⁻⁻[:] = iet⁺⁺  
        return nothing
    else 
        device = devi(Architectures.CPU())
        applyD_kernel_IE! = apply_D_IE!(device)
        event = applyD_kernel_IE!(RS_type, n_stokes, 
            ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻, ndrange=getKernelDim(RS_type, ier⁻⁺));
        wait(device, event);
        return nothing
    end
end

function apply_D_matrix_IE!(RS_type, n_stokes::Int, ier⁻⁺::CuArray{FT,4}, iet⁺⁺::CuArray{FT,4}, ier⁺⁻::CuArray{FT,4}, iet⁻⁻::CuArray{FT,4}) where {FT}
    if n_stokes == 1
        ier⁺⁻[:] = ier⁻⁺
        iet⁻⁻[:] = iet⁺⁺  
        return nothing
    else 
        device = devi(Architectures.GPU())
        applyD_kernel_IE! = apply_D_IE!(device)
        event = applyD_kernel_IE!(RS_type, n_stokes, 
            ier⁻⁺, iet⁺⁺, ier⁺⁻, iet⁻⁻, ndrange=getKernelDim(RS_type, ier⁻⁺));
        wait(device, event);
        synchronize();
        return nothing
    end
end

#function apply_D_matrix_SFI!(n_stokes::Int, J₀⁻::CuArray{FT,3}) where {FT}
#
#    n_stokes == 1 && return nothing
#    device = devi(architecture(J₀⁻)) #Suniti: how to do this so that ieJ₀⁻ can also be included?
#    applyD_kernel! = apply_D_SFI!(device)
#    event = applyD_kernel!(n_stokes, J₀⁻, ndrange=size(J₀⁻));
#    wait(device, event);
#    synchronize();
#    
#    return nothing
#end
    
function apply_D_matrix_SFI_IE!(RS_type, n_stokes::Int, ieJ₀⁻::CuArray{FT,4}) where {FT}
    
    n_stokes == 1 && return nothing

    device = devi(Architectures.GPU())
    applyD_kernel_IE! = apply_D_SFI_IE!(device)
    event = applyD_kernel_IE!(RS_type, n_stokes, 
                    ieJ₀⁻, ndrange=getKernelDimSFI(RS_type, ieJ₀⁻));
    wait(device, event);
    synchronize()
    return nothing
end

function apply_D_matrix_SFI_IE!(RS_type, n_stokes::Int, ieJ₀⁻::Array{FT,4}) where {FT}
    
    n_stokes == 1 && return nothing

    device = devi(Architectures.CPU())
    applyD_kernel_IE! = apply_D_SFI_IE!(device)
    event = applyD_kernel_IE!(RS_type, n_stokes, 
                    ieJ₀⁻, ndrange=getKernelDimSFI(RS_type, ieJ₀⁻));
    wait(device, event);
    return nothing
end

