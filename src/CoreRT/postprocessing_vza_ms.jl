#=

This file contains the function to perform azimuthal-weighting to the RT matrices after all 
kernel calculations. 

=#

"Perform post-processing to azimuthally-weight RT matrices"
function postprocessing_vza_ms!(RS_type::noRS, sensor_levels, 
        iμ₀, pol_type, 
        composite_layer, 
        vza, qp_μ, m, vaz, μ₀, weight, 
        nSpec, SFI, uwJ, dwJ, #ie_uwJ, ie_dwJ, 
        I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    
    # idfx of μ0 = cos(sza)
    st_iμ0, istart0, iend0 = get_indices(iμ₀, pol_type);
    
    # Convert these to Arrays (if CuArrays), so they can be accessed by index
    topR⁺⁻ = Array(composite_layer.topR⁺⁻);
    botR⁻⁺ = Array(composite_layer.botR⁻⁺);
    topJ₀⁺ = Array(composite_layer.topJ₀⁺);
    topJ₀⁻ = Array(composite_layer.topJ₀⁻);
    botJ₀⁻ = Array(composite_layer.botJ₀⁻);
    
    tuwJ = similar(topJ₀⁺)
    tdwJ = similar(topJ₀⁺)
    for ims = 1:length(sensor_levels)
        if(sensor_levels[ims]==0)
            tuwJ[ims,:] = topJ₀⁻[ims,:]
            tdwJ[ims,:] = topJ₀⁺[ims,:]
        else
            compute_interlayer_flux!(RS_type,
                                    I_static, 
                                    topR⁺⁻[ims,:], botR⁻⁺[ims,:], 
                                    topJ₀⁺[ims,:], botJ₀⁻[ims,:],
                                    tdwJ[ims,:], tuwJ[ims,:]);

            #tmpR = @timeit "postprocessing_vza inv1" batch_inv!(
            #        tmp_inv, I_static .- topR⁺⁻[ims,:] ⊠ botR⁻⁺[ims,:]) #Suniti
            #tdwJ[ims,:] = tmpR ⊠ (topJ₀⁺[ims,:] .+ topR⁺⁻[ims,:]⊠botJ₀⁻[ims,:])

            #tmpR = @timeit "postprocessing_vza inv1" batch_inv!(
            #        tmp_inv, I_static .- botR⁻⁺[ims,:] ⊠ topR⁺⁻[ims,:]) #Suniti
            #tuwJ[ims,:] = tmpR ⊠ (botJ₀⁻[ims,:] .+ botR⁻⁺[ims,:]⊠topJ₀⁺[ims,:])
        end
    end
    # Loop over all viewing zenith angles
    for i = 1:length(vza)

        # Find the nearest quadrature point idx
        iμ = nearest_point(qp_μ, cosd(vza[i]));
        st_iμ, istart, iend = get_indices(iμ, pol_type);
        
        # Compute bigCS
        cos_m_phi, sin_m_phi = (cosd(m * vaz[i]), sind(m * vaz[i]));
        bigCS = weight * Diagonal([cos_m_phi, cos_m_phi, sin_m_phi, sin_m_phi][1:pol_type.n]);

        # Accumulate Fourier moments after azimuthal weighting
        for ims = 1:length(sensor_levels)
            for s = 1:nSpec
                if true#SFI
                    uwJ[ims,i,:,s] += bigCS * tuwJ[ims,istart:iend,1, s];
                    dwJ[ims,i,:,s] += bigCS * tdwJ[ims,istart:iend,1, s];
                # else
                #     R[i,:,s] += bigCS * (R⁻⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                #     T[i,:,s] += bigCS * (T⁺⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                end
            end
        end
    end
end

function postprocessing_vza_ms!(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus}, 
        iμ₀, pol_type, 
        composite_layer, 
        vza, qp_μ, m, vaz, μ₀, weight, 
        nSpec, SFI, 
        uwJ, dwJ, uwieJ, dwieJ,
        I_static::AbstractArray{FT2}) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
        #R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)
    
    # idx of μ0 = cos(sza)
    st_iμ0, istart0, iend0 = get_indices(iμ₀, pol_type);

    # Convert these to Arrays (if CuArrays), so they can be accessed by index
    topR⁺⁻ = Array(composite_layer.topR⁺⁻);
    botR⁻⁺ = Array(composite_layer.botR⁻⁺);
    topJ₀⁺ = Array(composite_layer.topJ₀⁺);
    topJ₀⁻ = Array(composite_layer.topJ₀⁻);
    botJ₀⁻ = Array(composite_layer.botJ₀⁻);
    
    tuwJ = similar(topJ₀⁺)
    tdwJ = similar(topJ₀⁺)

    botieR⁻⁺ = Array(composite_layer.botieR⁻⁺);
    topieR⁺⁻ = Array(composite_layer.topieR⁺⁻);
    
    topieJ₀⁺ = Array(composite_layer.topieJ₀⁺);
    topieJ₀⁻ = Array(composite_layer.topieJ₀⁻);
    botieJ₀⁻ = Array(composite_layer.botieJ₀⁻);

    tuwieJ = similar(topieJ₀⁺)
    tdwieJ = similar(topieJ₀⁺)

    for ims = 1:length(sensor_levels)
        if(sensor_levels[ims]==0)
            # elastic
            tuwJ[ims,:] = topJ₀⁻[ims,:]
            tdwJ[ims,:] = topJ₀⁺[ims,:]
            # inelastic
            tuwieJ[ims,:] = topieJ₀⁻[ims,:]
            tdwieJ[ims,:] = topieJ₀⁺[ims,:]
        else
            compute_interlayer_flux!(RS_type,
                                    I_static, 
                                    topR⁺⁻[ims,:], botR⁻⁺[ims,:], 
                                    topJ₀⁺[ims,:], botJ₀⁻[ims,:],
                                    tdwJ[ims,:], tuwJ[ims,:],
                                    topieR⁺⁻[ims,:], botieR⁻⁺[ims,:], 
                                    topieJ₀⁺[ims,:], botieJ₀⁻[ims,:],
                                    tdwieJ[ims,:], tuwieJ[ims,:]);
        end
    end
    # Loop over all viewing zenith angles
    for i = 1:length(vza)

        # Find the nearest quadrature point idx
        iμ = nearest_point(qp_μ, cosd(vza[i]));
        st_iμ, istart, iend = get_indices(iμ, pol_type);
        
        # Compute bigCS
        cos_m_phi, sin_m_phi = (cosd(m * vaz[i]), sind(m * vaz[i]));
        bigCS = weight * Diagonal([cos_m_phi, cos_m_phi, sin_m_phi, sin_m_phi][1:pol_type.n]);

        # Accumulate Fourier moments after azimuthal weighting
        for ims = 1:length(sensor_levels)
            for s = 1:nSpec
                
                if true#SFI
                    uwJ[ims,i,:,s] += bigCS * tuwJ[ims,istart:iend,1, s];
                    dwJ[ims,i,:,s] += bigCS * tdwJ[ims,istart:iend,1, s];
                    #R_SFI[i,:,s] += bigCS * J₀⁻[istart:iend,1, s];
                    #T_SFI[i,:,s] += bigCS * J₀⁺[istart:iend,1, s];
                    #@show i, s, R_SFI[i,:,s]
                    #@show i, s, ieR_SFI[i,:,s]
                    
                    uwieJ[ims,i,:,s] += bigCS * tuwieJ[ims,istart:iend,1, s];
                    dwieJ[ims,i,:,s] += bigCS * tdwieJ[ims,istart:iend,1, s];
                        #@show i, s, t, ieR_SFI[i,:,s]
                    
                    #@show i, s, ieR_SFI[i,:,s]
                #else
                #    R[i,:,s] += bigCS * (R⁻⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                #    T[i,:,s] += bigCS * (T⁺⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                    #no modification for Raman because SFI will be the only mode used for inelastic computations
                end
                
            end
        end
    end
end