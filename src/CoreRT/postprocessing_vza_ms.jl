#=

This file contains the function to perform azimuthal-weighting to the RT matrices after all 
kernel calculations. 

=#

"Perform post-processing to azimuthally-weight RT matrices"
function postprocessing_vza_ms!(RS_type::noRS, 
        sensor_levels, 
        iμ₀, pol_type, 
        composite_layer, 
        vza, qp_μ, m, vaz, μ₀, weight, 
        nSpec, SFI, 
        uwJ, dwJ, uwieJ, dwieJ, 
        I_static::AbstractArray{FT2}, 
        arr_type) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    
    # idfx of μ0 = cos(sza)
    st_iμ0, istart0, iend0 = get_indices(iμ₀, pol_type);
    
    # Convert these to Arrays (if CuArrays), so they can be accessed by index
    topR⁺⁻ = composite_layer.topR⁺⁻;
    botR⁻⁺ = composite_layer.botR⁻⁺;
    topJ₀⁺ = composite_layer.topJ₀⁺;
    topJ₀⁻ = composite_layer.topJ₀⁻;
    botJ₀⁺ = composite_layer.botJ₀⁺;
    botJ₀⁻ = composite_layer.botJ₀⁻;
    
    tuwJ = [zeros(typeof(topJ₀⁺[1][1,1,1]), (length(topJ₀⁺[1][:,1,1]), 1, nSpec)) for i=1:length(sensor_levels)] #similar(topJ₀⁺); #deepcopy(topJ₀⁺)
    tdwJ = [zeros(typeof(topJ₀⁺[1][1,1,1]), (length(topJ₀⁺[1][:,1,1]), 1, nSpec)) for i=1:length(sensor_levels)]#similar(topJ₀⁺); #deepcopy(topJ₀⁺)
    #@show size(tuwJ[1])
    for ims = 1:length(sensor_levels)
        if(sensor_levels[ims]==0)
            tuwJ[ims] .= botJ₀⁻[ims]
            tdwJ[ims] .= botJ₀⁺[ims]
        else
            
            #@show size(topR⁺⁻[ims]), size(botR⁻⁺[ims]),size(topJ₀⁺[ims]),size(botJ₀⁻[ims])#,size(tdwJ[ims]),size(tuwJ[ims])
            compute_interlayer_flux!(RS_type,
                                    I_static, 
                                    topR⁺⁻[ims], botR⁻⁺[ims], 
                                    topJ₀⁺[ims], botJ₀⁻[ims],
                                    tdwJ[ims], tuwJ[ims], arr_type);

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

        #@show tuwJ[1][istart:iend,1, 1], tdwJ[1][istart:iend,1, 1]
        # Accumulate Fourier moments after azimuthal weighting
        for ims = 1:length(sensor_levels)
            for s = 1:nSpec
                #if true#SFI
                uwJ[ims][i,:,s] .+= bigCS * tuwJ[ims][istart:iend,1, s];
                dwJ[ims][i,:,s] .+= bigCS * tdwJ[ims][istart:iend,1, s];
                # else
                #     R[i,:,s] += bigCS * (R⁻⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                #     T[i,:,s] += bigCS * (T⁺⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                #end
            end
        end
    end
end

function postprocessing_vza_ms!(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus}, 
        sensor_levels,
        iμ₀, pol_type, 
        composite_layer, 
        vza, qp_μ, m, vaz, μ₀, weight, 
        nSpec, SFI, 
        uwJ, dwJ, uwieJ, dwieJ,
        I_static,
        arr_type) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
        #R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)
    
    # idx of μ0 = cos(sza)
    st_iμ0, istart0, iend0 = get_indices(iμ₀, pol_type);

    # Convert these to Arrays (if CuArrays), so they can be accessed by index
    topR⁺⁻ = (composite_layer.topR⁺⁻);
    botR⁻⁺ = (composite_layer.botR⁻⁺);
    topJ₀⁺ = (composite_layer.topJ₀⁺);
    topJ₀⁻ = (composite_layer.topJ₀⁻);
    botJ₀⁺ = (composite_layer.botJ₀⁺);
    botJ₀⁻ = (composite_layer.botJ₀⁻);
    
    tuwJ = [zeros(typeof(topJ₀⁺[1][1,1,1]), (length(topJ₀⁺[1][:,1,1]), 1, nSpec)) for i=1:length(sensor_levels)] #similar(topJ₀⁺); #deepcopy(topJ₀⁺)
    tdwJ = [zeros(typeof(topJ₀⁺[1][1,1,1]), (length(topJ₀⁺[1][:,1,1]), 1, nSpec)) for i=1:length(sensor_levels)]#similar(topJ₀⁺); #deepcopy(topJ₀⁺)
    #@show size(tuwJ), size(tdwJ)   
    botieR⁻⁺ = (composite_layer.botieR⁻⁺);
    topieR⁺⁻ = (composite_layer.topieR⁺⁻);
    
    topieJ₀⁺ = (composite_layer.topieJ₀⁺);
    #topieJ₀⁻ = (composite_layer.topieJ₀⁻);
    botieJ₀⁺ = (composite_layer.botieJ₀⁺);
    botieJ₀⁻ = (composite_layer.botieJ₀⁻);

    #tuwieJ = deepcopy(topieJ₀⁺)
    #tdwieJ = deepcopy(topieJ₀⁺)
    tuwieJ = [zeros(typeof(topJ₀⁺[1][1,1,1]), (length(topJ₀⁺[1][:,1,1]), 1, nSpec, length(topieJ₀⁺[1][1,1,1,:]))) for i=1:length(sensor_levels)] #similar(topJ₀⁺); #deepcopy(topJ₀⁺)
    tdwieJ = [zeros(typeof(topJ₀⁺[1][1,1,1]), (length(topJ₀⁺[1][:,1,1]), 1, nSpec, length(topieJ₀⁺[1][1,1,1,:]))) for i=1:length(sensor_levels)] #similar(topJ₀⁺); #deepcopy(topJ₀⁺)
    #@show size(tuwieJ), size(tdwieJ)  
    for ims = 1:length(sensor_levels)
        if(sensor_levels[ims]==0)
            # elastic
            tuwJ[ims] .= botJ₀⁻[ims]
            tdwJ[ims] .= botJ₀⁺[ims]
            # inelastic
            tuwieJ[ims] .= botieJ₀⁻[ims]
            tdwieJ[ims] .= botieJ₀⁺[ims]
        else
            compute_interlayer_flux!(RS_type,
                            I_static, 
                            topR⁺⁻[ims], botR⁻⁺[ims], 
                            topJ₀⁺[ims], botJ₀⁻[ims],
                            tdwJ[ims], tuwJ[ims],
                            topieR⁺⁻[ims], botieR⁻⁺[ims], 
                            topieJ₀⁺[ims], botieJ₀⁻[ims],
                            tdwieJ[ims], tuwieJ[ims], arr_type);
        
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
                    uwJ[ims][i,:,s] .+= bigCS * tuwJ[ims][istart:iend,1, s];
                    dwJ[ims][i,:,s] .+= bigCS * tdwJ[ims][istart:iend,1, s];
                    #R_SFI[i,:,s] += bigCS * J₀⁻[istart:iend,1, s];
                    #T_SFI[i,:,s] += bigCS * J₀⁺[istart:iend,1, s];
                    #@show i, s, R_SFI[i,:,s]
                    #@show i, s, ieR_SFI[i,:,s]
                    for t =1:length(topieJ₀⁺[1][1,1,1,:])
                        uwieJ[ims][i,:,s] .+= bigCS * tuwieJ[ims][istart:iend,1, s, t];
                        dwieJ[ims][i,:,s] .+= bigCS * tdwieJ[ims][istart:iend,1, s, t];
                        #@show i, s, t, ieR_SFI[i,:,s]
                    end
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

"Perform post-processing to azimuthally-weight RT matrices"
function postprocessing_vza_ms_canopy!(RS_type::noRS, 
        sensor_levels, 
        quad_points,
        iμ₀, pol_type, 
        composite_layer, 
        solJ₀,
        vza, qp_μ, m, vaz, μ₀, weight, 
        nSpec, SFI, 
        uwJ, dwJ, uwieJ, dwieJ, 
        hdr_J₀⁻, bhr_uw, bhr_dw,
        I_static::AbstractArray{FT2}, 
        arr_type) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}
    
    # idfx of μ0 = cos(sza)
    st_iμ0, istart0, iend0 = get_indices(iμ₀, pol_type);
    
    # Convert these to Arrays (if CuArrays), so they can be accessed by index
    topR⁺⁻ = composite_layer.topR⁺⁻;
    botR⁻⁺ = composite_layer.botR⁻⁺;
    topJ₀⁺ = composite_layer.topJ₀⁺;
    topJ₀⁻ = composite_layer.topJ₀⁻;
    botJ₀⁺ = composite_layer.botJ₀⁺;
    botJ₀⁻ = composite_layer.botJ₀⁻;
    # @show FT2
    tuwJ = [arr_type(zeros(FT2, (size(topJ₀⁺[1],1), 1, nSpec))) for i=1:length(sensor_levels)] #similar(topJ₀⁺); #deepcopy(topJ₀⁺)
    tdwJ = [arr_type(zeros(FT2, (size(topJ₀⁺[1],1), 1, nSpec))) for i=1:length(sensor_levels)]#similar(topJ₀⁺); #deepcopy(topJ₀⁺)
    #@show size(tuwJ[1])
    for ims = 1:length(sensor_levels)
        if(sensor_levels[ims]==0)
            tuwJ[ims] .= botJ₀⁻[ims]
            tdwJ[ims] .= botJ₀⁺[ims]
            #if(m==2)
            #    @show botJ₀⁺[ims][:,1,end]
            #    @show tdwJ[ims][:,1,end]
            #    showp
            #end
            
        else
            #@show size(topR⁺⁻[ims]), size(botR⁻⁺[ims]),size(topJ₀⁺[ims]),size(botJ₀⁻[ims])#,size(tdwJ[ims]),size(tuwJ[ims])
            compute_interlayer_flux!(RS_type,
                                    I_static, 
                                    topR⁺⁻[ims], botR⁻⁺[ims], 
                                    topJ₀⁺[ims], botJ₀⁻[ims],
                                    tdwJ[ims], tuwJ[ims], arr_type);

            #=if(m==0)
                @show topJ₀⁺[ims][:,1,end]
                @show tdwJ[ims][:,1,end]
                @show "========="
                @show botJ₀⁻[ims][:,1,end]
                @show tuwJ[ims][:,1,end]
                #showp
            end=#
                                    
        end

    end
    #@show size(bhr_uw), size(bhr_dw)
    interaction_hdrf_canopy!(SFI, tdwJ[2], tuwJ[2], solJ₀,
                            m, pol_type, quad_points, 
                            hdr_J₀⁻, bhr_uw, bhr_dw)
    
    # Loop over all viewing zenith angles
    for i in eachindex(vza)

        # Find the nearest quadrature point idx
        iμ = nearest_point(qp_μ, cosd(vza[i]));
        st_iμ, istart, iend = get_indices(iμ, pol_type);
        
        # Compute bigCS
        cos_m_phi, sin_m_phi = (cosd(m * vaz[i]), sind(m * vaz[i]));
        bigCS = weight * Diagonal([cos_m_phi, cos_m_phi, sin_m_phi, sin_m_phi][1:pol_type.n]);

        #@show tuwJ[1][istart:iend,1, 1], tdwJ[1][istart:iend,1, 1]
        # Accumulate Fourier moments after azimuthal weighting
        for ims = 1:length(sensor_levels)
            _tuwJ = Array(tuwJ[ims])
            _tdwJ = Array(tdwJ[ims])
            for s = 1:nSpec
                #if true#SFI
                #@show typeof(uwJ[ims]), typeof(tuwJ[ims]), typeof(bigCS)
                uwJ[ims][i,:,s] .+= bigCS * _tuwJ[istart:iend,1, s];
                dwJ[ims][i,:,s] .+= bigCS * _tdwJ[istart:iend,1, s];
                # else
                #     R[i,:,s] += bigCS * (R⁻⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                #     T[i,:,s] += bigCS * (T⁺⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                #end
            end
        end
    end
end