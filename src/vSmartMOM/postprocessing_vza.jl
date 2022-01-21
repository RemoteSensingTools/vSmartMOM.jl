#=

This file contains the function to perform azimuthal-weighting to the RT matrices after all 
kernel calculations. 

=#

"Perform post-processing to azimuthally-weight RT matrices"
function postprocessing_vza!(RS_type::noRS, iμ₀, pol_type, composite_layer, vza, qp_μ, m, vaz, μ₀, weight, nSpec, SFI, R, R_SFI, T, T_SFI)
    
    # idx of μ0 = cos(sza)
    st_iμ0, istart0, iend0 = get_indices(iμ₀, pol_type);

    # Convert these to Arrays (if CuArrays), so they can be accessed by index
    R⁻⁺ = Array(composite_layer.R⁻⁺);
    T⁺⁺ = Array(composite_layer.T⁺⁺);
    J₀⁺ = Array(composite_layer.J₀⁺);
    J₀⁻ = Array(composite_layer.J₀⁻);

    # Loop over all viewing zenith angles
    for i = 1:length(vza)

        # Find the nearest quadrature point idx
        iμ = nearest_point(qp_μ, cosd(vza[i]));
        st_iμ, istart, iend = get_indices(iμ, pol_type);
        
        # Compute bigCS
        cos_m_phi, sin_m_phi = (cosd(m * vaz[i]), sind(m * vaz[i]));
        bigCS = weight * Diagonal([cos_m_phi, cos_m_phi, sin_m_phi, sin_m_phi][1:pol_type.n]);

        # Accumulate Fourier moments after azimuthal weighting
        for s = 1:nSpec
            
            if SFI
                R_SFI[i,:,s] += bigCS * J₀⁻[istart:iend,1, s];
                T_SFI[i,:,s] += bigCS * J₀⁺[istart:iend,1, s];
            else
                R[i,:,s] += bigCS * (R⁻⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                T[i,:,s] += bigCS * (T⁺⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
            end
            
        end
    end
end

function postprocessing_vza!(RS_type::Union{RRS, VS_0to1, VS_1to0}, 
        iμ₀, pol_type, composite_layer, 
        vza, qp_μ, m, vaz, μ₀, weight, 
        nSpec, SFI, R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)
    
    # idx of μ0 = cos(sza)
    st_iμ0, istart0, iend0 = get_indices(iμ₀, pol_type);

    # Convert these to Arrays (if CuArrays), so they can be accessed by index
    R⁻⁺ = Array(composite_layer.R⁻⁺);
    T⁺⁺ = Array(composite_layer.T⁺⁺);
    J₀⁺ = Array(composite_layer.J₀⁺);
    J₀⁻ = Array(composite_layer.J₀⁻);

    #ieR⁻⁺ = Array(composite_layer.ieR⁻⁺);
    #ieT⁺⁺ = Array(composite_layer.ieT⁺⁺);
    ieJ₀⁺ = Array(composite_layer.ieJ₀⁺);
    ieJ₀⁻ = Array(composite_layer.ieJ₀⁻);
    # Loop over all viewing zenith angles
    for i = 1:length(vza)

        # Find the nearest quadrature point idx
        iμ = nearest_point(qp_μ, cosd(vza[i]));
        st_iμ, istart, iend = get_indices(iμ, pol_type);
        
        # Compute bigCS
        cos_m_phi, sin_m_phi = (cosd(m * vaz[i]), sind(m * vaz[i]));
        bigCS = weight * Diagonal([cos_m_phi, cos_m_phi, sin_m_phi, sin_m_phi][1:pol_type.n]);

        # Accumulate Fourier moments after azimuthal weighting
        for s = 1:nSpec
            
            if SFI
                R_SFI[i,:,s] += bigCS * J₀⁻[istart:iend,1, s];
                T_SFI[i,:,s] += bigCS * J₀⁺[istart:iend,1, s];
                for t =1:size(ieJ₀⁺,4)# in eachindex ieJ₀⁺[1,1,1,:]
                    ieR_SFI[i,:,s] += bigCS * ieJ₀⁻[istart:iend,1, s, t];
                    ieT_SFI[i,:,s] += bigCS * ieJ₀⁺[istart:iend,1, s, t];
                end
            else
                R[i,:,s] += bigCS * (R⁻⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                T[i,:,s] += bigCS * (T⁺⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
                #no modification for Raman because SFI will be the only mode used for inelastic computations
            end
            
        end
    end
end