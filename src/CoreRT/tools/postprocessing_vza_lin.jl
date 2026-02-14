#=

This file contains the function to perform azimuthal-weighting to the RT matrices after all 
kernel calculations. 

=#

"Perform post-processing to azimuthally-weight RT matrices"
function postprocessing_vza!(RS_type::noRS, 
                    iμ₀, pol_type, 
                    composite_layer, 
                    composite_layer_lin,
                    vza, qp_μ, m, vaz, μ₀, 
                    weight, nSpec, 
                    SFI, 
                    R_SFI, T_SFI, 
                    Ṙ_SFI, Ṫ_SFI) 
    
    # idx of μ0 = cos(sza)
    st_iμ0, istart0, iend0 = get_indices(iμ₀, pol_type);
    
    
    # Convert these to Arrays (if CuArrays), so they can be accessed by index
    #R⁻⁺ = Array(composite_layer.R⁻⁺);
    #T⁺⁺ = Array(composite_layer.T⁺⁺);
    J₀⁺ = Array(composite_layer.J₀⁺);
    J₀⁻ = Array(composite_layer.J₀⁻);

    #Ṙ⁻⁺ = Array(composite_layer_lin.Ṙ⁻⁺);
    #Ṫ⁺⁺ = Array(composite_layer_lin.Ṫ⁺⁺);
    J̇₀⁺ = Array(composite_layer_lin.J̇₀⁺);
    J̇₀⁻ = Array(composite_layer_lin.J̇₀⁻);
    
    Nparams = size(J̇₀⁻,1)
    # Loop over all viewing zenith angles
    for i = 1:length(vza)

        # Find the nearest quadrature point idx
        iμ = nearest_point(qp_μ, cosd(vza[i]));
        st_iμ, istart, iend = get_indices(iμ, pol_type);
        
        # Compute bigCS
        cos_m_phi, sin_m_phi = (cosd(m * vaz[i]), sind(m * vaz[i]));
        bigCS = weight * Diagonal([cos_m_phi, cos_m_phi, sin_m_phi, sin_m_phi][1:pol_type.n]);

        #@show J₀⁻[istart:iend,1, 1], J₀⁺[istart:iend,1, 1]
        # Accumulate Fourier moments after azimuthal weighting
        #@show vza[i], vaz[i], J₀⁻[istart:iend,1, 1], J₀⁺[istart:iend,1, 1];
        for s = 1:nSpec
            
            #if SFI
            R_SFI[i,:,s] .+= bigCS * J₀⁻[istart:iend,1, s];
            T_SFI[i,:,s] .+= bigCS * J₀⁺[istart:iend,1, s];

            for iparam = 1:Nparams
                Ṙ_SFI[iparam,i,:,s] .+= bigCS * J̇₀⁻[iparam, istart:iend,1, s];
                Ṫ_SFI[iparam,i,:,s] .+= bigCS * J̇₀⁺[iparam, istart:iend,1, s];
            end
            #else
            #    R[i,:,s] .+= bigCS * (R⁻⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
            #    T[i,:,s] .+= bigCS * (T⁺⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀;
            #end
            
        end
    end
end

