# Create a simple Lambertian albedo layer
function create_surface_layer!(lambertian::LambertianSurfaceScalar{FT}, 
                        added_layer::AddedLayer{FT}, 
                        SFI,
                        m::Int,    # Fourier Moment
                        pol_type,  # 
                        iμ0,
                        qp_μN, wt_qN,
                        τ_sum) where {FT}
    
    if m == 0
        # Albedo normalized by π (and factor 2 for 0th Fourier Moment)
        ρ = 2lambertian.albedo/FT(π)
        # Get size of added layer
        dim = size(added_layer.r⁻⁺)
        Nquad = dim[1] ÷ pol_type.n
        # Ensure matrices are zero:
        tmp = zeros(pol_type.n)
        tmp[1] = ρ   
        R_surf = Array(Diagonal(tmp))
        R_surf = repeat(R_surf',Nquad)
        R_surf = repeat(R_surf',Nquad)
        if SFI
            I₀_NquadN = zeros(FT,dim[1]);
            i_start   = pol_type.n*(iμ0-1) + 1 
            i_end     = pol_type.n*iμ0;
            I₀_NquadN[i_start:i_end] = pol_type.I₀;
        
            added_layer.J₀⁺[:] .= 0
            added_layer.J₀⁻[:,1,:] = qp_μN[i_start]*(R_surf*I₀_NquadN) .* exp.(-τ_sum/qp_μN[i_start])';  
            #repeat(repeat(tmp .* exp.(-τ_sum/qp_μN[i_start]), Nquad), 1,1,dim[3])
        end
        R_surf = R_surf * Diagonal(qp_μN.*wt_qN)
        tmp = ones(pol_type.n*Nquad)
        T_surf = Array(Diagonal(tmp))

        
        added_layer.r⁻⁺[:] = repeat(R_surf,1,1,dim[3]);
        added_layer.r⁻⁺[:] .= 0;
        added_layer.t⁺⁺[:] = repeat(T_surf,1,1,dim[3]);
        added_layer.t⁻⁻[:] = repeat(T_surf,1,1,dim[3]);

        
    else
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁻⁺[:] .= 0;
        added_layer.t⁺⁺[:] .= 0;
        added_layer.t⁻⁻[:] .= 0;
        added_layer.J₀⁺[:] .= 0;
        added_layer.J₀⁻[:] .= 0;
    end


end