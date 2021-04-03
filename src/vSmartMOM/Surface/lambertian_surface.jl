# Create a simple Lambertian albedo layer
function create_surface_layer!(lambertian::LambertianSurfaceScalar{FT}, 
                        added_layer::AddedLayer,#{FT}, 
                        SFI,
                        m::Int,    # Fourier Moment
                        pol_type,  # 
                        quadPoints,
                        τ_sum,
                        architecture) where {FT}
    
    if m == 0
        @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀, μ₀ = quadPoints
        arr_type = array_type(architecture)
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
        # Move to architecture:
        R_surf = arr_type(R_surf)
        if SFI
            I₀_NquadN = similar(qp_μN);
            I₀_NquadN[:] .=0;
            i_start   = pol_type.n*(iμ₀-1) + 1 
            i_end     = pol_type.n*iμ₀;
            I₀_NquadN[i_start:i_end] = pol_type.I₀;
        
            added_layer.J₀⁺[:] .= 0
            added_layer.J₀⁻[:,1,:] = μ₀*(R_surf*I₀_NquadN) .* exp.(-τ_sum/μ₀)';  
            #repeat(repeat(tmp .* exp.(-τ_sum/qp_μN[i_start]), Nquad), 1,1,dim[3])
        end
        R_surf = R_surf * Diagonal(qp_μN.*wt_μN)
        tmp    = ones(pol_type.n*Nquad)
        T_surf = arr_type(Diagonal(tmp))

        
        added_layer.r⁻⁺ .= R_surf;
        added_layer.r⁻⁺ .= 0;
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;

        
    else
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁻⁺[:] .= 0;
        added_layer.t⁺⁺[:] .= 0;
        added_layer.t⁻⁻[:] .= 0;
        added_layer.J₀⁺[:] .= 0;
        added_layer.J₀⁻[:] .= 0;
    end
end