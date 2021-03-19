function get_scattering_interface(scattering_interface, scatter, iz)

    # First layer (TOA)
    if (iz == 1)

        # If scattering, 4. If non-scattering, 1. 
        scattering_interface = scatter ? ScatteringInterface_11() : ScatteringInterface_00()
    
    # Not the first layer (not TOA)
    elseif !(scattering_interface isa ScatteringInterface_00)

        # If kn was 1, then toggle between 0/0 and 0/1 
        # Else, toggle between 1/0 and 1/1
        scattering_interface = (scattering_interface isa ScatteringInterface_00) ? 
                                    (!scatter ? ScatteringInterface_00() : ScatteringInterface_01()) : 
                                    (!scatter ? ScatteringInterface_10() : ScatteringInterface_11())
    end

    return scattering_interface
end

"minimum number of doublings needed to reach an optical depth τ_end, starting with an optical depth dτ.
#The starting optical depth dτ is also determined from its maximum possible value, τ"
function doubling_number(dτ_max, τ_end) # check if τ_end can be replaced by τ_end*ϖ for absorbing atmospheres
    FT = eltype(dτ_max)
    # @show FT, eltype(τ_end), dτ_max, τ_end
    # minimum number of doublings needed to reach an optical depth τ_end, starting with an optical depth dτ.
    # The starting optical depth dτ is also determined from its maximum possible value, dτ_max
    if τ_end <= dτ_max
        dτ = τ_end
        ndoubl = 0
        return dτ, ndoubl
    else
        q1 = log10(2.0)
        q2 = log10(dτ_max)
        q3 = log10(τ_end)
        tlimit = (q3 - q2) / q1
        nlimit = floor(Int, tlimit)
        diff = tlimit - nlimit
        if diff < eps(FT)
            dτ = dτ_max
            ndoubl = nlimit
        else
            ndoubl = nlimit + 1       
            x = q3 - q1 * ndoubl
            dτ = 10.0^x
        end 
        return dτ, ndoubl
    end
end

# Finds index i of f_array (i) which is nearest point to f
function nearest_point(f_array, f)
    return argmin(abs.(f_array.-f))
end

function get_indices(iμ::Integer, pol_type::AbstractPolarizationType) 

    st_iμ = (iμ - 1) * pol_type.n
    istart = st_iμ + 1
    iend   = st_iμ + pol_type.n

    return st_iμ, istart, iend
end

default_matrix(FT, arr_type, dims, nSpec)   = arr_type(zeros(FT, tuple(dims[1], dims[2], nSpec)))
default_J_matrix(FT, arr_type, dims, nSpec) = arr_type(zeros(FT, tuple(dims[1], nSpec)))
#default_test_vector(FT, arr_type,3) = arr_type(zeros(FT, tuple(3)))
make_added_layer(FT, arr_type, dims, nSpec)     = AddedLayer(default_matrix(FT, arr_type, dims, nSpec), 
                                                         default_matrix(FT, arr_type, dims, nSpec), 
                                                         default_matrix(FT, arr_type, dims, nSpec),
                                                         default_matrix(FT, arr_type, dims, nSpec),
                                                         default_J_matrix(FT, arr_type, dims, nSpec),
                                                         default_J_matrix(FT, arr_type, dims, nSpec)
                                                         )

make_composite_layer(FT, arr_type, dims, nSpec) = CompositeLayer(default_matrix(FT, arr_type, dims, nSpec), 
                                                                 default_matrix(FT, arr_type, dims, nSpec), 
                                                                 default_matrix(FT, arr_type, dims, nSpec),
                                                                 default_matrix(FT, arr_type, dims, nSpec),
                                                                 default_J_matrix(FT, arr_type, dims, nSpec),
                                                                 default_J_matrix(FT, arr_type, dims, nSpec)
                                                            )