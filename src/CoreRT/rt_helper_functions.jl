#=

This file contains helper functions that are used throughout the vSmartMOM module

=#

"Given the previous scattering interface and current layer information, return what type of scattering interface is nexts"
function get_scattering_interface(scattering_interface,scatter, iz)

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

"Minimum number of doublings needed to reach an optical depth τ_end, starting with an optical depth dτ.
The starting optical depth dτ is also determined from its maximum possible value, τ"
function doubling_number(dτ_max, τ_end) # check if τ_end can be replaced by τ_end*ϖ for absorbing atmospheres
    FT = eltype(dτ_max)
    #@show dτ_max, τ_end
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
        #@show q1, q2, q3
        tlimit = (q3 - q2) / q1
        nlimit = floor(Int, tlimit)
#@show tlimit, nlimit 
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

"Finds index i of f_array (i) which is nearest point to f"
nearest_point(f_array, f) = argmin(abs.(f_array.-f))

"Get indices scaled according to pol_type"
function get_indices(iμ::Integer, pol_type::AbstractPolarizationType) 

    st_iμ = (iμ - 1) * pol_type.n
    istart = st_iμ + 1
    iend   = st_iμ + pol_type.n

    return st_iμ, istart, iend
end

"Default matrix in RT calculation (zeros)"
default_matrix(FT, arr_type, dims, nSpec)   = arr_type(zeros(FT, tuple(dims[1], dims[2], nSpec)))
"Default matrix in ieRT calculation (zeros)"
default_matrix_ie(FT, arr_type, dims, nSpec, nRaman)   = arr_type(zeros(FT, tuple(dims[1], dims[2], nSpec, nRaman)))
"Default matrix in linRT calculation (zeros)"
default_matrix_lin(FT, arr_type, Nparams, dims, nSpec)   = arr_type(zeros(FT, tuple(Nparams, dims[1], dims[2], nSpec)))

"Default J matrix in RT calculation (zeros)"
default_J_matrix(FT, arr_type, dims, nSpec) = arr_type(zeros(FT, tuple(dims[1], 1, nSpec)))
"Default J matrix in ieRT calculation (zeros)"
default_J_matrix_ie(FT, arr_type, dims, nSpec, nRaman) = arr_type(zeros(FT, tuple(dims[1], 1, nSpec, nRaman)))
"Default J matrix in linRT calculation (zeros)"
default_J_matrix_lin(FT, arr_type, Nparams, dims, nSpec) = arr_type(zeros(FT, tuple(Nparams, dims[1], 1, nSpec)))

"Default matrix in RT calculation (zeros)"
default_matrix(FT, arr_type, NSens, dims, nSpec)   = [zeros(FT, (dims[1], dims[2], nSpec)) for i=1:NSens]
#arr_type(zeros(FT, tuple(NSens, dims[1], dims[2], nSpec)))
"Default matrix in ieRT calculation (zeros)"
default_matrix_ie(FT, arr_type, NSens, dims, nSpec, nRaman)   = [zeros(FT, (dims[1], dims[2], nSpec, nRaman)) for i=1:NSens]
#zeros(FT, tuple(NSens, dims[1], dims[2], nSpec, nRaman)))

"Default J matrix in RT calculation (zeros)"
default_J_matrix(FT, arr_type, NSens, dims, nSpec) = [zeros(FT, (dims[1], 1, nSpec)) for i=1:NSens]
#arr_type(zeros(FT, tuple(NSens, dims[1], 1, nSpec)))
"Default J matrix in ieRT calculation (zeros)"
default_J_matrix_ie(FT, arr_type, NSens, dims, nSpec, nRaman) = [zeros(FT, (dims[1], 1, nSpec, nRaman)) for i=1:NSens]
#arr_type(zeros(FT, tuple(NSens, dims[1], 1, nSpec, nRaman)))

##### Only for testing, random matrices:
"Default matrix in RT calculation (random)"
default_matrix_rand(FT, arr_type, dims, nSpec)   = arr_type(randn(FT, tuple(dims[1], dims[2], nSpec)))

"Default J matrix in RT calculation (random)"
default_J_matrix_rand(FT, arr_type, dims, nSpec) = arr_type(randn(FT, tuple(dims[1], 1, nSpec)))


"Make an added layer, supplying all default matrices"
make_added_layer(RS_type::Union{noRS, noRS_plus}, FT, arr_type, dims, nSpec) = AddedLayer(
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec)
                                                        )
"Make an added layer and its linearized counterpart, supplying all default matrices"
make_added_layer(RS_type::Union{noRS, noRS_plus}, FT, arr_type, Nparams, dims, nSpec) = 
                                                    AddedLayer(
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec)
                                                    ), 
                                                    AddedLayerLin(
                                                        # derivatives wrt τ, ϖ and Z
                                                        default_matrix_lin(FT, arr_type, 3, dims, nSpec), 
                                                        default_matrix_lin(FT, arr_type, 3, dims, nSpec), 
                                                        default_matrix_lin(FT, arr_type, 3, dims, nSpec),
                                                        default_matrix_lin(FT, arr_type, 3, dims, nSpec),
                                                        default_J_matrix_lin(FT, arr_type, 3, dims, nSpec),
                                                        default_J_matrix_lin(FT, arr_type, 3, dims, nSpec),
                                                        # derivatives wrt all parameters
                                                        default_matrix_lin(FT, arr_type, Nparams, dims, nSpec), 
                                                        default_matrix_lin(FT, arr_type, Nparams, dims, nSpec), 
                                                        default_matrix_lin(FT, arr_type, Nparams, dims, nSpec),
                                                        default_matrix_lin(FT, arr_type, Nparams, dims, nSpec),
                                                        default_J_matrix_lin(FT, arr_type, Nparams, dims, nSpec),
                                                        default_J_matrix_lin(FT, arr_type, Nparams, dims, nSpec)
                                                    )

"Make an added layer, supplying all default matrices"
make_added_layer(RS_type::Union{RRS, RRS_plus,VS_0to1_plus, VS_1to0_plus}, FT, arr_type, dims, nSpec)  = AddedLayerRS(
                                                default_matrix(FT, arr_type, dims, nSpec), 
                                                default_matrix(FT, arr_type, dims, nSpec), 
                                                default_matrix(FT, arr_type, dims, nSpec),
                                                default_matrix(FT, arr_type, dims, nSpec),
                                                default_J_matrix(FT, arr_type, dims, nSpec),
                                                default_J_matrix(FT, arr_type, dims, nSpec),
                                                default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman), 
                                                default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman), 
                                                default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                default_J_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                default_J_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman)
                                                )
                                                         

"Make a random added layer, supplying all random matrices"
make_added_layer_rand(RS_type::Union{noRS, noRS_plus}, FT, arr_type, dims, nSpec)  = AddedLayer(
                                                        default_matrix_rand(FT, arr_type, dims, nSpec), 
                                                        default_matrix_rand(FT, arr_type, dims, nSpec), 
                                                        default_matrix_rand(FT, arr_type, dims, nSpec),
                                                        default_matrix_rand(FT, arr_type, dims, nSpec),
                                                        default_J_matrix_rand(FT, arr_type, dims, nSpec),
                                                        default_J_matrix_rand(FT, arr_type, dims, nSpec)
                                                        )
                                                         
"Make a composite layer, supplying all default matrices"
make_composite_layer(RS_type::Union{noRS, noRS_plus}, 
    FT, arr_type, dims, nSpec) = CompositeLayer(
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec)
                                                        )

"Make a composite layer and its linearized counterpart, supplying all default matrices"
make_composite_layer(RS_type::Union{noRS, noRS_plus}, FT, arr_type, Nparams, dims, nSpec) = 
                                                    CompositeLayer(
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec)
                                                    ), 
                                                    CompositeLayerLin(
                                                        # derivatives wrt all parameters
                                                        default_matrix_lin(FT, arr_type, Nparams, dims, nSpec), 
                                                        default_matrix_lin(FT, arr_type, Nparams, dims, nSpec), 
                                                        default_matrix_lin(FT, arr_type, Nparams, dims, nSpec),
                                                        default_matrix_lin(FT, arr_type, Nparams, dims, nSpec),
                                                        default_J_matrix_lin(FT, arr_type, Nparams, dims, nSpec),
                                                        default_J_matrix_lin(FT, arr_type, Nparams, dims, nSpec)
                                                    )

"Make a composite layer, supplying all default matrices"
make_composite_layer(RS_type::Union{RRS, RRS_plus,VS_0to1_plus, VS_1to0_plus}, 
    FT, arr_type, dims, nSpec) = CompositeLayerRS(
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec),
                                                        default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman), 
                                                        default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman), 
                                                        default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                        default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                        default_J_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                        default_J_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman)
                                                        )
                                                    
"Make a composite layer, supplying all default matrices"
make_composite_layer(RS_type::Union{noRS, noRS_plus}, 
    FT, arr_type, NSens, dims, nSpec) = 
    CompositeLayerMS(default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec))
"Make a composite layer, supplying all default matrices"
make_composite_layer(RS_type::Union{RRS, RRS_plus, VS_0to1_plus, VS_1to0_plus}, 
    FT, arr_type, NSens, dims, nSpec) = 
    CompositeLayerMSRS(default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman), 
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman), 
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_J_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_J_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman), 
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman), 
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_J_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_J_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman)
                    )
"Given a ComputedAtmosphereProperties object, extract a ComputedLayerProperties object using data from the iz index of all arrays in the ComputedAtmosphereProperties"
function get_layer_properties(computed_atmospheric_properties::ComputedAtmosphereProperties, iz, arr_type)
     @unpack τ_λ_all, ϖ_λ_all, τ_all, ϖ_all, Z⁺⁺_all, Z⁻⁺_all , dτ_max_all, dτ_all, ndoubl_all, dτ_λ_all, expk_all, scatter_all, τ_sum_all, fscattRayl_all,  scattering_interfaces_all = computed_atmospheric_properties

    τ_λ = arr_type(τ_λ_all[:, iz])
    ϖ_λ = arr_type(ϖ_λ_all[:, iz])
    τ   = τ_all[iz]
    ϖ   = ϖ_all[iz]
    Z⁺⁺ = arr_type(Z⁺⁺_all[:,:,iz])
    Z⁻⁺ = arr_type(Z⁻⁺_all[:,:,iz])

    dτ_max = dτ_max_all[iz]
    dτ     = dτ_all[iz]
    ndoubl = ndoubl_all[iz]
    dτ_λ   = arr_type(dτ_λ_all[:, iz])
    expk   = arr_type(expk_all[:, iz])
    scatter = scatter_all[iz]
    τ_sum = arr_type(τ_sum_all[:,iz])
    scattering_interface = scattering_interfaces_all[iz]
    fscattRayl = fscattRayl_all[iz]
    #ϖ_Cabannes = ϖ_Cabannes_all[iz]
    # τ * ϖ should remain constant even though they individually change over wavelength
    # @assert all(i -> (i ≈ τ * ϖ), τ_λ .* ϖ_λ)

    return ComputedLayerProperties(τ_λ, ϖ_λ, τ, ϖ, 
        Z⁺⁺, Z⁻⁺, 
        dτ_max, dτ, 
        ndoubl, 
        dτ_λ, 
        expk, 
        scatter, 
        τ_sum, 
        fscattRayl, 
        scattering_interface)
end

