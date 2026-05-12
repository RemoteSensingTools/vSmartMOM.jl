#=

This file contains helper functions that are used throughout the linearized vSmartMOM module

=#
"Default matrix in linRT calculation (zeros)"
default_matrix(FT, lin::LinMode, arr_type, Nparams, dims, nSpec)   = arr_type(zeros(FT, tuple(dims[1], dims[2], nSpec, Nparams)))

"Default J matrix in linRT calculation (zeros)"
default_J_matrix(FT, lin::LinMode, arr_type, Nparams, dims, nSpec) = arr_type(zeros(FT, tuple(dims[1], 1, nSpec, Nparams)))



"Make an added layer and its linearized counterpart, supplying all default matrices"
function make_added_layer(lin::LinMode, RS_type::Union{noRS, noRS_plus}, FT, arr_type, Nparams, dims, nSpec) 
    t1 = default_matrix(FT, arr_type, dims, nSpec)
    t2 = default_matrix(FT, arr_type, dims, nSpec)
    t1_ptr = batched_pointer_cache(t1)
    t2_ptr = batched_pointer_cache(t2)
    return AddedLayer(
        r⁻⁺ = default_matrix(FT, arr_type, dims, nSpec), 
        t⁺⁺ = default_matrix(FT, arr_type, dims, nSpec), 
        r⁺⁻ = default_matrix(FT, arr_type, dims, nSpec),
        t⁻⁻ = default_matrix(FT, arr_type, dims, nSpec),
        j₀⁺ = default_J_matrix(FT, arr_type, dims, nSpec),
        j₀⁻ = default_J_matrix(FT, arr_type, dims, nSpec),
        temp1 = t1, temp2 = t2, temp1_ptr = t1_ptr, temp2_ptr = t2_ptr,
        dbl_gp_refl = default_matrix(FT, arr_type, dims, nSpec),
        dbl_j₁⁺ = default_J_matrix(FT, arr_type, dims, nSpec),
        dbl_j₁⁻ = default_J_matrix(FT, arr_type, dims, nSpec),
    ), 
    AddedLayerLin(
        # derivatives wrt τ, ϖ and Z
        ṙ⁻⁺ = default_matrix(FT, lin, arr_type, 3, dims, nSpec),
        ṫ⁺⁺ = default_matrix(FT, lin, arr_type, 3, dims, nSpec),
        ṙ⁺⁻ = default_matrix(FT, lin, arr_type, 3, dims, nSpec),
        ṫ⁻⁻ = default_matrix(FT, lin, arr_type, 3, dims, nSpec),
        J̇₀⁺ = default_J_matrix(FT, lin, arr_type, 3, dims, nSpec),
        J̇₀⁻ = default_J_matrix(FT, lin, arr_type, 3, dims, nSpec),
        # derivatives wrt all parameters
        ap_ṙ⁻⁺ = default_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        ap_ṫ⁺⁺ = default_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        ap_ṙ⁺⁻ = default_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        ap_ṫ⁻⁻ = default_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        ap_J̇₀⁺ = default_J_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        ap_J̇₀⁻ = default_J_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        # Doubling workspace (pre-allocated to avoid per-call allocations)
        dbl_gp_refl_lin    = default_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        dbl_tt_gp_refl_lin = default_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        dbl_ap_expk_lin    = arr_type(zeros(FT, nSpec, Nparams)),
        dbl_J₁⁺            = arr_type(zeros(FT, dims[1], 1, nSpec)),
        dbl_J₁⁻            = arr_type(zeros(FT, dims[1], 1, nSpec)),
        dbl_ap_J̇₁⁺         = default_J_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        dbl_ap_J̇₁⁻         = default_J_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        dbl_gp_refl        = arr_type(zeros(FT, dims[1], dims[2], nSpec)),
        dbl_tt_gp_refl     = arr_type(zeros(FT, dims[1], dims[2], nSpec)),
    )
end

"Make a composite layer and its linearized counterpart, supplying all default matrices"
make_composite_layer(lin::LinMode, RS_type::Union{noRS, noRS_plus}, FT, arr_type, Nparams, dims, nSpec) =
                                        CompositeLayer(
                                            R⁻⁺ = default_matrix(FT, arr_type, dims, nSpec),
                                            R⁺⁻ = default_matrix(FT, arr_type, dims, nSpec),
                                            T⁺⁺ = default_matrix(FT, arr_type, dims, nSpec),
                                            T⁻⁻ = default_matrix(FT, arr_type, dims, nSpec),
                                            J₀⁺ = default_J_matrix(FT, arr_type, dims, nSpec),
                                            J₀⁻ = default_J_matrix(FT, arr_type, dims, nSpec),
                                            # v0.7 Phase A.2a — lin path keeps J₀_by_src empty (no per-source slots
                                            # for now; lin support for thermal lands in Phase A.3).
                                        ),
                                        CompositeLayerLin(
                                            # derivatives wrt all parameters
                                            default_matrix(FT, lin, arr_type, Nparams, dims, nSpec), 
                                            default_matrix(FT, lin,  arr_type, Nparams, dims, nSpec), 
                                            default_matrix(FT, lin,  arr_type, Nparams, dims, nSpec),
                                            default_matrix(FT, lin,  arr_type, Nparams, dims, nSpec),
                                            default_J_matrix(FT, lin,  arr_type, Nparams, dims, nSpec),
                                            default_J_matrix(FT, lin,  arr_type, Nparams, dims, nSpec)
                                        )
