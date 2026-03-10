#=

This file contains helper functions that are used throughout the linearized vSmartMOM module

=#
"Default matrix in linRT calculation (zeros)"
default_matrix(FT, lin::LinMode, arr_type, Nparams, dims, nSpec)   = arr_type(zeros(FT, tuple(dims[1], dims[2], nSpec, Nparams)))

"Default J matrix in linRT calculation (zeros)"
default_J_matrix(FT, lin::LinMode, arr_type, Nparams, dims, nSpec) = arr_type(zeros(FT, tuple(dims[1], 1, nSpec, Nparams)))



"Make an added layer and its linearized counterpart, supplying all default matrices"
make_added_layer(lin::LinMode, RS_type::Union{noRS, noRS_plus}, FT, arr_type, Nparams, dims, nSpec) = 
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
        default_matrix(FT, lin, arr_type, 3, dims, nSpec), 
        default_matrix(FT, lin, arr_type, 3, dims, nSpec), 
        default_matrix(FT, lin, arr_type, 3, dims, nSpec),
        default_matrix(FT, lin, arr_type, 3, dims, nSpec),
        default_J_matrix(FT, lin, arr_type, 3, dims, nSpec),
        default_J_matrix(FT, lin, arr_type, 3, dims, nSpec),
        # derivatives wrt all parameters
        default_matrix(FT, lin, arr_type, Nparams, dims, nSpec), 
        default_matrix(FT, lin, arr_type, Nparams, dims, nSpec), 
        default_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        default_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        default_J_matrix(FT, lin, arr_type, Nparams, dims, nSpec),
        default_J_matrix(FT, lin, arr_type, Nparams, dims, nSpec)
    )

"Make a composite layer and its linearized counterpart, supplying all default matrices"
make_composite_layer(lin::LinMode, RS_type::Union{noRS, noRS_plus}, FT, arr_type, Nparams, dims, nSpec) = 
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
                                            default_matrix(FT, lin, arr_type, Nparams, dims, nSpec), 
                                            default_matrix(FT, lin,  arr_type, Nparams, dims, nSpec), 
                                            default_matrix(FT, lin,  arr_type, Nparams, dims, nSpec),
                                            default_matrix(FT, lin,  arr_type, Nparams, dims, nSpec),
                                            default_J_matrix(FT, lin,  arr_type, Nparams, dims, nSpec),
                                            default_J_matrix(FT, lin,  arr_type, Nparams, dims, nSpec)
                                        )