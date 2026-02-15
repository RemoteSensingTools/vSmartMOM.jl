"""
    lin_added_layer_all_params_helper!(RS_type, pol_type, SFI, quad_points,
                                       computed_layer_properties_lin,
                                       added_layer_lin, architecture, ndoubl)

Apply the **chain rule** to map derivatives from the 3 core optical parameters
``(\\tau, \\varpi, \\mathbf{Z})`` to ``N_\\text{params}`` physical state vector parameters.

After `elemental!` and `doubling!`, the `added_layer_lin` contains derivatives of the
RT matrices ``\\mathbf{r}, \\mathbf{t}, \\mathbf{J}_0`` with respect to 3 core parameters:
1. ``\\dot{\\mathbf{M}}[1]``: ``\\partial\\mathbf{M}/\\partial(d\\tau)`` — elemental optical depth
2. ``\\dot{\\mathbf{M}}[2]``: ``\\partial\\mathbf{M}/\\partial\\varpi`` — single-scattering albedo
3. ``\\dot{\\mathbf{M}}[3]``: ``\\partial\\mathbf{M}/\\partial\\mathbf{Z}`` — scattering phase matrix

This function applies the chain rule to obtain derivatives w.r.t. physical parameters ``p_j``:

```math
\\frac{\\partial \\mathbf{r}}{\\partial p_j} = 
  \\frac{\\partial \\mathbf{r}}{\\partial(d\\tau)} \\cdot \\frac{\\partial(d\\tau)}{\\partial p_j} +
  \\frac{\\partial \\mathbf{r}}{\\partial \\varpi} \\cdot \\frac{\\partial \\varpi}{\\partial p_j} +
  \\frac{\\partial \\mathbf{r}}{\\partial \\mathbf{Z}} \\cdot \\frac{\\partial \\mathbf{Z}}{\\partial p_j}
```

!!! note "Elemental optical depth scaling (Bug 14 fix)"
    The core derivatives are w.r.t. the **elemental** optical depth ``d\\tau = \\tau/2^{n_d}``,
    not the full layer ``\\tau``. Therefore the chain rule uses:
    ```math
    \\frac{\\partial(d\\tau)}{\\partial p_j} = \\frac{1}{2^{n_d}} \\frac{\\partial\\tau}{\\partial p_j}
    ```
    The ``\\varpi`` and ``\\mathbf{Z}`` derivatives are used directly since these intensive
    properties do not depend on the layer subdivision.

# Arguments
- `RS_type::noRS`: Raman scattering type (elastic only).
- `pol_type`: Polarization type.
- `SFI`: Source Function Integration flag.
- `quad_points`: Quadrature points (includes ``\\mu_0``, ``i_{\\mu_0}``).
- `computed_layer_properties_lin`: Layer-level ``(\\dot{\\tau}, \\dot{\\varpi}, \\dot{\\mathbf{Z}})``
  with shape `[Nparams × nSpec]`.
- `added_layer_lin`: In/out — core derivatives on input (3 × ...), physical parameter
  derivatives on output (`Nparams × ...`).
- `architecture`: CPU or GPU.
- `ndoubl::Int`: Number of doublings, needed for ``d\\tau = \\tau/2^{n_d}`` scaling.
"""
function lin_added_layer_all_params_helper!(RS_type::noRS{FT}, 
                                pol_type, SFI, quad_points, 
                                computed_layer_properties_lin, 
                                added_layer_lin::AddedLayerLin{FT},
                                architecture, ndoubl::Int) where {FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    @unpack ap_ṙ⁺⁻, ap_ṙ⁻⁺, ap_ṫ⁻⁻, ap_ṫ⁺⁺, ap_J̇₀⁺, ap_J̇₀⁻ = added_layer_lin
    @unpack τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺ = computed_layer_properties_lin
    @unpack D, n = pol_type
    @unpack qp_μ, μ₀, Nquad, iμ₀Nstart = quad_points
    @unpack F₀ = RS_type

    arr_type = array_type(architecture)

    nD=Int(size(Ż⁺⁺,2)/n)
    D_diag = repeat(arr_type(D), nD)             # full diagonal entries
    bigD = Diagonal(D_diag)                     # D-matrix
    
    nparams = size(computed_layer_properties_lin.τ̇)[1]
    nspec = size(computed_layer_properties_lin.τ̇)[2]
    nbigD = size(bigD,1)
    i₀ = iμ₀Nstart:iμ₀Nstart+n-1

    # Compute elemental τ̇: core derivatives are w.r.t. dτ = τ/2^ndoubl
    dτ̇ = τ̇ ./ FT(2^ndoubl)

    Ż⁺⁺_I₀ = arr_type(zeros(nbigD, nspec))
    Ż⁻⁺_I₀ = arr_type(zeros(nbigD, nspec))
    Ż⁺⁺ = arr_type(Ż⁺⁺)
    Ż⁻⁺ = arr_type(Ż⁻⁺)
    F₀ = arr_type(F₀)
    for iparam=1:nparams 
        for ii = 1:nspec
            Ż⁺⁺_I₀[:,ii] = Ż⁺⁺[iparam,:,i₀,ii] * F₀[:,ii]
            Ż⁻⁺_I₀[:,ii] = Ż⁻⁺[iparam,:,i₀,ii] * F₀[:,ii]
        end
        @views ap_ṫ⁺⁺[iparam,:,:,:] .= added_layer_lin.ṫ⁺⁺[1,:,:,:].*reshape(dτ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṫ⁺⁺[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṫ⁺⁺[3,:,:,:].*Ż⁺⁺[iparam,:,:,:] 
        @views ap_ṫ⁻⁻[iparam,:,:,:] .= added_layer_lin.ṫ⁻⁻[1,:,:,:].*reshape(dτ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṫ⁻⁻[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṫ⁻⁻[3,:,:,:].*(reshape(bigD,nbigD,nbigD,1).*Ż⁺⁺[iparam,:,:,:].*reshape(bigD,nbigD,nbigD,1))

        @views ap_ṙ⁻⁺[iparam,:,:,:] .= added_layer_lin.ṙ⁻⁺[1,:,:,:].*reshape(dτ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṙ⁻⁺[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṙ⁻⁺[3,:,:,:].*Ż⁻⁺[iparam,:,:,:]  
        @views ap_ṙ⁺⁻[iparam,:,:,:] .= added_layer_lin.ṙ⁺⁻[1,:,:,:].*reshape(dτ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṙ⁺⁻[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                            added_layer_lin.ṙ⁺⁻[3,:,:,:].*(reshape(bigD,nbigD,nbigD,1).*Ż⁻⁺[iparam,:,:,:].*reshape(bigD,nbigD,nbigD,1))
        if SFI
            @views ap_J̇₀⁺[iparam,:,1,:] .= added_layer_lin.J̇₀⁺[1,:,1,:].*reshape(dτ̇[iparam,:],1,nspec) + 
                                added_layer_lin.J̇₀⁺[2,:,1,:].*reshape(ϖ̇[iparam,:],1,nspec) + 
                                added_layer_lin.J̇₀⁺[3,:,1,:].*Ż⁺⁺_I₀
            @views ap_J̇₀⁻[iparam,:,1,:] .= added_layer_lin.J̇₀⁻[1,:,1,:].*reshape(dτ̇[iparam,:],1,nspec) + 
                                added_layer_lin.J̇₀⁻[2,:,1,:].*reshape(ϖ̇[iparam,:],1,nspec) + 
                                added_layer_lin.J̇₀⁻[3,:,1,:].*Ż⁻⁺_I₀ 
        end
    end
end

"Compute interaction between composite and added layers"
function lin_added_layer_all_params!(RS_type::noRS{FT}, 
    pol_type, SFI, quad_points,
    computed_layer_properties_lin, 
    added_layer_lin::AddedLayerLin{FT}, architecture, ndoubl::Int) where {FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    
    lin_added_layer_all_params_helper!(RS_type, pol_type, 
                    SFI, quad_points,
                    computed_layer_properties_lin, 
                    added_layer_lin, architecture, ndoubl)
    synchronize_if_gpu()
end
