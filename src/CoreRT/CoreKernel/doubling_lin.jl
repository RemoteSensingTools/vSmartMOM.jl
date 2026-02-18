#=
 
This file contains RT doubling-related functions
 
=#

"""
    doubling_helper!(pol_type, SFI, expk, expk_lin, ndoubl, added_layer, 
                     added_layer_lin, I_static, architecture)

Compute full homogeneous layer matrices from the elemental layer using the **Doubling Method**
(de Haan, Bosma & Hovenier 1987), and simultaneously propagate derivatives with respect to
the 3 core optical parameters ``(\\tau, \\varpi, \\mathbf{Z})``.

Starting from the elemental layer (optical depth ``d\\tau = \\tau/2^{n_d}``), this function
doubles the layer ``n_d`` times. After each doubling step, the optical depth doubles and
the reflection/transmission matrices are updated.

# Doubling formulas (de Haan et al. 1987, Eq. 25)

For a homogeneous layer with reflection ``\\mathbf{R}`` and transmission ``\\mathbf{T}``,
the doubled layer has:
```math
\\mathbf{G} = (\\mathbf{I} - \\mathbf{R} \\mathbf{R})^{-1}
```
```math
\\mathbf{R}_{2\\tau} = \\mathbf{R} + \\mathbf{T} \\, \\mathbf{G} \\, \\mathbf{R} \\, \\mathbf{T}
```
```math
\\mathbf{T}_{2\\tau} = \\mathbf{T} \\, \\mathbf{G} \\, \\mathbf{T}
```

# Linearized doubling

The derivatives propagate through the doubling via the product/chain rule. For each core
parameter ``c \\in \\{\\tau, \\varpi, \\mathbf{Z}\\}``:
```math
\\dot{\\mathbf{G}}_c = \\mathbf{G} (\\dot{\\mathbf{R}}_c \\mathbf{R} + 
  \\mathbf{R} \\dot{\\mathbf{R}}_c) \\mathbf{G}
```
```math
\\dot{\\mathbf{R}}_{2\\tau,c} = \\dot{\\mathbf{R}}_c + 
  \\dot{\\mathbf{T}}_c \\mathbf{G} \\mathbf{R} \\mathbf{T} + \\ldots
```

The source function vectors ``\\mathbf{J}_0^\\pm`` are also doubled when `SFI=true`,
with the beam attenuation factor ``e^{-\\tau/\\mu_0}`` applied between doublings.

# Arguments
- `pol_type`: Polarization type.
- `SFI`: Source Function Integration flag.
- `expk`: Beam attenuation factor ``e^{-d\\tau/\\mu_0}`` `[nSpec]`.
- `expk_lin`: Its derivative ``-e^{-d\\tau/\\mu_0}/\\mu_0`` `[nSpec]`.
- `ndoubl::Int`: Number of doubling iterations.
- `added_layer::AddedLayer`: Forward RT matrices (modified in-place).
- `added_layer_lin::AddedLayerLin`: Linearized RT matrices (modified in-place).
- `I_static`: Identity matrix for batched operations.
- `architecture`: CPU or GPU.
"""
function doubling_helper!(pol_type, 
                          SFI, 
                          expk, expk_lin,
                          ndoubl::Int, 
                          added_layer::AddedLayer,
                          added_layer_lin::AddedLayerLin,
                          I_static::AbstractArray{FT}, 
                          architecture) where {FT}

    # Unpack the added layer
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻ = added_layer
    @unpack ṙ⁺⁻, ṙ⁻⁺, ṫ⁻⁻, ṫ⁺⁺, J̇₀⁺, J̇₀⁻ = added_layer_lin
    # Device architecture
    dev = devi(architecture)
    arr_type = array_type(architecture)

    # Note: short-circuit evaluation => return nothing evaluated iff ndoubl == 0 
    ndoubl == 0 && return nothing
    
    # Geometric progression of reflections (1-RR)⁻¹
    gp_refl      = similar(t⁺⁺)
    tt⁺⁺_gp_refl = similar(t⁺⁺)
    gp_refl_lin       = arr_type(zeros(3, size(t⁺⁺)[1], size(t⁺⁺)[2], size(t⁺⁺)[3]))
    tt⁺⁺_gp_refl_lin  = arr_type(zeros(3, size(t⁺⁺)[1], size(t⁺⁺)[2], size(t⁺⁺)[3]))
    if SFI
        # Dummy for source 
        J₁⁺ = similar(j₀⁺)
        J̇₁⁺ = similar(J̇₀⁺)
        # Dummy for J
        J₁⁻ = similar(j₀⁻)
        J̇₁⁻ = similar(J̇₀⁻)
    end

    # Loop over number of doublings
    for n = 1:ndoubl
        
        # T⁺⁺(λ)[I - R⁺⁻(λ)R⁻⁺(λ)]⁻¹, for doubling R⁺⁻,R⁻⁺ and T⁺⁺,T⁻⁻ is identical
        batch_inv!(gp_refl, I_static .- r⁻⁺ ⊠ r⁻⁺)
        tt⁺⁺_gp_refl[:] = t⁺⁺ ⊠ gp_refl
        for iparam = 1:3
            @views gp_refl_lin[iparam,:,:,:] .= gp_refl ⊠ (ṙ⁻⁺[iparam,:,:,:] ⊠ r⁻⁺ .+ r⁻⁺ ⊠ ṙ⁻⁺[iparam,:,:,:]) ⊠ gp_refl 
            @views tt⁺⁺_gp_refl_lin[iparam,:,:,:] .= ṫ⁺⁺[iparam,:,:,:] ⊠ gp_refl .+ t⁺⁺ ⊠ gp_refl_lin[iparam,:,:,:]
        end
        if SFI
            # J⁺₂₁(λ) = J⁺₁₀(λ).exp(-τ(λ)/μ₀)
            @views J₁⁺[:,1,:] = j₀⁺[:,1,:] .* expk'
            # J⁻₁₂(λ)  = J⁻₀₁(λ).exp(-τ(λ)/μ₀)
            @views J₁⁻[:,1,:] = j₀⁻[:,1,:] .* expk'
            for iparam = 1:3
                if iparam == 1
                    @views J̇₁⁺[iparam,:,1,:] .= J̇₀⁺[iparam,:,1,:] .* expk' .+ j₀⁺[:,1,:] .* expk_lin'        
                    @views J̇₁⁻[iparam,:,1,:] .= J̇₀⁻[iparam,:,1,:] .* expk' .+ j₀⁻[:,1,:] .* expk_lin'
                    
                    @views expk_lin .= 2*expk .* expk_lin
                else
                    @views J̇₁⁺[iparam,:,1,:] .= J̇₀⁺[iparam,:,1,:] .* expk'         
                    @views J̇₁⁻[iparam,:,1,:] .= J̇₀⁻[iparam,:,1,:] .* expk' 
                end
                @views J̇₀⁻[iparam,:,:,:] .= J̇₀⁻[iparam,:,:,:] .+ 
                        (tt⁺⁺_gp_refl_lin[iparam,:,:,:] ⊠ (J₁⁻ .+ r⁻⁺ ⊠ j₀⁺)) .+
                        (tt⁺⁺_gp_refl ⊠ (J̇₁⁻[iparam,:,:,:] .+ ṙ⁻⁺[iparam,:,:,:] ⊠ j₀⁺ .+ r⁻⁺ ⊠ J̇₀⁺[iparam,:,:,:]))  
                @views J̇₀⁺[iparam,:,:,:] .= J̇₁⁺[iparam,:,:,:] .+ 
                    (tt⁺⁺_gp_refl_lin[iparam,:,:,:] ⊠ (j₀⁺ .+ r⁻⁺ ⊠ J₁⁻)) .+
                    (tt⁺⁺_gp_refl ⊠ (J̇₀⁺[iparam,:,:,:] .+ ṙ⁻⁺[iparam, :,:,:] ⊠ J₁⁻ .+ r⁻⁺ ⊠ J̇₁⁻[iparam, :,:,:]))
            end

            # J⁻₀₂(λ) = J⁻₀₁(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹[J⁻₁₂(λ) + R⁻⁺₂₁(λ)J⁺₁₀(λ)] (see Eqs.8 in Raman paper draft)
            j₀⁻[:] = j₀⁻ .+ (tt⁺⁺_gp_refl ⊠ (J₁⁻ .+ r⁻⁺ ⊠ j₀⁺)) 
            # J⁺₂₀(λ) = J⁺₂₁(λ) + T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹[J⁺₁₀(λ) + R⁺⁻₀₁(λ)J⁻₁₂(λ)] (see Eqs.8 in Raman paper draft)
            j₀⁺[:] = J₁⁺ .+ (tt⁺⁺_gp_refl ⊠ (j₀⁺ .+ r⁻⁺ ⊠ J₁⁻))
            expk[:] = expk.^2
        end  

        for iparam = 1:3
            ṙ⁻⁺[iparam, :,:,:] .= ṙ⁻⁺[iparam, :,:,:] .+ 
                        tt⁺⁺_gp_refl_lin[iparam, :,:,:] ⊠ r⁻⁺ ⊠ t⁺⁺ .+
                        tt⁺⁺_gp_refl ⊠ (ṙ⁻⁺[iparam,:,:,:] ⊠ t⁺⁺ .+
                        r⁻⁺ ⊠ ṫ⁺⁺[iparam, :,:,:])
            ṫ⁺⁺[iparam, :,:,:]  = tt⁺⁺_gp_refl_lin[iparam, :,:,:] ⊠ t⁺⁺ .+ 
                        tt⁺⁺_gp_refl ⊠ ṫ⁺⁺[iparam, :,:,:]
        end
        # R⁻⁺₂₀(λ) = R⁻⁺₁₀(λ) + T⁻⁻₀₁(λ)[I - R⁻⁺₂₁(λ)R⁺⁻₀₁(λ)]⁻¹R⁻⁺₂₁(λ)T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        r⁻⁺[:]  = r⁻⁺ .+ (tt⁺⁺_gp_refl ⊠ r⁻⁺ ⊠ t⁺⁺)

        # T⁺⁺₂₀(λ) = T⁺⁺₂₁(λ)[I - R⁺⁻₀₁(λ)R⁻⁺₂₁(λ)]⁻¹T⁺⁺₁₀(λ) (see Eqs.8 in Raman paper draft)
        t⁺⁺[:]  = tt⁺⁺_gp_refl ⊠ t⁺⁺
    end

    # After doubling, revert D(DR)->R, where D = Diagonal{1,1,-1,-1}
    # For SFI, after doubling, revert D(DJ₀⁻)->J₀⁻

    synchronize_if_gpu()

    apply_D_matrix!(pol_type.n, added_layer.r⁻⁺, added_layer.t⁺⁺, added_layer.r⁺⁻, added_layer.t⁻⁻)

    SFI && apply_D_matrix_SFI!(pol_type.n, added_layer.j₀⁻)

    return nothing 

end

function doubling!(pol_type, SFI, expk, expk_lin,
                    ndoubl::Int, 
                    added_layer::AddedLayer,#{FT},
                    added_layer_lin::AddedLayerLin,
                    I_static::AbstractArray, 
                    architecture) #where {FT}

    doubling_helper!(pol_type, SFI, expk, expk_lin, 
        ndoubl, added_layer, added_layer_lin, I_static, architecture)
    synchronize_if_gpu()
end

"""
    doubling_allparams_helper!(pol_type, SFI, expk, ndoubl, added_layer, 
                               added_layer_lin, I_static, architecture, dτ̇, μ₀)

Propagate **N physical-parameter** derivatives through the doubling method (Bug 19 fix).

Unlike `doubling_helper!` which propagates only the 3 core derivatives (τ, ϖ, Z),
this function propagates the `ap_` (all-params) fields through doubling. This is 
necessary because the Z chain rule must be applied at the **elemental** level 
(where it is correctly element-wise), not after doubling (where matrix products 
have mixed the Z indices).

The chain rule (`lin_added_layer_all_params!`) should be called BEFORE this function
to fill the `ap_ṙ⁻⁺`, `ap_ṫ⁺⁺`, `ap_J̇₀⁺`, `ap_J̇₀⁻` fields.

For SFI, the beam attenuation derivative `d(e^{-τ/μ₀})/dp_j = -e^{-τ/μ₀}/μ₀ ⋅ ∂τ/∂p_j`
is per-parameter, handled via `dτ̇` (the elemental τ derivative per parameter).
"""
function doubling_allparams_helper!(pol_type, 
                          SFI, 
                          expk,
                          ndoubl::Int, 
                          added_layer::AddedLayer,
                          added_layer_lin::AddedLayerLin,
                          I_static::AbstractArray{FT}, 
                          architecture,
                          dτ̇::AbstractArray,
                          μ₀::FT) where {FT}

    # Unpack the added layer (forward)
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻ = added_layer
    # Use the all-params derivatives
    @unpack ap_ṙ⁺⁻, ap_ṙ⁻⁺, ap_ṫ⁻⁻, ap_ṫ⁺⁺, ap_J̇₀⁺, ap_J̇₀⁻ = added_layer_lin

    dev = devi(architecture)
    arr_type = array_type(architecture)
    
    ndoubl == 0 && return nothing
    
    Nparams = size(ap_ṙ⁻⁺, 1)
    nμ   = size(t⁺⁺, 1)
    nSpec = size(t⁺⁺, 3)
    
    # Temporaries for the geometric progression (forward)
    gp_refl      = similar(t⁺⁺)
    tt⁺⁺_gp_refl = similar(t⁺⁺)
    
    # Temporaries for linearized geometric progression (N params)
    gp_refl_lin      = arr_type(zeros(Nparams, nμ, nμ, nSpec))
    tt⁺⁺_gp_refl_lin = arr_type(zeros(Nparams, nμ, nμ, nSpec))
    
    # Per-parameter beam attenuation derivatives for SFI
    if SFI
        J₁⁺ = similar(j₀⁺)
        ap_J̇₁⁺ = similar(ap_J̇₀⁺)
        J₁⁻ = similar(j₀⁻)
        ap_J̇₁⁻ = similar(ap_J̇₀⁻)
        # Per-parameter expk_lin: d(exp(-dτ/μ₀))/dp_j = -exp(-dτ/μ₀)/μ₀ * dτ̇_j
        ap_expk_lin = arr_type(zeros(Nparams, nSpec))
        for iparam = 1:Nparams
            ap_expk_lin[iparam,:] .= -expk ./ μ₀ .* dτ̇[iparam,:]
        end
    end

    # Loop over number of doublings
    for n = 1:ndoubl
        
        # Forward: geometric progression (1-RR)⁻¹
        batch_inv!(gp_refl, I_static .- r⁻⁺ ⊠ r⁻⁺)
        tt⁺⁺_gp_refl[:] = t⁺⁺ ⊠ gp_refl
        
        # Linearized geometric progression for all N params
        for iparam = 1:Nparams
            @views gp_refl_lin[iparam,:,:,:] .= gp_refl ⊠ (ap_ṙ⁻⁺[iparam,:,:,:] ⊠ r⁻⁺ .+ r⁻⁺ ⊠ ap_ṙ⁻⁺[iparam,:,:,:]) ⊠ gp_refl
            @views tt⁺⁺_gp_refl_lin[iparam,:,:,:] .= ap_ṫ⁺⁺[iparam,:,:,:] ⊠ gp_refl .+ t⁺⁺ ⊠ gp_refl_lin[iparam,:,:,:]
        end
        
        if SFI
            # Forward source doubling
            @views J₁⁺[:,1,:] = j₀⁺[:,1,:] .* expk'
            @views J₁⁻[:,1,:] = j₀⁻[:,1,:] .* expk'
            
            for iparam = 1:Nparams
                # Each parameter has its own beam attenuation derivative
                @views ap_J̇₁⁺[iparam,:,1,:] .= ap_J̇₀⁺[iparam,:,1,:] .* expk' .+ j₀⁺[:,1,:] .* ap_expk_lin[iparam,:]'
                @views ap_J̇₁⁻[iparam,:,1,:] .= ap_J̇₀⁻[iparam,:,1,:] .* expk' .+ j₀⁻[:,1,:] .* ap_expk_lin[iparam,:]'
                
                # Update per-param expk_lin for next doubling: d(expk²)/dp = 2*expk*d(expk)/dp
                @views ap_expk_lin[iparam,:] .= 2 .* expk .* ap_expk_lin[iparam,:]
                
                # Source function doubling (same structure as core, but with ap_ fields)
                @views ap_J̇₀⁻[iparam,:,:,:] .= ap_J̇₀⁻[iparam,:,:,:] .+ 
                        (tt⁺⁺_gp_refl_lin[iparam,:,:,:] ⊠ (J₁⁻ .+ r⁻⁺ ⊠ j₀⁺)) .+
                        (tt⁺⁺_gp_refl ⊠ (ap_J̇₁⁻[iparam,:,:,:] .+ ap_ṙ⁻⁺[iparam,:,:,:] ⊠ j₀⁺ .+ r⁻⁺ ⊠ ap_J̇₀⁺[iparam,:,:,:]))
                @views ap_J̇₀⁺[iparam,:,:,:] .= ap_J̇₁⁺[iparam,:,:,:] .+ 
                    (tt⁺⁺_gp_refl_lin[iparam,:,:,:] ⊠ (j₀⁺ .+ r⁻⁺ ⊠ J₁⁻)) .+
                    (tt⁺⁺_gp_refl ⊠ (ap_J̇₀⁺[iparam,:,:,:] .+ ap_ṙ⁻⁺[iparam,:,:,:] ⊠ J₁⁻ .+ r⁻⁺ ⊠ ap_J̇₁⁻[iparam,:,:,:]))
            end
            
            # Forward source function updates
            j₀⁻[:] = j₀⁻ .+ (tt⁺⁺_gp_refl ⊠ (J₁⁻ .+ r⁻⁺ ⊠ j₀⁺))
            j₀⁺[:] = J₁⁺ .+ (tt⁺⁺_gp_refl ⊠ (j₀⁺ .+ r⁻⁺ ⊠ J₁⁻))
            expk[:] = expk.^2
        end
        
        # Linearized R and T doubling (N params)
        for iparam = 1:Nparams
            ap_ṙ⁻⁺[iparam,:,:,:] .= ap_ṙ⁻⁺[iparam,:,:,:] .+ 
                        tt⁺⁺_gp_refl_lin[iparam,:,:,:] ⊠ r⁻⁺ ⊠ t⁺⁺ .+
                        tt⁺⁺_gp_refl ⊠ (ap_ṙ⁻⁺[iparam,:,:,:] ⊠ t⁺⁺ .+
                        r⁻⁺ ⊠ ap_ṫ⁺⁺[iparam,:,:,:])
            ap_ṫ⁺⁺[iparam,:,:,:] = tt⁺⁺_gp_refl_lin[iparam,:,:,:] ⊠ t⁺⁺ .+ 
                        tt⁺⁺_gp_refl ⊠ ap_ṫ⁺⁺[iparam,:,:,:]
        end
        
        # Forward R and T doubling
        r⁻⁺[:]  = r⁻⁺ .+ (tt⁺⁺_gp_refl ⊠ r⁻⁺ ⊠ t⁺⁺)
        t⁺⁺[:]  = tt⁺⁺_gp_refl ⊠ t⁺⁺
    end

    # After doubling, apply D matrix to forward quantities
    synchronize_if_gpu()
    apply_D_matrix!(pol_type.n, added_layer.r⁻⁺, added_layer.t⁺⁺, added_layer.r⁺⁻, added_layer.t⁻⁻)
    SFI && apply_D_matrix_SFI!(pol_type.n, added_layer.j₀⁻)
    
    # Apply D matrix to all-params derivatives  
    # For n_stokes=1: ap_ṙ⁺⁻ = ap_ṙ⁻⁺, ap_ṫ⁻⁻ = ap_ṫ⁺⁺
    # For n_stokes>1: need proper D transformation (sign flips based on Stokes indices)
    if pol_type.n == 1
        ap_ṙ⁺⁻[:] = ap_ṙ⁻⁺
        ap_ṫ⁻⁻[:] = ap_ṫ⁺⁺
    else
        # General Stokes case: apply D transformation per parameter
        n_stokes = pol_type.n
        nD = div(nμ, n_stokes)
        for iparam = 1:Nparams
            for iSpec = 1:nSpec
                for jμ = 1:nμ
                    j_s = mod1(jμ, n_stokes)
                    for iμ = 1:nμ
                        i_s = mod1(iμ, n_stokes)
                        # First negate r⁻⁺ for rows with i_s > 2
                        r_val = ap_ṙ⁻⁺[iparam, iμ, jμ, iSpec]
                        if i_s > 2
                            r_val = -r_val
                        end
                        # Set r⁺⁻ and t⁻⁻ with appropriate sign
                        if (i_s <= 2 && j_s <= 2) || (i_s > 2 && j_s > 2)
                            ap_ṙ⁺⁻[iparam, iμ, jμ, iSpec] = r_val
                            ap_ṫ⁻⁻[iparam, iμ, jμ, iSpec] = ap_ṫ⁺⁺[iparam, iμ, jμ, iSpec]
                        else
                            ap_ṙ⁺⁻[iparam, iμ, jμ, iSpec] = -r_val
                            ap_ṫ⁻⁻[iparam, iμ, jμ, iSpec] = -ap_ṫ⁺⁺[iparam, iμ, jμ, iSpec]
                        end
                        # Also update ap_ṙ⁻⁺ with the negation for i_s > 2
                        if i_s > 2
                            ap_ṙ⁻⁺[iparam, iμ, jμ, iSpec] = -ap_ṙ⁻⁺[iparam, iμ, jμ, iSpec]
                        end
                    end
                end
            end
        end
        # SFI: apply D to ap_J̇₀⁻
        if SFI
            for iparam = 1:Nparams
                for iSpec = 1:nSpec
                    for iμ = 1:nμ
                        i_s = mod1(iμ, n_stokes)
                        if i_s > 2
                            ap_J̇₀⁻[iparam, iμ, 1, iSpec] = -ap_J̇₀⁻[iparam, iμ, 1, iSpec]
                        end
                    end
                end
            end
        end
    end

    return nothing
end

"""
    doubling_allparams!(pol_type, SFI, expk, ndoubl, added_layer, added_layer_lin,
                        I_static, architecture, dτ̇, μ₀)

Wrapper for `doubling_allparams_helper!`. See that function for documentation.
"""
function doubling_allparams!(pol_type, SFI, expk,
                    ndoubl::Int, 
                    added_layer::AddedLayer,
                    added_layer_lin::AddedLayerLin,
                    I_static::AbstractArray, 
                    architecture, dτ̇, μ₀)

    doubling_allparams_helper!(pol_type, SFI, expk,
        ndoubl, added_layer, added_layer_lin, I_static, architecture, dτ̇, μ₀)
    synchronize_if_gpu()
end

# WARNING: make sure the linearized version does not clash with the Raman version
@kernel function apply_D!(n_stokes::Int,  
                        r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻,
                        ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻)
    iμ, jμ, n = @index(Global, NTuple)
    i = mod(iμ, n_stokes)
    j = mod(jμ, n_stokes)
    i12 = (i == 1) || (i == 2)
    j12 = (j == 1) || (j == 2)
    nparams = size(ṙ⁻⁺, 1)

    if !i12
        r⁻⁺[iμ, jμ, n] = -r⁻⁺[iμ, jμ, n]
        for iparam = 1:nparams
            ṙ⁻⁺[iparam, iμ, jμ, n] = -ṙ⁻⁺[iparam, iμ, jμ, n]
        end
    end

    same_block = (i12 && j12) || (!i12 && !j12)
    s = ifelse(same_block, one(eltype(r⁻⁺)), -one(eltype(r⁻⁺)))

    r⁺⁻[iμ, jμ, n] = s * r⁻⁺[iμ, jμ, n]
    t⁻⁻[iμ, jμ, n] = s * t⁺⁺[iμ, jμ, n]
    for iparam = 1:nparams
        ṙ⁺⁻[iparam, iμ, jμ, n] = s * ṙ⁻⁺[iparam, iμ, jμ, n]
        ṫ⁻⁻[iparam, iμ, jμ, n] = s * ṫ⁺⁺[iparam, iμ, jμ, n]
    end

end

@kernel function apply_D_SFI!(n_stokes::Int, J₀⁻, J̇₀⁻)
    iμ, _, n = @index(Global, NTuple)
    i = mod(iμ, n_stokes)
    i12 = (i == 1) || (i == 2)
    if !i12
        J₀⁻[iμ, 1, n] = - J₀⁻[iμ, 1, n] 
        for iparam = 1:size(J̇₀⁻, 1)
            J̇₀⁻[iparam, iμ, 1, n] = -J̇₀⁻[iparam, iμ, 1, n]
        end
    end
end

function apply_D_matrix!(n_stokes::Int, 
        r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, 
        r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3},
        ṙ⁻⁺::AbstractArray{FT,4}, ṫ⁺⁺::AbstractArray{FT,4}, 
        ṙ⁺⁻::AbstractArray{FT,4}, ṫ⁻⁻::AbstractArray{FT,4}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺  
        ṙ⁺⁻[:] = ṙ⁻⁺
        ṫ⁻⁻[:] = ṫ⁺⁺    
        return nothing
    else 
        device = devi(architecture(r⁻⁺))
        applyD_kernel! = apply_D!(device)
        event = applyD_kernel!(n_stokes, 
                                r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, 
                                ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻, 
                                ndrange=size(r⁻⁺));
        ##wait(device, event);
        synchronize_if_gpu();
        return nothing
    end
end

#=function apply_D_matrix!(n_stokes::Int, r⁻⁺::Array{FT,3}, t⁺⁺::Array{FT,3}, r⁺⁻::Array{FT,3}, t⁻⁻::Array{FT,3}) where {FT}
    if n_stokes == 1
        r⁺⁻[:] = r⁻⁺
        t⁻⁻[:] = t⁺⁺
        
        return nothing
    else 
        device = devi(Architectures.CPU())
        applyD_kernel! = apply_D!(device)
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
        #wait(device, event);
        return nothing
    end
end=#

function apply_D_matrix_SFI!(n_stokes::Int, 
                    J₀⁻::AbstractArray{FT,3}, 
                    J̇₀⁻::AbstractArray{FT,4}) where {FT}
    n_stokes == 1 && return nothing
    device = devi(architecture(J₀⁻))
    applyD_kernel! = apply_D_SFI!(device)
    event = applyD_kernel!(n_stokes, J₀⁻, J̇₀⁻, ndrange=size(J₀⁻));
    ##wait(device, event);
    synchronize_if_gpu();
    nothing
end

#=
function apply_D_matrix_SFI!(n_stokes::Int, J₀⁻::Array{FT,3}) where {FT}
    
    n_stokes == 1 && return nothing

    device = devi(architecture(J₀⁻))
    applyD_kernel! = apply_D_SFI!(device)
    event = applyD_kernel!(n_stokes, J₀⁻, ndrange=size(J₀⁻));
    #wait(device, event);
    
    return nothing
end=#
