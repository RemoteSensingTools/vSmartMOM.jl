#=
 
This file implements rt_kernel!, which performs the core RT routines (elemental, doubling, interaction)
 
=#

"""
    rt_kernel!(RS_type, pol_type, SFI, added_layer, added_layer_lin, 
               composite_layer, composite_layer_lin, 
               computed_layer_properties, computed_layer_properties_lin,
               scattering_interface, τ_sum, τ̇_sum, m, quad_points, 
               I_static, architecture, qp_μN, iz)

Core RT kernel for a single atmospheric layer `iz` in the linearized mode.

Orchestrates the four fundamental steps of the Matrix Operator Method for layer `iz`:
1. **Elemental**: Compute single-scattering reflection/transmission matrices and their
   derivatives w.r.t. the 3 core optical parameters ``(\\tau, \\varpi, \\mathbf{Z})``.
2. **Chain Rule** (`lin_added_layer_all_params!`): Map the 3 core derivatives to
   `Nparams` physical parameter derivatives at the **elemental** level, where the
   Z chain rule is correctly element-wise (Bug 19 fix).
3. **Doubling** (`doubling_allparams!`): Build full-layer matrices from the elemental 
   layer by repeated doubling, propagating all `Nparams` derivatives using proper 
   matrix products.
4. **Interaction** (Adding): Combine the current layer with the composite layer above 
   using the appropriate `ScatteringInterface` dispatch (00, 01, 10, or 11).

At the TOA layer (`iz == 1`), the composite layer is initialized directly from the
doubled `ap_*` derivatives. For subsequent layers, the adding method accumulates results
from TOA downward.

!!! note "Bug 19: Z chain rule must be applied before doubling"
    The Z derivative ``\\partial\\mathbf{r}/\\partial\\mathbf{Z}`` is a 4th-rank tensor that is
    diagonal at the elemental level (each ``r_{ij}`` depends only on ``Z_{ij}``). After doubling,
    matrix products mix all indices, so the element-wise chain rule ``\\dot{r}[3] .* \\dot{Z}``
    becomes incorrect. The fix applies the chain rule BEFORE doubling.

# Dispatch
- `noRS`: Elastic scattering only (no Raman).
- Future: `RRS`, `VRS` for inelastic scattering (not yet linearized).

See Sanghavi & Stephens (2013), Eqs. 19–34 for the elemental formulas, and
Sanghavi, Davis & Eldering (2014) for the full linearization framework.
"""
### New update: including towers/airborne sensors
function rt_kernel!(RS_type::noRS{FT}, 
                    pol_type, SFI, 
                    added_layer, added_layer_lin,
                    composite_layer, composite_layer_lin,
                    computed_layer_properties::CoreScatteringOpticalProperties, 
                    computed_layer_properties_lin::CoreScatteringOpticalPropertiesLin,
                    scattering_interface, 
                    τ_sum, τ̇_sum, 
                    m, quad_points, 
                    I_static, 
                    architecture, 
                    qp_μN, iz) where {FT}
    @unpack qp_μ, μ₀, Nquad, iμ₀Nstart = quad_points
    @unpack F₀ = RS_type
    # Just unpack core optical properties from 
    @unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties
    @unpack τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺ = computed_layer_properties_lin
    @unpack D, n = pol_type

    arr_type = array_type(architecture)
    
    if any(isnan, τ̇) || any(isnan, ϖ̇) || any(isnan, Ż⁺⁺) || any(isnan, Ż⁻⁺)
    end
    if any(isinf, τ̇) || any(isinf, ϖ̇) 
    end
    
    nD=Int(size(added_layer.t⁺⁺,1)/n)
    D_diag = repeat(arr_type(D), nD)             # full diagonal entries
    bigD = Diagonal(D_diag)                     # D-matrix
    dτ_max = minimum([maximum(τ .* ϖ), FT(0.001) * minimum(qp_μ)])

    _, ndoubl = doubling_number(dτ_max, maximum(τ .* ϖ))
    scatter = true # edit later
    
    # Compute dτ vector
    dτ = τ ./ 2^ndoubl
    dτ̇ = τ̇ ./ 2^ndoubl

    expk = arr_type(exp.(-dτ /μ₀))
    expk_lin = arr_type(exp.(-dτ /μ₀)*(-1/μ₀))
    
    # If there is scattering, perform the elemental and doubling steps
    if scatter
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        arr_type(τ_sum), arr_type(τ̇_sum), 
                                        dτ, arr_type(F₀),
                                        computed_layer_properties,
                                        m, ndoubl, scatter, quad_points,  
                                        added_layer,  
                                        added_layer_lin,
                                        architecture)
        
        
        # Bug 19 fix: Apply chain rule BEFORE doubling.
        # At the elemental level, ∂r/∂Z is element-wise (diagonal 4th-rank tensor),
        # so the element-wise chain rule Ż.*ṫ[3] is correct here.
        # After doubling, matrix products mix all (i,j) indices, making
        # element-wise multiplication wrong for the Z term.
        @timeit "chain_rule" lin_added_layer_all_params!(
                    RS_type, pol_type,
                    SFI, quad_points, 
                    computed_layer_properties_lin, 
                    added_layer_lin, architecture, ndoubl)

        # Bug 22 fix: Add per-parameter τ̇_sum beam attenuation contribution to SFI.
        # The beam attenuation exp(-τ_sum/μ₀) derivative w.r.t. physical parameter p_j is:
        #   J₀⁺ * (-τ̇_sum[j,:] / μ₀)
        # This was previously (incorrectly) baked into J̇₀⁺[1] using only τ̇_sum[1,:].
        # Now it's applied per-parameter after the chain rule.
        if SFI
            nparams_τ_sum = size(τ̇_sum, 1)
            nspec_here = size(added_layer.j₀⁺, 3)
            for iparam = 1:nparams_τ_sum
                @views added_layer_lin.ap_J̇₀⁺[iparam,:,1,:] .+= 
                    added_layer.j₀⁺[:,1,:] .* reshape(-τ̇_sum[iparam,:] ./ μ₀, 1, nspec_here)
                @views added_layer_lin.ap_J̇₀⁻[iparam,:,1,:] .+= 
                    added_layer.j₀⁻[:,1,:] .* reshape(-τ̇_sum[iparam,:] ./ μ₀, 1, nspec_here)
            end
        end

        # Compute per-parameter elemental τ̇ for beam attenuation in doubling.
        # Pad to full Nparams (ap_* array size) — surface params have zero dτ̇.
        Nparams_ap = size(added_layer_lin.ap_ṙ⁻⁺, 1)
        nparams_τ  = size(τ̇, 1)
        nspec_τ    = size(τ̇, 2)
        dτ̇_allparams = arr_type(zeros(FT, Nparams_ap, nspec_τ))
        dτ̇_allparams[1:nparams_τ, :] .= τ̇ ./ FT(2^ndoubl)
        
        @timeit "doubling" doubling_allparams!(pol_type, SFI, 
                                        expk,
                                        ndoubl, 
                                        added_layer, 
                                        added_layer_lin, 
                                        I_static, architecture,
                                        dτ̇_allparams, μ₀)
        
        if any(isnan, added_layer_lin.ap_ṙ⁻⁺) || any(isnan, added_layer_lin.ap_ṫ⁺⁺)
        end
        
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.j₀⁻[:] .= 0;
        _set_transmission_noscat!(added_layer.t⁺⁺, added_layer.t⁻⁻, τ, qp_μN)
        added_layer_lin.ṙ⁻⁺[:] .= 0;
        added_layer_lin.ṙ⁺⁻[:] .= 0;
        added_layer_lin.J̇₀⁻[:] .= 0;
        temp_lin = collect(exp.(-τ./qp_μN') .* (-1 ./ qp_μN'))
        for iλ = 1:length(τ)
            added_layer_lin.ṫ⁺⁺[1,:,:,iλ] = Diagonal(temp_lin[iλ,:])
            added_layer_lin.ṫ⁻⁻[1,:,:,iλ] = Diagonal(temp_lin[iλ,:])
        end
    end
    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        composite_layer.T⁺⁺ .= added_layer.t⁺⁺
        composite_layer.T⁻⁻ .= added_layer.t⁻⁻
        composite_layer.R⁻⁺ .= added_layer.r⁻⁺
        composite_layer.R⁺⁻ .= added_layer.r⁺⁻
        composite_layer.J₀⁺ .= added_layer.j₀⁺
        composite_layer.J₀⁻ .= added_layer.j₀⁻

        # Bug 19 fix: Copy N-parameter derivatives directly from ap_* fields.
        # Chain rule was applied BEFORE doubling (where Z term is element-wise correct),
        # and doubling_allparams! propagated them using proper matrix products.
        # ap_* arrays include all Nparams (layer + surface); they were zeroed before
        # chain rule, so surface params are correctly zero.
        Nparams_comp = size(composite_layer_lin.Ṫ⁺⁺, 1)
        for iparam = 1:Nparams_comp
            @views composite_layer_lin.Ṫ⁺⁺[iparam,:,:,:] .= added_layer_lin.ap_ṫ⁺⁺[iparam,:,:,:]
            @views composite_layer_lin.Ṫ⁻⁻[iparam,:,:,:] .= added_layer_lin.ap_ṫ⁻⁻[iparam,:,:,:]
            @views composite_layer_lin.Ṙ⁻⁺[iparam,:,:,:] .= added_layer_lin.ap_ṙ⁻⁺[iparam,:,:,:]
            @views composite_layer_lin.Ṙ⁺⁻[iparam,:,:,:] .= added_layer_lin.ap_ṙ⁺⁻[iparam,:,:,:]
            if SFI
                @views composite_layer_lin.J̇₀⁺[iparam,:,:,:] .= added_layer_lin.ap_J̇₀⁺[iparam,:,:,:]
                @views composite_layer_lin.J̇₀⁻[iparam,:,:,:] .= added_layer_lin.ap_J̇₀⁻[iparam,:,:,:]
            end
        end
    # If this is not the TOA, perform the interaction step
    else
        # Bug 19 fix: Chain rule already applied before doubling, so ap_* fields
        # are ready for the interaction step. No need to call lin_added_layer_all_params! here.
        @timeit "interaction" interaction!(scattering_interface, 
                    SFI, 
                    composite_layer, composite_layer_lin, 
                    added_layer, added_layer_lin, 
                    I_static)
    end
end
