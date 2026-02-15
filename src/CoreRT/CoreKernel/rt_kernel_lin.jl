#=
 
This file implements rt_kernel!, which performs the core RT routines (elemental, doubling, interaction)
 
=#
#No Raman (default)
# Perform the Core RT routines (elemental, doubling, interaction)
#=
function rt_kernel!(RS_type::noRS, pol_type, SFI, 
        added_layer, added_layer_lin, 
        composite_layer, composite_layer_lin, 
        computed_layer_properties, computed_layer_properties_lin, 
        m, quad_points, I_static, architecture, qp_μN, iz) 

    @unpack τ_λ, ϖ_λ, τ, ϖ, Z⁺⁺, Z⁻⁺, dτ_max, dτ, ndoubl, dτ_λ, expk, scatter, τ_sum, scattering_interface = computed_layer_properties
    # Downselect the following parameters as appropriate
    @unpack τ̇_λ, ϖ̇_λ, τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺, dτ̇, dτ̇_λ, expk_lin, τ̇_sum = computed_layer_properties_lin
    @unpack F₀ = RS_type

    Nparams = size(τ̇_λ)[1]
    #@show τ, ϖ, dτ_max, ndoubl
    # If there is scattering, perform the elemental and doubling steps
    if scatter  
        @timeit "elemental" elemental!(pol_type, SFI, 
                                        τ_sum, dτ_λ, dτ, 
                                        ϖ_λ, ϖ, 
                                        Z⁺⁺, Z⁻⁺, 
                                        τ̇_sum, dτ̇_λ, dτ̇, 
                                        ϖ̇_λ, ϖ̇, 
                                        Ż⁺⁺, Ż⁻⁺, 
                                        F₀,
                                        m, ndoubl, 
                                        scatter, 
                                        quad_points,  
                                        added_layer,  
                                        added_layer_lin,  
                                        I_static, 
                                        architecture)
        #println("Elemental done...")
        @timeit "doubling"   doubling!(pol_type, SFI, expk, expk_lin, ndoubl, added_layer, added_layer_lin, I_static, architecture)
        #println("Doubling done...")
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.j₀⁻[:] .= 0;
        added_layer_lin.ṙ⁻⁺[:] .= 0;
        added_layer_lin.ṙ⁺⁻[:] .= 0;
        added_layer_lin.J̇₀⁻[:] .= 0;
        temp = Array(exp.(-τ_λ./qp_μN'))
        temp_lin = Array(exp.(-τ_λ./qp_μN') * (-1 ./ qp_μN))
        #added_layer.t⁺⁺, added_layer.t⁻⁻ = (Diagonal(exp(-τ_λ / qp_μN)), Diagonal(exp(-τ_λ / qp_μN)))   
        for iλ = 1:length(τ_λ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);
            
            # let ṫ, ṙ, snf J̇ in each layer be functions only of τ*, ϖ* and Z*, which in turn are composite functions of Nparams individual state parameters   
            added_layer_lin.ṫ⁺⁺[1,:,:,iλ] = Diagonal(temp_lin[iλ,:])
            added_layer_lin.ṫ⁻⁻[1,:,:,iλ] = Diagonal(temp_lin[iλ,:])
            
        end
    end
    #M1 = Array(added_layer.t⁺⁺)
    #M2 = Array(added_layer.r⁺⁻)
    #M3 = Array(added_layer.j₀⁻)
    #M4 = Array(added_layer.j₀⁺)
    #@show M1[1,1,1], M2[1,1,1], M3[1,1,1], M4[1,1,1]
    # @assert !any(isnan.(added_layer.t⁺⁺))
    
    # If this TOA, just copy the added layer into the composite layer
    if (iz == 1)
        composite_layer.T⁺⁺[:], composite_layer.T⁻⁻[:] = (added_layer.t⁺⁺, added_layer.t⁻⁻)
        composite_layer.R⁻⁺[:], composite_layer.R⁺⁻[:] = (added_layer.r⁻⁺, added_layer.r⁺⁻)
        composite_layer.J₀⁺[:], composite_layer.J₀⁻[:] = (added_layer.j₀⁺, added_layer.j₀⁻ )
        
        # zero composite variables first
        composite_layer_lin.Ṫ⁺⁺[:] .= 0
        composite_layer_lin.Ṫ⁻⁻[:] .= 0
        composite_layer_lin.Ṙ⁻⁺[:] .= 0
        composite_layer_lin.Ṙ⁺⁻[:] .= 0
        composite_layer_lin.J̇₀⁺[:] .= 0
        composite_layer_lin.J̇₀⁻[:] .= 0

        for iparam = 1:Nparams
            # the following is placeholder code: check later for 
            # 1. use of dτ̇_λ/dϖ̇_λ vs. dτ̇/dϖ̇
            # 2. dimensions
            composite_layer_lin.Ṫ⁺⁺[iparam,:] += added_layer_lin.ṫ⁺⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁺⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁺⁺[3,:,:,:].*dŻ⁺⁺[iparam] 
            composite_layer_lin.Ṫ⁻⁻[iparam,:] += added_layer_lin.ṫ⁻⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁻⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṫ⁻⁻[3,:,:,:].*dŻ⁻⁻[iparam] 

            composite_layer_lin.Ṙ⁻⁺[iparam,:] += added_layer_lin.ṙ⁻⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁻⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁻⁺[3,:,:,:].*dŻ⁻⁺[iparam]  
            composite_layer_lin.Ṙ⁺⁻[iparam,:] += added_layer_lin.ṙ⁺⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁺⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.ṙ⁺⁻[3,:,:,:].*dŻ⁺⁻[iparam] 

            composite_layer_lin.J̇₀⁺[iparam,:] += added_layer_lin.J̇₀⁺[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁺[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁺[3,:,:,:].*dŻ⁺⁺[iparam] 
            composite_layer_lin.J̇₀⁻[iparam,:] += added_layer_lin.J̇₀⁻[1,:,:,:].*dτ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁻[2,:,:,:].*dϖ̇_λ[iparam] + 
                                                added_layer_lin.J̇₀⁻[3,:,:,:].*dŻ⁻⁺[iparam] 
        end
    # If this is not the TOA, perform the interaction step
    else
        @timeit "lin_added_layer_all_params" lin_added_layer_all_params!(SFI, 
                    computed_layer_properties_lin, 
                    added_layer_lin)
        @timeit "interaction" interaction!(RS_type, scattering_interface, SFI, 
                            computed_layer_properties, computed_layer_properties_lin, 
                            composite_layer, composite_layer_lin, 
                            added_layer, added_layer_lin, 
                            I_static)
    end
end
=#

"""
    rt_kernel!(RS_type, pol_type, SFI, added_layer, added_layer_lin, 
               composite_layer, composite_layer_lin, 
               computed_layer_properties, computed_layer_properties_lin,
               scattering_interface, τ_sum, τ̇_sum, m, quad_points, 
               I_static, architecture, qp_μN, iz)

Core RT kernel for a single atmospheric layer `iz` in the linearized mode.

Orchestrates the three fundamental steps of the Matrix Operator Method for layer `iz`:
1. **Elemental**: Compute single-scattering reflection/transmission matrices and their
   derivatives w.r.t. the 3 core optical parameters ``(\\tau, \\varpi, \\mathbf{Z})``.
2. **Doubling**: Build full-layer matrices from the elemental layer by repeated doubling.
3. **Chain Rule** (`lin_added_layer_all_params!`): Map the 3 core derivatives to
   `Nparams` physical parameter derivatives via:
   ```math
   \\frac{\\partial \\mathbf{r}}{\\partial p_j} = 
     \\dot{\\mathbf{r}}_\\tau \\cdot \\dot{\\tau}_j +
     \\dot{\\mathbf{r}}_\\varpi \\cdot \\dot{\\varpi}_j +
     \\dot{\\mathbf{r}}_Z \\cdot \\dot{\\mathbf{Z}}_j
   ```
4. **Interaction** (Adding): Combine the current layer with the composite layer above 
   using the appropriate `ScatteringInterface` dispatch (00, 01, 10, or 11).

At the TOA layer (`iz == 1`), the chain rule is applied inline and the composite layer
is initialized directly. For subsequent layers, the adding method accumulates results
from TOA downward.

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
        
        if any(isnan, added_layer_lin.ṙ⁻⁺) || any(isnan, added_layer_lin.ṫ⁺⁺) || any(isnan, added_layer_lin.J̇₀⁺) || any(isnan, added_layer_lin.J̇₀⁻)
            for ic in 1:3
                println("    core[$ic]: ṙ⁻⁺=$(any(isnan, added_layer_lin.ṙ⁻⁺[ic,:,:,:])), ṫ⁺⁺=$(any(isnan, added_layer_lin.ṫ⁺⁺[ic,:,:,:])), J̇₀⁺=$(any(isnan, added_layer_lin.J̇₀⁺[ic,:,:,:])), J̇₀⁻=$(any(isnan, added_layer_lin.J̇₀⁻[ic,:,:,:]))")
            end
        end
        
        @timeit "doubling"   doubling!(pol_type, SFI, 
                                        expk, expk_lin, 
                                        ndoubl, 
                                        added_layer, 
                                        added_layer_lin, 
                                        I_static, architecture)
        
        if any(isnan, added_layer_lin.ṙ⁻⁺) || any(isnan, added_layer_lin.ṫ⁺⁺) || any(isnan, added_layer_lin.J̇₀⁺) || any(isnan, added_layer_lin.J̇₀⁻)
        end
    else # This might not work yet on GPU!
        # If not, there is no reflectance. Assign r/t appropriately
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁺⁻[:] .= 0;
        added_layer.j₀⁻[:] .= 0;
        temp = Array(exp.(-τ./qp_μN'))
        added_layer_lin.ṙ⁻⁺[:] .= 0;
        added_layer_lin.ṙ⁺⁻[:] .= 0;
        added_layer_lin.J̇₀⁻[:] .= 0;
        temp_lin = Array(exp.(-τ./qp_μN') .* (-1 ./ qp_μN'))
        for iλ = 1:length(τ)
            added_layer.t⁺⁺[:,:,iλ] = Diagonal(temp[iλ,:]);
            added_layer.t⁻⁻[:,:,iλ] = Diagonal(temp[iλ,:]);

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

        fill!(composite_layer_lin.Ṫ⁺⁺, 0)
        fill!(composite_layer_lin.Ṫ⁻⁻, 0)
        fill!(composite_layer_lin.Ṙ⁻⁺, 0)
        fill!(composite_layer_lin.Ṙ⁺⁻, 0)
        fill!(composite_layer_lin.J̇₀⁺, 0)
        fill!(composite_layer_lin.J̇₀⁻, 0)

        nspec = size(added_layer.t⁺⁺,3)
        nparams = size(τ̇)[1]
        nbigD = size(bigD,1)
        i₀ = iμ₀Nstart:iμ₀Nstart+n-1

        Ż⁺⁺_I₀ = arr_type(zeros(nbigD, nspec))
        Ż⁻⁺_I₀ = arr_type(zeros(nbigD, nspec))
        Ż⁺⁺ = arr_type(Ż⁺⁺)
        Ż⁻⁺ = arr_type(Ż⁻⁺)
        F₀ = arr_type(F₀)
        for iparam = 1:nparams
            for ii = 1:nspec
                Ż⁺⁺_I₀[:,ii] = Ż⁺⁺[iparam,:,i₀,ii] * F₀[:,ii]
                Ż⁻⁺_I₀[:,ii] = Ż⁻⁺[iparam,:,i₀,ii] * F₀[:,ii]
            end

            @views composite_layer_lin.Ṫ⁺⁺[iparam,:,:,:] .= added_layer_lin.ṫ⁺⁺[1,:,:,:].*reshape(dτ̇[iparam,:],1,1,nspec) .+ 
                                                added_layer_lin.ṫ⁺⁺[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                                                added_layer_lin.ṫ⁺⁺[3,:,:,:].*Ż⁺⁺[iparam,:,:,:] 
            @views composite_layer_lin.Ṫ⁻⁻[iparam,:,:,:] .= added_layer_lin.ṫ⁻⁻[1,:,:,:].*reshape(dτ̇[iparam,:],1,1,nspec) .+ 
                                                added_layer_lin.ṫ⁻⁻[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                                                added_layer_lin.ṫ⁻⁻[3,:,:,:].*(reshape(bigD,nbigD,nbigD,1).*Ż⁺⁺[iparam,:,:,:].*reshape(bigD,nbigD,nbigD,1))

            @views composite_layer_lin.Ṙ⁻⁺[iparam,:,:,:] .= added_layer_lin.ṙ⁻⁺[1,:,:,:].*reshape(dτ̇[iparam,:],1,1,nspec) .+ 
                                                added_layer_lin.ṙ⁻⁺[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                                                added_layer_lin.ṙ⁻⁺[3,:,:,:].*Ż⁻⁺[iparam,:,:,:]  
            @views composite_layer_lin.Ṙ⁺⁻[iparam,:,:,:] .= added_layer_lin.ṙ⁺⁻[1,:,:,:].*reshape(dτ̇[iparam,:],1,1,nspec) .+ 
                                                added_layer_lin.ṙ⁺⁻[2,:,:,:].*reshape(ϖ̇[iparam,:],1,1,nspec) .+ 
                                                added_layer_lin.ṙ⁺⁻[3,:,:,:].*(reshape(bigD,nbigD,nbigD,1).*Ż⁻⁺[iparam,:,:,:].*reshape(bigD,nbigD,nbigD,1)) 

            @views composite_layer_lin.J̇₀⁺[iparam,:,1,:] .= added_layer_lin.J̇₀⁺[1,:,1,:].*reshape(dτ̇[iparam,:],1,nspec) .+ 
                                                added_layer_lin.J̇₀⁺[2,:,1,:].*reshape(ϖ̇[iparam,:],1,nspec) .+ 
                                                added_layer_lin.J̇₀⁺[3,:,1,:].*Ż⁺⁺_I₀
            @views composite_layer_lin.J̇₀⁻[iparam,:,1,:] .= added_layer_lin.J̇₀⁻[1,:,1,:].*reshape(dτ̇[iparam,:],1,nspec) .+ 
                                                added_layer_lin.J̇₀⁻[2,:,1,:].*reshape(ϖ̇[iparam,:],1,nspec) .+ 
                                                added_layer_lin.J̇₀⁻[3,:,1,:].*Ż⁻⁺_I₀
        end
        if any(isnan, composite_layer_lin.Ṙ⁻⁺) || any(isnan, composite_layer_lin.Ṫ⁺⁺)
            for ip in 1:nparams
                rn = any(isnan, @view composite_layer_lin.Ṙ⁻⁺[ip,:,:,:])
                tn = any(isnan, @view composite_layer_lin.Ṫ⁺⁺[ip,:,:,:])
                if rn || tn
                    println("    param $ip: Ṙ⁻⁺_nan=$rn, Ṫ⁺⁺_nan=$tn")
                end
            end
        end
    # If this is not the TOA, perform the interaction step
    else
        @timeit "lin_added_layer_all_params" lin_added_layer_all_params!(
                    RS_type::noRS, pol_type,
                    SFI, quad_points, 
                    computed_layer_properties_lin, 
                    added_layer_lin, architecture, ndoubl)

        if any(isnan, added_layer_lin.ap_ṙ⁻⁺) || any(isnan, added_layer_lin.ap_ṫ⁺⁺) || any(isnan, added_layer_lin.ap_J̇₀⁺)
            nparams_check = size(added_layer_lin.ap_ṙ⁻⁺, 1)
            for ip in 1:nparams_check
                rn = any(isnan, @view added_layer_lin.ap_ṙ⁻⁺[ip,:,:,:])
                tn = any(isnan, @view added_layer_lin.ap_ṫ⁺⁺[ip,:,:,:])
                if rn || tn
                    println("    param $ip: ap_ṙ⁻⁺_nan=$rn, ap_ṫ⁺⁺_nan=$tn")
                end
            end
        end

        @timeit "interaction" interaction!(scattering_interface, 
                    SFI, 
                    composite_layer, composite_layer_lin, 
                    added_layer, added_layer_lin, 
                    I_static)
        
        if any(isnan, composite_layer_lin.Ṙ⁻⁺) || any(isnan, composite_layer_lin.Ṙ⁺⁻) || any(isnan, composite_layer_lin.Ṫ⁺⁺)
        end
    end
end
