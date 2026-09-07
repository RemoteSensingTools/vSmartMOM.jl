"""
    contract_local_jacobian!(out, basis_tangent, C)

Apply `out[i,j,s,p] = Σ_b basis_tangent[i,j,s,b] C[s,b,p]`. Spectral points remain
independent; matrix indices are preserved. Unlike elemental phase partials,
the local operator tangents already include all angular coupling in doubling.
Any trailing surface/source columns in `out` are set to zero.
"""
function contract_local_jacobian!(out, basis_tangent, C)
    size(basis_tangent,3) == size(C,1) == size(out,3) || throw(DimensionMismatch("spectral basis size"))
    size(basis_tangent,4) == size(C,2) || throw(DimensionMismatch("local basis size"))
    size(out,4) >= size(C,3) || throw(DimensionMismatch("retrieval output size"))
    backend = KernelAbstractions.get_backend(out)
    _contract_local_jacobian!(backend)(out, basis_tangent, C, Val(size(C,2)); ndrange=size(out))
    return out
end

@kernel function _contract_local_jacobian!(out, @Const(basis_tangent), @Const(C), ::Val{NB}) where {NB}
    i,j,s,p = @index(Global, NTuple)
    value = zero(eltype(out))
    if p <= size(C,3)
        @inbounds for b in 1:NB
            coefficient = C[s,b,p]
            # Gas/profile layouts contain exact structural zeros. Skip their
            # matrix reads; this is exact sparsity, with no numerical cutoff.
            if !iszero(coefficient)
                value += basis_tangent[i,j,s,b] * coefficient
            end
        end
    end
    @inbounds out[i,j,s,p] = value
end

"""
    build_doubled_layer_lin!(..., optics, jacobian, ..., workspace)

Construct a layer with the production finite-thickness elemental formulas
(Sanghavi & Frankenberg 2023-II, (10)–(11)) and propagate its supplied tangent
directions through every doubling. Multiple dispatch chooses direct physical
directions or a local optical basis; the forward solver is shared.
"""
function build_doubled_layer_lin!(pol, SFI, above, dabove, dτ, F₀,
        optics, jac::CoreScatteringOpticalPropertiesLin, m, nd, quad,
        added, tangent, arch, identity, workspace=nothing)
    @timeit "elemental" elemental!(pol, SFI, above, dabove, dτ, F₀,
        optics, jac, m, nd, true, quad, added, tangent, arch)
    # Local direct-beam transmission e=exp(-δτ/μ₀), de=-e dδτ/μ₀.
    # Surface/source slots have no atmospheric optical-depth perturbation.
    dt = _zero_tangent(jac.τ̇, length(dτ), size(tangent.ap_ṙ⁻⁺,4))
    dt[:,1:size(jac.τ̇,2)] .= jac.τ̇ ./ eltype(dτ)(2^nd)
    e = exp.(-dτ ./ quad.μ₀)
    @timeit "doubling" doubling_allparams!(pol, SFI, e, nd, added,
        tangent, identity, arch, dt, quad.μ₀; N_active=size(jac.τ̇,2))
    return nothing
end

function build_doubled_layer_lin!(pol, SFI, above, dabove, dτ, F₀,
        optics, jac::LocalOpticalJacobian, m, nd, quad,
        added, tangent, arch, identity, workspace::LocalJacobianWorkspace)
    basis_tangent = workspace.layer
    build_doubled_layer_lin!(pol, SFI, above, workspace.above, dτ, F₀,
        optics, jac.basis, m, nd, quad, added, basis_tangent, arch, identity)
    @timeit "Jacobian contraction" begin
        for name in (:ap_ṙ⁻⁺,:ap_ṙ⁺⁻,:ap_ṫ⁺⁺,:ap_ṫ⁻⁻,:ap_J̇₀⁺,:ap_J̇₀⁻)
            contract_local_jacobian!(getproperty(tangent,name),
                getproperty(basis_tangent,name), jac.coefficients)
        end
        if SFI
            # E=exp(-τ_above/μ₀) is held fixed while doubling local seeds.
            # J=E Ĵ is linear in E, hence dJ=E dĴ-J dτ_above/μ₀.
            # Append this physical-coordinate term exactly once, after local
            # contraction. Local dτ attenuation was already differentiated.
            np = size(dabove,2)
            attenuation = reshape(dabove,1,1,size(dabove,1),np) ./ quad.μ₀
            @views tangent.ap_J̇₀⁺[:,:,:,1:np] .-= added.j₀⁺ .* attenuation
            @views tangent.ap_J̇₀⁻[:,:,:,1:np] .-= added.j₀⁻ .* attenuation
        end
    end
    return nothing
end
