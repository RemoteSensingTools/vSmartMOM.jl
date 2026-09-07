"Number of local directions: τ, ϖ, then an aerosol/Rayleigh phase difference and four dZ per aerosol."
local_jacobian_size(n_aerosols::Integer) = 2 + 5n_aerosols

"Select the local basis only when it reduces the atmospheric tangent dimension."
function use_local_jacobian(mode::Symbol, layout, bands, n_aerosols)
    mode in (:auto,:local,:physical) || throw(ArgumentError(
        "jacobian_basis must be :auto, :local, or :physical"))
    mode == :physical && return false
    single_band = bands isa Integer || length(bands) == 1
    mode == :local && !single_band && throw(ArgumentError(
        "local Jacobian assembly currently requires one spectral band per solve"))
    return single_band && (mode == :local || local_jacobian_size(n_aerosols) < n_layer_params(layout))
end

"Prepare a local tangent workspace once, shared across layers and Fourier orders."
function make_local_jacobian_workspace(rs, FT, AT, nb, dims, ns, quad, pol)
    _, tangent = make_added_layer(LinMode(), rs, FT, AT, nb, dims, ns;
        external_solar=quad.external_solar, nStokes=pol.n)
    LocalJacobianWorkspace(tangent, _zero_tangent(AT(FT[]),ns,nb))
end

"""
    local_phase_basis(rayleigh, aerosol_blocks)

Shared matrix directions `[0, 0, Z₁-Zᵣ, dZ₁/dnᵣ₁, …]`, with the four
microphysical directions ordered as `(nᵣ, nᵢ, μ_logr, σ_logr)`, where the size
coordinates are those of `LogNormal(μ_logr, σ_logr)`. The Rayleigh reference
is independent of layer height, so this array is built only once per Fourier
order. See [`build_local_jacobian_cache`](@ref) for the mixture derivative.
Each microphysical tangent differentiates the **truncated** phase matrix.
"""
function local_phase_basis(rayleigh, aerosol_blocks)
    rayleigh === nothing && return nothing
    ray = _ensure_3d(rayleigh)
    ns = isempty(aerosol_blocks) ? size(ray,3) : maximum(size(a[1],3) for a in aerosol_blocks)
    out = _zero_tangent(ray,size(ray,1),size(ray,2),ns,local_jacobian_size(length(aerosol_blocks)))
    for (i,(phase,tangent)) in enumerate(aerosol_blocks)
        b = 3 + 5(i-1)
        out[:,:,:,b] .= _ensure_3d(phase) .- ray
        if tangent !== nothing
            # Scattering: (micro,row,col[,wavelength]); RT: (row,col,wavelength,direction).
            t = ndims(tangent) == 3 ? reshape(tangent,size(tangent)...,1) : tangent
            out[:,:,:,b+1:b+4] .= permutedims(t,(2,3,4,1))
        end
    end
    return out
end

"""
    mix_local_forward_phase(rayleigh, aerosol_phases, mixing_weights)

Build the forward phase with the same successive weighted averages as optical
property `+` (S2014 (C.24)). The derivative basis still uses `Zᵢ-Zᵣ`; using that
subtraction to reconstruct the *forward* mixture introduces avoidable Float32
roundoff. Cached weights let this pass reuse one output array per layer/block.
"""
function mix_local_forward_phase(rayleigh, aerosol_phases, mixing_weights)
    rayleigh === nothing && return nothing
    ray = _ensure_3d(rayleigh)
    isempty(mixing_weights) && return ray
    ns = length(first(mixing_weights)[1])
    out = similar(ray,size(ray,1),size(ray,2),ns)
    out .= ray
    for (aerosol,(wx,wy)) in zip(aerosol_phases,mixing_weights)
        out .= reshape(wx,1,1,ns) .* out .+
               reshape(wy,1,1,ns) .* _ensure_3d(aerosol)
    end
    return out
end

"""
    construct_local_optical_jacobians(rs, band, m, model, lin_model, cache)

Attach one shared phase basis to the Fourier-independent scalar chain rule.
All layers reuse the same local optical seeds. No phase tensor is allocated
with a retrieval-column dimension. After doubling, matrix adding contracts
complete operator tangents; source adding contracts equivalent-source vectors.
"""
function construct_local_optical_jacobians(rs, band, m, model, lin_model,
                                         cache::LocalOpticalJacobianCache)
    InelasticScattering.has_inelastic(rs) && throw(ArgumentError(
        "Linearized Raman-active optical properties are intentionally unsupported."))
    iB = band isa Integer ? band : only(band)
    AT = array_type(model)
    ray = _compute_phase_blocks(model,model.greek_rayleigh[iB],m,AT)
    aeros = [_compute_aerosol_phase_blocks_lin(model,model.aerosol_optics[iB][i],
        lin_model.lin_aerosol_optics[iB][i],get_spec_bands(model)[iB],m,AT)
        for i in eachindex(model.aerosol_optics[iB])]
    phase_basis(k,dk,rk) = local_phase_basis(ray[rk],[(a[k],a[dk]) for a in aeros])
    phase = (phase_basis(1,3,1),phase_basis(2,4,2),phase_basis(5,7,3),phase_basis(6,8,4))
    basis = CoreScatteringOpticalPropertiesLin(cache.basis_tau,cache.basis_omega,phase...)
    forward_phases = map(k->[a[k] for a in aeros],(1,2,5,6))
    layers = map(cache.layers) do layer
        # Aerosol tuples interleave forward/tangent diffuse blocks before the
        # external-solar blocks: (Z++, Z-+, dZ++, dZ-+, Z₀+, Z₀-, ...).
        blocks = ntuple(k->mix_local_forward_phase(ray[k],
            forward_phases[k],layer.mixing_weights),4)
        CoreScatteringOpticalProperties(layer.τ,layer.ϖ,blocks...)
    end
    derivatives = [LocalOpticalJacobian(l.τ̇,l.ϖ̇,basis,l.coefficients) for l in cache.layers]
    fractions = [[l.rayleigh_fraction] for l in cache.layers]
    return layers,derivatives,fractions
end

# Uncached entry point for diagnostics and direct callers.
construct_local_optical_jacobians(rs,band,m,model,lin_model,cache::LinMInvariantCache) =
    construct_local_optical_jacobians(rs,band,m,model,lin_model,
        build_local_jacobian_cache(band,model,lin_model,cache))

expandOpticalProperties(optics::CoreScatteringOpticalProperties,
    jac::LocalOpticalJacobian, AT; expand_Z=false) = (optics,jac)
