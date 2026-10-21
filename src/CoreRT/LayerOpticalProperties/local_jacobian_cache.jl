"Fourier-independent optical values, retrieval coefficients, and shared scalar seeds."
struct LocalOpticalJacobianCache{L,T}
    layers::L
    basis_tau::T
    basis_omega::T
    microphysics_columns::Vector{Vector{Int}}
end

"""
    build_local_jacobian_cache(band, model, lin_model, cache)

Prepare the scalar chain rule once per solve. Write Sᵢ=τᵢϖᵢ, S=ΣSᵢ,
αᵢ=Sᵢ/S for aerosol i, and use Rayleigh as the fixed phase reference:

    Z = Zᵣ + Σᵢ αᵢ (Zᵢ-Zᵣ)
    dZ = Σᵢ [dαᵢ (Zᵢ-Zᵣ) + αᵢ dZᵢ]
    dαᵢ = (dSᵢ-αᵢ dS)/S,  dϖ = (dS-ϖ dτ)/τ.

These follow by differentiating the mixture definitions S2014 (C.22)–(C.24)
and applying (C.25)–(C.26). The reference choice makes the phase directions
independent of layer height: all layers share one phase basis per Fourier
order. Only small scalar coefficients carry retrieval columns. Gas columns
have dS=dSᵢ=0, hence exactly zero phase coefficients.

The active selection determines the retained microphysical phase directions
for each species. Mixture directions stay present for every species, including
fixed or zero-loading aerosols; their weights can respond to other parameters.

The aerosol cache supplies δ-M-modified values and derivatives, including
both dfᵗ terms in τ/ϖ. Truncated phase tangents supply the remaining dfᵗ
normalization. No truncation chain term is dropped at this boundary.
"""
function build_local_jacobian_cache(band, model, lin_model, cache::LinMInvariantCache)
    iB = band isa Integer ? band : only(band)
    AT = array_type(model)
    FT = eltype(model.τ_rayl[iB])
    na = length(cache.aerosol[1])
    nz = length(cache.rayl_τ_dev[1])
    selection = cache.selection
    pressure = selection === nothing || selection.include_pressure
    columns = [selection === nothing ? collect(1:7) : selection.aerosol_columns[i] for i in 1:na]
    ng = size(cache.lin_gas[1][1].τ̇,2)
    np = Int(pressure) + sum(length,columns; init=0) + ng
    microphysics_columns = local_microphysics_columns(na,selection)
    nb = local_jacobian_size(na,selection)
    mixture_columns = 3 .+ cumsum(vcat(0, [1+length(c) for c in microphysics_columns]))[1:na]
    layers = map(1:nz) do z
        ray_tau = cache.rayl_τ_dev[1][z]
        components = [cache.aerosol[1][i][z] for i in 1:na]
        # Preserve the forward optical algebra's evaluation order. Rebuilding
        # S as a flat sum is mathematically equivalent, but changes Float32 ϖ
        # by an ulp; repeated doubling can amplify that into a radiance change.
        # These spectral weights also reproduce the ordinary successive phase
        # mixtures without constructing retrieval-sized phase tangents.
        τ_scat = copy(ray_tau)
        ϖ_scat = one.(ray_tau)
        mixing_weights = Tuple{typeof(ray_tau),typeof(ray_tau)}[]
        for a in components
            wx = τ_scat .* ϖ_scat
            wy = a.τ .* a.ϖ
            τ_scat = τ_scat .+ a.τ
            S_mix = wx .+ wy
            ϖ_scat = S_mix ./ τ_scat
            # Use division, as in the forward + operator, rather than a
            # reciprocal multiply with a different rounding sequence.
            denominator = ifelse.(S_mix .> zero(FT), S_mix, one(FT))
            push!(mixing_weights,(wx ./ denominator, wy ./ denominator))
        end
        S = τ_scat .* ϖ_scat
        τ = τ_scat .+ cache.gas[1][z].τ
        ϖ = S ./ τ
        invS = ifelse.(S .> zero(FT), one(FT) ./ S, zero(FT))
        weights = [a.τ .* a.ϖ .* invS for a in components]
        dt = _zero_tangent(τ,length(S),np)
        dS = zero(dt)
        component_dS = [zero(dt) for _ in 1:na]
        C = _zero_tangent(τ,length(S),nb,np)
        if pressure
            raydot = _to_device(AT,lin_model.τ̇_rayl_psurf[iB][:,z])
            dt[:,1] .= raydot .+ _to_device(AT,lin_model.τ̇_abs_psurf[iB][:,z])
            dS[:,1] .= raydot
        end
        offset = Int(pressure)
        for (i,a) in enumerate(components)
            ix = offset .+ (1:length(columns[i]))
            b = mixture_columns[i]
            scatterdot = a.τ̇ .* a.ϖ .+ a.τ .* a.ϖ̇
            dt[:,ix] .= a.τ̇
            dS[:,ix] .= scatterdot
            component_dS[i][:,ix] .= scatterdot
            for (j,native) in enumerate(columns[i])
                if 2 <= native <= 5
                    micro = findfirst(==(native-1),microphysics_columns[i])
                    C[:,b+micro,offset+j] .= weights[i]
                end
            end
            if pressure
                optics = model.aerosol_optics[iB][i]
                factor = one(FT) .- optics.fᵗ .* optics.ω̃
                factor = factor isa Number ? factor : _to_device(AT,factor)
                rawdot = _to_device(AT,lin_model.τ̇_aer_psurf[iB][i,:,z])
                tdot = factor .* rawdot
                sdot = tdot .* a.ϖ
                dt[:,1] .+= tdot
                dS[:,1] .+= sdot
                component_dS[i][:,1] .= sdot
            end
            offset += length(columns[i])
        end
        dt[:,offset+1:np] .= cache.lin_gas[1][z].τ̇
        dw = (dS .- ϖ .* dt) ./ τ
        C[:,1,:] .= dt
        C[:,2,:] .= dw
        for i in 1:na
            C[:,mixture_columns[i],:] .= (component_dS[i] .- weights[i] .* dS) .* invS
        end
        (;τ,ϖ,τ̇=dt,ϖ̇=dw,coefficients=C,weights,mixing_weights,rayleigh_fraction=ray_tau ./ τ)
    end
    seed = first(layers).τ
    basis_tau = _zero_tangent(seed,length(seed),nb)
    basis_omega = zero(basis_tau)
    basis_tau[:,1] .= one(FT)
    basis_omega[:,2] .= one(FT)
    return LocalOpticalJacobianCache(layers,basis_tau,basis_omega,microphysics_columns)
end
