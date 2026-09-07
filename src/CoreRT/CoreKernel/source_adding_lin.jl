"""
    use_source_adding(mode, model, sources, SFI, local_basis, brdf)

Select the opt-in equivalent-source Jacobian solve. It currently supports
solar illumination with optional surface SIF, Lambertian surfaces and endpoint radiances. Matrix adding
remains the default and the reference implementation for broader source and
observer configurations. Unsupported explicit requests fail before layer-workspace allocation.
"""
function use_source_adding(mode, model, sources, SFI, local_basis, brdf)
    mode in (:matrix, :source) || throw(ArgumentError(
        "jacobian_adding must be :matrix or :source"))
    mode === :matrix && return false
    local_basis && SFI && _source_adding_supported(sources) &&
        isempty(model.obs_geom.sensor_levels) &&
        brdf isa Union{LambertianSurfaceScalar,LambertianSurfaceLegendre} ||
        throw(ArgumentError("jacobian_adding=:source requires a local optical basis, " *
            "SolarBeam/SurfaceSIF sources with SFI, a scalar or Legendre Lambertian surface, and endpoint observers"))
    return true
end

_source_adding_supported(::AbstractSource) = false
_source_adding_supported(::Union{SolarBeam,SurfaceSIF,NoSource}) = true
_source_adding_supported(s::SourceSet) = all(_source_adding_supported, s.sources)

# Snapshots contain only the six completed operators, without doubling/adding
# scratch. Prefix snapshots need only the three fields used in interface balance.
const _SOURCE_LAYER_FIELDS = (:r⁻⁺, :r⁺⁻, :t⁺⁺, :t⁻⁻, :j₀⁺, :j₀⁻)
const _SOURCE_TANGENT_FIELDS = (:ap_ṙ⁻⁺, :ap_ṙ⁺⁻, :ap_ṫ⁺⁺, :ap_ṫ⁻⁻, :ap_J̇₀⁺, :ap_J̇₀⁻)
_source_snapshot(layer, fields) = NamedTuple{fields}(map(n -> similar(getproperty(layer,n)), fields))
function _copy_source_snapshot!(destination, source)
    for name in keys(destination)
        copyto!(getproperty(destination,name), getproperty(source,name))
    end
    return nothing
end

"""
    make_source_adding_workspace(added, tangent, nparams, nsurf, nsif, nz, AT)

Own the local layer tangents, forward prefixes, interface inverses and vector
scratch for one spectral batch. Atmospheric matrix storage scales with the
local optical basis, not retrieval length. Coefficients, forcing vectors and
outputs still scale with `nparams`. Reuse the workspace across Fourier orders.
"""
function make_source_adding_workspace(added, tangent, nparams, nsurf, nsif, nz, AT)
    matrix() = similar(added.r⁻⁺)
    vector() = similar(added.j₀⁺)
    directions(n) = similar(added.j₀⁺, size(added.j₀⁺)..., n)
    prefix() = (; R⁺⁻=matrix(), T⁻⁻=matrix(), J₀⁺=vector())
    prefixes = [prefix() for _ in 1:nz+1]
    fill!(prefixes[1].R⁺⁻,0)
    fill!(prefixes[1].J₀⁺,0)
    # The TOA prefix is vacuum. Install its identity on the caller's backend.
    FT = eltype(added.r⁻⁺)
    prefixes[1].T⁻⁻ .= AT(Matrix{FT}(I,size(added.r⁻⁺,1),size(added.r⁻⁺,1)))
    nb = size(tangent.ap_ṙ⁻⁺,4)
    (; layers=[_source_snapshot(added,_SOURCE_LAYER_FIELDS) for _ in 1:nz],
       tangents=[_source_snapshot(tangent,_SOURCE_TANGENT_FIELDS) for _ in 1:nz],
       prefixes, G=[matrix() for _ in 1:nz+1], scratch=matrix(),
       down=[vector() for _ in 1:nz+1], up=[vector() for _ in 1:nz+1],
       current_up=vector(), v1=vector(), v2=vector(),
       local_up=directions(nb), local_down=directions(nb),
       surface_up=directions(nsurf), surface_down=directions(nsurf),
       surface_solar_up=vector(), surface_solar_down=vector(),
       surface_emission=directions(nsif),
       surface_zero_above=_zero_tangent(added.r⁻⁺,size(added.r⁻⁺,3),0),
       force_up=directions(nparams), force_down=directions(nparams),
       result=(;J̇₀⁻=directions(nparams), J̇₀⁺=directions(nparams)),
       t1=directions(nparams), t2=directions(nparams),
       zero_above=_zero_tangent(added.r⁻⁺,size(added.r⁻⁺,3),nb))
end

"""
    prepare_source_adding_surface!(w, surface, sources, brdf, m, pol, arch)

Retain the reflected direct-solar source before adding surface emission. The
boundary balance (S2014 (33)–(37)) is affine: `j⁻ = j_solar⁻ + 2 SIF₀` at
`m=0`. Only `j_solar⁻` carries `exp(-τ_column/μ₀)`. The SIF amplitude/slope
derivatives enter as vectors `2 ∂SIF₀/∂p`, with the normalization defined in
`surface_source_contribute!`; they need no surface matrix tangents. Clear them
each Fourier order, since isotropic emission contributes only at `m=0`.
"""
function prepare_source_adding_surface!(w, surface, sources, brdf, m, pol, arch)
    copyto!(w.surface_solar_up, surface.j₀⁻)
    copyto!(w.surface_solar_down, surface.j₀⁺)
    fill!(w.surface_emission, 0)
    surface_source_contribute_lin!(sources, brdf, surface,
        (; ap_J̇₀⁻=w.surface_emission), m, pol, arch, axes(w.surface_emission,4))
    return nothing
end

"Double a local layer, retain its complete tangents, and extend the forward prefix."
function source_adding_layer!(w, z, composite, added, tangent, rs, pol,
        optics, jac::LocalOpticalJacobian, above, m, q, identity, arch, interface;
        dτ_max_threshold=nothing, dτ_min_floor=nothing)
    dτ, nd = get_dtau_ndoubl(optics,q; dτ_max_threshold, dτ_min_floor)
    AT = array_type(arch)
    build_doubled_layer_lin!(pol,true,_to_device(AT,above),w.zero_above,dτ,
        _to_device(AT,rs.F₀),optics,jac.basis,m,nd,q,added,tangent,arch,identity)
    _copy_source_snapshot!(w.layers[z],added)
    _copy_source_snapshot!(w.tangents[z],tangent)
    if z == 1
        composite.R⁻⁺ .= added.r⁻⁺
        composite.R⁺⁻ .= added.r⁺⁻
        composite.T⁺⁺ .= added.t⁺⁺
        composite.T⁻⁻ .= added.t⁻⁻
        composite.J₀⁺ .= added.j₀⁺
        composite.J₀⁻ .= added.j₀⁻
    else
        interaction!(interface,true,composite,added,identity)
    end
    _copy_source_snapshot!(w.prefixes[z+1],composite)
    return nothing
end

"""
    source_incident_fields!(w, surface, identity)

Recover the fixed diffuse illumination D,U on each layer, including the
surface as the last affine operator. For prefix P above layer L:

    G = (I - r_L⁻⁺ R_P⁺⁻)⁻¹
    u = G (r_L⁻⁺ J_P⁺ + t_L⁻⁻ U + j_L⁻)
    D = R_P⁺⁻ u + J_P⁺.

The backward pass starts with U=0 below the surface. The surface's t⁺⁺=I
preserves BOA downwelling. These are the interface balances underlying S2014
(23)–(28); this ordering is our derivation, not an algorithm attributed to the
paper. The same G is reused by the tangent source solve.
"""
function source_incident_fields!(w, surface, identity)
    fill!(w.current_up,0)
    nz = length(w.layers)
    for z in nz+1:-1:1
        L = z > nz ? surface : w.layers[z]
        P = w.prefixes[z]
        _jacobian_geometric_inverse!(w.G[z],w.scratch,L.r⁻⁺,P.R⁺⁻,identity)
        copyto!(w.up[z],w.current_up)
        _bmm!(w.v1,L.r⁻⁺,P.J₀⁺)
        _bmm!(w.v2,L.t⁻⁻,w.current_up)
        w.v1 .+= w.v2 .+ L.j₀⁻
        _bmm!(w.current_up,w.G[z],w.v1)
        _bmm!(w.down[z],P.R⁺⁻,w.current_up)
        w.down[z] .+= P.J₀⁺
    end
    return nothing
end

"Apply S2014 (C.6) at fixed incident fields: f↑=dR D+dT U+dj, likewise f↓."
function equivalent_source_forcing!(up, down, tangent, D, U)
    _jmul!(up,tangent.ap_ṙ⁻⁺,D)
    _jmul!(up,tangent.ap_ṫ⁻⁻,U,one(eltype(up)))
    up .+= tangent.ap_J̇₀⁻
    _jmul!(down,tangent.ap_ṫ⁺⁺,D)
    _jmul!(down,tangent.ap_ṙ⁺⁻,U,one(eltype(down)))
    down .+= tangent.ap_J̇₀⁺
    return nothing
end

# The local seeds hold E_above=exp(-τ_above/μ₀) fixed. Contract whole forcing
# vectors, then restore dE_above exactly once. Trailing surface slots are zero.
@kernel function _contract_source_forcing!(out, @Const(local_force), @Const(C),
                                          @Const(source), @Const(above), μ₀)
    i, _, s, p = @index(Global,NTuple)
    value = zero(eltype(out))
    if p <= size(C,3)
        @inbounds for b in axes(C,2)
            c = C[s,b,p]
            iszero(c) || (value += local_force[i,1,s,b]*c)
        end
        @inbounds value -= source[i,1,s]*above[s,p]/μ₀
    end
    @inbounds out[i,1,s,p] = value
end

# Surface construction holds atmospheric attenuation fixed and differentiates
# only albedo coefficients. Scatter those small directions into retrieval slots
# and restore exp(-τ_column/μ₀)'s derivative for atmospheric columns. `source`
# is the reflected solar source saved BEFORE SIF injection. Differentiating
# the full surface source here would spuriously attenuate emitted fluorescence.
@kernel function _expand_surface_forcing!(out, @Const(surface_force), @Const(source),
                                         @Const(above), μ₀, first_surface)
    i, _, s, p = @index(Global,NTuple)
    value = zero(eltype(out))
    if p <= size(above,2)
        @inbounds value = -source[i,1,s]*above[s,p]/μ₀
    end
    b = p - first_surface + 1
    if 1 <= b <= size(surface_force,4)
        @inbounds value += surface_force[i,1,s,b]
    end
    @inbounds out[i,1,s,p] = value
end

@kernel function _append_surface_emission!(out, @Const(emission), first_sif)
    i, _, s, b = @index(Global,NTuple)
    @inbounds out[i,1,s,first_sif+b-1] += emission[i,1,s,b]
end

"""
    source_adding_tangents!(w, surface, surface_lin, jacobians, above, q, AT, surface_columns, sif_columns)

Solve the equivalent-source column with the forward operators fixed. Applying
S2014 (27)–(28) to the forcing vectors gives

    v = G (r_L⁻⁺ δJ_P⁺ + f⁻)
    δJ_new⁻ = δJ_P⁻ + T_P⁻⁻ v
    δJ_new⁺ = f⁺ + t_L⁺⁺ (δJ_P⁺ + R_P⁺⁻ v).

Here δJ is the response to all accumulated forcing evaluated on the full-column
incident fields. Atmospheric forcing is contracted from the local optical basis;
surface matrix forcing uses only its own albedo columns. Its solar attenuation
is appended once during surface-vector expansion, followed by the independent
SIF emission derivatives. Atmospheric tangents see the full solar + SIF incident
fields, so they include attenuation and multiple scattering of SIF. No retrieval-sized matrix
tangents are formed in either atmospheric adding or surface closure.
"""
function source_adding_tangents!(w, surface, surface_lin, jacobians, above, q, AT, surface_columns, sif_columns)
    Jup, Jdown = w.result.J̇₀⁻, w.result.J̇₀⁺
    fill!(Jup,0)
    fill!(Jdown,0)
    nz = length(w.layers)
    backend = KernelAbstractions.get_backend(Jup)
    for z in 1:nz+1
        L = z > nz ? surface : w.layers[z]
        P = w.prefixes[z]
        if z <= nz
            equivalent_source_forcing!(w.local_up,w.local_down,w.tangents[z],w.down[z],w.up[z])
            C = jacobians[z].coefficients
            # Optical-depth prefix arrays can reside on the host.
            dabove = _to_device(AT,@view(above[:,:,z]))
            _contract_source_forcing!(backend)(w.force_up,w.local_up,C,L.j₀⁻,dabove,q.μ₀;ndrange=size(w.force_up))
            _contract_source_forcing!(backend)(w.force_down,w.local_down,C,L.j₀⁺,dabove,q.μ₀;ndrange=size(w.force_down))
        else
            equivalent_source_forcing!(w.surface_up,w.surface_down,surface_lin,w.down[z],w.up[z])
            dabove = _to_device(AT,@view(above[:,:,z]))
            _expand_surface_forcing!(backend)(w.force_up,w.surface_up,w.surface_solar_up,
                dabove,q.μ₀,first(surface_columns);ndrange=size(w.force_up))
            _expand_surface_forcing!(backend)(w.force_down,w.surface_down,w.surface_solar_down,
                dabove,q.μ₀,first(surface_columns);ndrange=size(w.force_down))
            isempty(sif_columns) || _append_surface_emission!(backend)(
                w.force_up,w.surface_emission,first(sif_columns);ndrange=size(w.surface_emission))
        end
        _jmul!(w.t1,L.r⁻⁺,Jdown)
        w.t1 .+= w.force_up
        _jmul!(w.t2,w.G[z],w.t1)
        _jmul!(w.t1,P.T⁻⁻,w.t2)
        Jup .+= w.t1
        _jmul!(w.t1,P.R⁺⁻,w.t2)
        w.t1 .+= Jdown
        _jmul!(Jdown,L.t⁺⁺,w.t1)
        Jdown .+= w.force_down
    end
    return nothing
end
