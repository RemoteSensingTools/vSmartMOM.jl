# Development experiment: equivalent-source tangent solve with one backward
# incident-field pass and one forward source pass. No suffix operators needed.

function source_sweep_workspace(layer, nb, np, nz)
    matrix() = similar(layer.R⁻⁺)
    vector() = similar(layer.J₀⁺)
    tangent(n) = similar(layer.J₀⁺,size(layer.J₀⁺)...,n)
    (; G=[matrix() for _ in 1:nz], scratch=matrix(),
       down=[vector() for _ in 1:nz], up=[vector() for _ in 1:nz],
       current_up=vector(), v1=vector(), v2=vector(),
       local_up=tangent(nb), local_down=tangent(nb),
       force_up=tangent(np), force_down=tangent(np),
       Jup=tangent(np), Jdown=tangent(np), t1=tangent(np), t2=tangent(np))
end

"""
    sweep_incident_fields!(w, identity)

Recover the fixed diffuse illumination on each layer from the bottom upward.
For prefix P above layer L and incident upwelling U at L's lower face,

    G = (I - R_L↑ R_P↓)⁻¹
    u = G (R_L↑ J_P↓ + T_L↑ U + j_L↑)
    D = R_P↓ u + J_P↓.

Here D and U illuminate L, and u becomes the next upward interface field.
The black boundary gives U=0 at the bottom. These are the interface balances
underlying S2014 (23)–(28); G also serves the tangent source pass below.
The ordering is an implementation derivation, not a paper attribution.
"""
function sweep_incident_fields!(w,identity)
    v=w.sweep
    fill!(v.current_up,0)
    for z in length(w.layers):-1:1
        L=w.layers[z]
        P=z==1 ? w.vacuum : w.prefix[z-1]
        C._jacobian_geometric_inverse!(v.G[z],v.scratch,L.R⁻⁺,P.R⁺⁻,identity)
        copyto!(v.up[z],v.current_up)
        C._bmm!(v.v1,L.R⁻⁺,P.J₀⁺)
        C._bmm!(v.v2,L.T⁻⁻,v.current_up)
        v.v1 .+= v.v2 .+ L.J₀⁻
        C._bmm!(v.current_up,v.G[z],v.v1)
        C._bmm!(v.down[z],P.R⁺⁻,v.current_up)
        v.down[z] .+= P.J₀⁺
    end
    return nothing
end

# Contract equivalent local forcing vectors, not doubled operator matrices.
# The local source tangent holds above-layer attenuation fixed; restore its
# physical-column derivative -j_L dτ_above/μ₀ once here (solar sources only).
@kernel function sweep_contract_forcing!(out,@Const(local_force),@Const(coeff),
                                         @Const(source),@Const(above),μ₀)
    i,s,p=@index(Global,NTuple)
    value=zero(eltype(out))
    @inbounds for b in axes(coeff,2)
        c=coeff[s,b,p]
        iszero(c) || (value+=local_force[i,1,s,b]*c)
    end
    @inbounds out[i,1,s,p]=value-source[i,1,s]*above[s,p]/μ₀
end

@kernel function sweep_project!(out,@Const(source),@Const(rows),@Const(weights))
    v,j,s,p=@index(Global,NTuple)
    @inbounds out[v,j,s,p]+=weights[v,j]*source[rows[v]+j-1,1,s,p]
end

"""
    sweep_source_tangents!(w, data, m, q, rows, weights)

At the fixed incident fields D,U, S2014 (C.6) converts each complete local
operator perturbation into equivalent sources:

    f↑ = dR_L↑ D + dT_L↑ U + dj_L↑
    f↓ = dT_L↓ D + dR_L↓ U + dj_L↓.

Contract these vectors to retrieval columns and solve their source-only column
with the forward operators held fixed. For the current prefix source tangents
δJ_P↑,δJ_P↓, the adding equations (S2014 23–28; SF2023-II 12) reduce to

    v = G (R_L↑ δJ_P↓ + f↑)
    δJ_new↑ = δJ_P↑ + T_P↑ v
    δJ_new↓ = f↓ + T_L↓ (δJ_P↓ + R_P↓ v).

Here δJ denotes the accumulated equivalent-source response evaluated with the
full-column incident fields. Matrix tangents occur only in local doubling and the forcing evaluation;
adding carries physical-column vectors. All old prefix sources are consumed
before writeback. The prototype assumes solar-only sources and a black surface.
"""
function sweep_source_tangents!(w,data,m,q,rows,weights)
    v=w.sweep
    fill!(v.Jup,0); fill!(v.Jdown,0)
    for z in eachindex(w.layers)
        L=w.layers[z]; dL=w.dlayers[z]
        P=z==1 ? w.vacuum : w.prefix[z-1]
        jac=data.local_tangents[m+1][z]
        C._jmul!(v.local_up,dL.Ṙ⁻⁺,v.down[z])
        C._jmul!(v.local_up,dL.Ṫ⁻⁻,v.up[z],one(eltype(v.local_up)))
        v.local_up .+= dL.J̇₀⁻
        C._jmul!(v.local_down,dL.Ṫ⁺⁺,v.down[z])
        C._jmul!(v.local_down,dL.Ṙ⁺⁻,v.up[z],one(eltype(v.local_down)))
        v.local_down .+= dL.J̇₀⁺
        backend=KernelAbstractions.get_backend(v.Jup)
        dims=(size(v.Jup,1),size(v.Jup,3),size(v.Jup,4))
        sweep_contract_forcing!(backend)(v.force_up,v.local_up,jac.coefficients,
            L.J₀⁻,data.dsums[z],q.μ₀;ndrange=dims)
        sweep_contract_forcing!(backend)(v.force_down,v.local_down,jac.coefficients,
            L.J₀⁺,data.dsums[z],q.μ₀;ndrange=dims)

        C._jmul!(v.t1,L.R⁻⁺,v.Jdown); v.t1 .+= v.force_up
        C._jmul!(v.t2,v.G[z],v.t1)
        C._jmul!(v.t1,P.T⁻⁻,v.t2); v.Jup .+= v.t1
        C._jmul!(v.t1,P.R⁺⁻,v.t2); v.t1 .+= v.Jdown
        C._jmul!(v.Jdown,L.T⁺⁺,v.t1); v.Jdown .+= v.force_down
    end
    backend=KernelAbstractions.get_backend(w.dR)
    sweep_project!(backend)(w.dR,v.Jup,rows,weights;ndrange=size(w.dR))
    sweep_project!(backend)(w.dT,v.Jdown,rows,weights;ndrange=size(w.dT))
    return nothing
end
