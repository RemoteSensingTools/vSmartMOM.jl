# Algebra experiment only: associative adding with forward prefixes/suffixes.
# Solar-only sources, black lower boundary. No production dispatch is changed.
using KernelAbstractions
include("source_response_operators.jl")
include("source_sweep_operators.jl")
const RESPONSE_FWD = (:R⁻⁺,:R⁺⁻,:T⁺⁺,:T⁻⁻,:J₀⁺,:J₀⁻)
const RESPONSE_ADDED = (:r⁻⁺,:r⁺⁻,:t⁺⁺,:t⁻⁻,:j₀⁺,:j₀⁻)
const RESPONSE_DOT = (:Ṙ⁻⁺,:Ṙ⁺⁻,:Ṫ⁺⁺,:Ṫ⁻⁻,:J̇₀⁺,:J̇₀⁻)
const RESPONSE_ADOT = (:ap_ṙ⁻⁺,:ap_ṙ⁺⁻,:ap_ṫ⁺⁺,:ap_ṫ⁻⁻,:ap_J̇₀⁺,:ap_J̇₀⁻)
function response_copy!(dest,src,dnames=RESPONSE_FWD,snames=RESPONSE_FWD)
    for (d,s) in zip(dnames,snames)
        copyto!(getproperty(dest,d),getproperty(src,s))
    end
end
response_zero!(dest,names=RESPONSE_DOT) = foreach(n->fill!(getproperty(dest,n),0),names)
function response_tangent(layer,nb)
    arrays=map(n->C._zero_tangent(getproperty(layer,n),size(getproperty(layer,n))...,nb),RESPONSE_FWD)
    C.CompositeLayerLin(;NamedTuple{RESPONSE_DOT}(arrays)...)
end
function response_workspace(rs,FT,AT,n,ns,nz,nb,shape,np)
    a,da=C.make_added_layer(LinMode(),rs,FT,AT,nb,(n,n),ns)
    c,dc=C.make_composite_layer(LinMode(),rs,FT,AT,nb,(n,n),ns)
    af=C.make_added_layer(rs,FT,AT,(n,n),ns)
    cf=C.make_composite_layer(rs,FT,AT,(n,n),ns)
    make_layers()=[C.make_composite_layer(rs,FT,AT,(n,n),ns) for _ in 1:nz]
    layers=make_layers(); dlayers=[response_tangent(l,nb) for l in layers]
    vacuum=C.make_composite_layer(rs,FT,AT,(n,n),ns)
    vacuum.T⁺⁺ .= AT(Matrix{FT}(I,n,n)); vacuum.T⁻⁻ .= vacuum.T⁺⁺
    sweep_mode=get(ENV,"AUDIT_RESPONSE_MODE","source")=="sweep"
    source_response=sweep_mode ? nothing : source_response_workspace(first(layers),nb)
    sweep=sweep_mode ? source_sweep_workspace(first(layers),nb,np,nz) : nothing
    (;a,da,c,dc,af,cf,layers,dlayers,vacuum,source_response,sweep,prefix=make_layers(),
      suffix=sweep_mode ? nothing : make_layers(),
      zero_above=AT(zeros(FT,ns,nb)),dR=AT(zeros(FT,shape...,np)),dT=AT(zeros(FT,shape...,np)))
end

@kernel function response_project!(out,@Const(source),@Const(coeff),@Const(rows),@Const(weights))
    v,j,s,p = @index(Global,NTuple)
    value = zero(eltype(out))
    @inbounds for b in axes(coeff,2)
        value += source[rows[v]+j-1,1,s,b]*coeff[s,b,p]
    end
    @inbounds out[v,j,s,p] += weights[v,j]*value
end

"Accumulate independently propagated layer responses in the physical output columns."
function accumulate_layer_responses!(w,data,m,q,identity,rows,weights,mode)
    (;a,da,c,dc,layers,dlayers,prefix,suffix)=w
    nz=length(layers); ns=size(w.dR,3); nb=size(w.zero_above,2)
    interface=C.ScatteringInterface_11()
    for z in 1:nz
        jac=data.local_tangents[m+1][z]
        if mode=="source"
            local_source_response!(w,z,jac,q,identity)
        else
            response_copy!(a,layers[z],RESPONSE_ADDED,RESPONSE_FWD)
            response_copy!(da,dlayers[z],RESPONSE_ADOT,RESPONSE_DOT)
            if z==1
                C.seed_composite_from_added!(c,dc,a,da)
            else
                response_copy!(c,prefix[z-1]);response_zero!(dc)
                C.interaction!(interface,true,c,dc,a,da,identity)
            end
            if z<nz
                response_copy!(a,suffix[z+1],RESPONSE_ADDED,RESPONSE_FWD)
                response_zero!(da,RESPONSE_ADOT)
                attenuation=reshape(jac.basis.τ̇,1,1,ns,nb)./q.μ₀
                da.ap_J̇₀⁺ .= -a.j₀⁺.*attenuation
                da.ap_J̇₀⁻ .= -a.j₀⁻.*attenuation
                C.interaction!(interface,true,c,dc,a,da,identity)
            end
        end
        backend=KernelAbstractions.get_backend(w.dR)
        response_project!(backend)(w.dR,dc.J̇₀⁻,jac.coefficients,rows,weights;ndrange=size(w.dR))
        response_project!(backend)(w.dT,dc.J̇₀⁺,jac.coefficients,rows,weights;ndrange=size(w.dT))
    end
    return nothing
end

"""
Solve the controlled column using local doubled operators and fixed forward
prefixes. `AUDIT_RESPONSE_MODE` selects independent matrix responses, equivalent
source responses with suffixes, or backward/forward source sweeps. Each mode
evaluates the same derivative of the S2014 (23)–(28) affine layer system;
the source helpers document their interface-balance derivations.

All modes include local doubling, forward cache construction and endpoint
projection in the timed call. Workspace allocation and supplied optical inputs
are prepared outside it. These are development experiments with solar-only
sources and a black lower boundary, not production RT entry points.
"""
function layer_responses!(R,T,dR,dT,data,w,rs,pol,q,identity,arch,model,AT)
    (;a,da,af,cf,layers,dlayers,prefix,suffix)=w
    nz=length(layers); ns=size(R,3)
    mode=get(ENV,"AUDIT_RESPONSE_MODE","source")
    mode in ("source","matrix","sweep") || error("Unknown response mode: $mode")
    interface=C.ScatteringInterface_11(); F₀=AT(rs.F₀)
    fill!(R,0);fill!(T,0);fill!(w.dR,0);fill!(w.dT,0)
    for m in 0:2
        for z in 1:nz
            optics=data.values[m+1][z]; jac=data.local_tangents[m+1][z]
            dt,nd=C.get_dtau_ndoubl(optics,q)
            C.build_doubled_layer_lin!(pol,true,data.sums[z],w.zero_above,dt,F₀,
                optics,jac.basis,m,nd,q,a,da,arch,identity)
            response_copy!(layers[z],a,RESPONSE_FWD,RESPONSE_ADDED)
            response_copy!(dlayers[z],da,RESPONSE_DOT,RESPONSE_ADOT)
            if z==1
                response_copy!(cf,layers[z])
            else
                C.interaction!(interface,true,cf,a,identity)
            end
            response_copy!(prefix[z],cf)
        end
        if mode=="sweep"
            sweep_incident_fields!(w,identity)
        else
            response_copy!(suffix[nz],layers[nz])
            for z in nz-1:-1:1
                response_copy!(cf,layers[z]);response_copy!(af,suffix[z+1],RESPONSE_ADDED,RESPONSE_FWD)
                C.interaction!(interface,true,cf,af,identity)
                response_copy!(suffix[z],cf)
            end
        end
        weight=m==0 ? 0.5/π : 1/π
        info=C._precompute_vza_weights(model.obs_geom.vza,model.obs_geom.vaz,q.qp_μ,pol,m,weight)
        rows=AT([i[1] for i in info])
        weights=AT([pol.n==1 ? info[v][3] : info[v][3].diag[j] for v in eachindex(info),j in 1:pol.n])
        if mode=="sweep"
            sweep_source_tangents!(w,data,m,q,rows,weights)
        else
            accumulate_layer_responses!(w,data,m,q,identity,rows,weights,mode)
        end
        C.postprocessing_vza!(rs,q.iμ₀,pol,prefix[end],model.obs_geom.vza,q.qp_μ,m,
            model.obs_geom.vaz,q.μ₀,weight,ns,true,nothing,R,nothing,T,nothing,nothing)
    end
    dR .= Array(w.dR);dT .= Array(w.dT)
    C.reset_timer!()
    nothing
end
