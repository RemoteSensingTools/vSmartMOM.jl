# Development experiment: convert a local operator perturbation into equivalent
# source vectors, then propagate those vectors through fixed forward subcolumns.
# S2014 (C.6) applied to the layer's affine input/output map gives the forcing;
# the two geometric inverses below are the ordinary (23)–(28) adding resolvents.

function source_response_workspace(layer, nb)
    matrix() = similar(layer.R⁻⁺)
    vector() = similar(layer.J₀⁺)
    tangent() = similar(layer.J₀⁺,size(layer.J₀⁺)...,nb)
    (; G=matrix(), H=matrix(), scratch=matrix(),
       down=vector(), up=vector(), incident_down=vector(), incident_up=vector(),
       v1=vector(), v2=vector(),
       force_up=tangent(), force_down=tangent(), below_up=tangent(), below_down=tangent(),
       rhs=tangent(), u=tangent(), d=tangent(), t1=tangent(), t2=tangent())
end

"""
Evaluate one layer's endpoint source response without matrix tangents in adding.
Let P be the prefix above L, S the suffix below L, and B=P⊕L. At the fixed
forward solution, D and U are the incident down/up fields on L. The local
operator perturbation is equivalent to two source perturbations:

    f↑ = dR_L D + dT_L↑ U + dj_L↑
    f↓ = dT_L↓ D + dR_L↓ U + dj_L↓.

The below-layer solar source perturbations are s±=-J_S± dτ_L/μ₀. Define
G=(I-R_L↑ R_P↓)⁻¹ and H=(I-R_B↓ R_S↑)⁻¹. The interface responses satisfy

    r = f↑ + T_L↑ s↑
    d = H [T_L↓ R_P↓ G r + f↓ + R_L↓ s↑]
    u = G [r + T_L↑ R_S↑ d].

The endpoint responses are T_P↑ u and T_S↓ d+s↓. These formulas follow by
substitution in the two interface balance equations; all inverses and all
matrix products without dots are forward-only. Only matrix-vector products
carry local directions. No claim is made that the papers use this ordering.
"""
function local_source_response!(w,z,jac,q,identity)
    L=w.layers[z]; dL=w.dlayers[z]; B=w.prefix[z]
    P=z==1 ? w.vacuum : w.prefix[z-1]
    S=z==length(w.layers) ? w.vacuum : w.suffix[z+1]
    v=w.source_response
    C._jacobian_geometric_inverse!(v.G,v.scratch,L.R⁻⁺,P.R⁺⁻,identity)
    C._jacobian_geometric_inverse!(v.H,v.scratch,B.R⁺⁻,S.R⁻⁺,identity)

    # Solve the fixed incident fields at the two faces of L.
    C._bmm!(v.v1,B.R⁺⁻,S.J₀⁻); v.v1 .+= B.J₀⁺
    C._bmm!(v.down,v.H,v.v1)
    C._bmm!(v.incident_up,S.R⁻⁺,v.down); v.incident_up .+= S.J₀⁻
    C._bmm!(v.v1,L.R⁻⁺,P.J₀⁺)
    C._bmm!(v.v2,L.T⁻⁻,v.incident_up); v.v1 .+= v.v2 .+ L.J₀⁻
    C._bmm!(v.up,v.G,v.v1)
    C._bmm!(v.incident_down,P.R⁺⁻,v.up); v.incident_down .+= P.J₀⁺

    C._jmul!(v.force_up,dL.Ṙ⁻⁺,v.incident_down)
    C._jmul!(v.force_up,dL.Ṫ⁻⁻,v.incident_up,one(eltype(v.force_up)))
    v.force_up .+= dL.J̇₀⁻
    C._jmul!(v.force_down,dL.Ṫ⁺⁺,v.incident_down)
    C._jmul!(v.force_down,dL.Ṙ⁺⁻,v.incident_up,one(eltype(v.force_down)))
    v.force_down .+= dL.J̇₀⁺
    attenuation=reshape(jac.basis.τ̇,1,1,size(jac.basis.τ̇)...)./q.μ₀
    v.below_up .= -S.J₀⁻.*attenuation
    v.below_down .= -S.J₀⁺.*attenuation

    C._jmul!(v.rhs,L.T⁻⁻,v.below_up); v.rhs .+= v.force_up
    C._jmul!(v.t1,v.G,v.rhs)
    C._jmul!(v.t2,P.R⁺⁻,v.t1)
    C._jmul!(v.d,L.T⁺⁺,v.t2); v.d .+= v.force_down
    C._jmul!(v.d,L.R⁺⁻,v.below_up,one(eltype(v.d)))
    C._jmul!(v.t1,v.H,v.d); copyto!(v.d,v.t1)
    C._jmul!(v.t1,S.R⁻⁺,v.d)
    C._jmul!(v.t2,L.T⁻⁻,v.t1); v.t2 .+= v.rhs
    C._jmul!(v.u,v.G,v.t2)
    C._jmul!(w.dc.J̇₀⁻,P.T⁻⁻,v.u)
    C._jmul!(w.dc.J̇₀⁺,S.T⁺⁺,v.d); w.dc.J̇₀⁺ .+= v.below_down
    return nothing
end
