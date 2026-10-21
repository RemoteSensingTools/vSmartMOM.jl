"""
    LocalOpticalJacobian(τ̇, ϖ̇, basis, coefficients)

Factor a layer's optical Jacobian as `dc/dx = (dc/dq) C`, where `basis`
stores the complete optical perturbations `dc/dq` and `C[s,b,p]` maps local
direction b to retrieval column p at wavelength s. `τ̇` and `ϖ̇` retain the
small physical-coordinate derivatives for column attenuation and diagnostics.

Sanghavi et al. (2014), (C.25)–(C.26), supply this chain rule. Each phase
basis direction is an entire matrix, including external-solar columns when
present. It is not an elementwise phase partial after multiple scattering.
"""
struct LocalOpticalJacobian{T,W,B,C} <: AbstractOpticalPropertiesLin
    τ̇::T
    ϖ̇::W
    basis::B
    coefficients::C
end

"Scratch for local elemental/doubling tangents and zero above-layer seeds."
struct LocalJacobianWorkspace{L,A}
    layer::L
    above::A
end
