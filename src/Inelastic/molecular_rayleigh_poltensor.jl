#=
Polarizability-tensor backend for `Scattering.molecular_rayleigh_properties(
UsingPolTen(), ...)`.

The public method, its argument checks and its documentation live in
`src/Scattering/molecular_rayleigh.jl`. This file only adapts the existing
Raman-module constants (`getMolecularConstants`) and the existing effective-
coefficient routines (`compute_effective_coefficents!`, `compute_σ_Rayl_coeff!`)
so that the pure-gas API uses exactly the same polarizability tensor as the RT
model: `ConstantAnisotropy()` is the tensor of the Raman cross sections
(`getRamanAtmoConstants`), `DispersiveAnisotropy()` the one of the automatic
(`depol: -1`) Rayleigh depolarization. Each call builds a fresh
`MolecularConstants` object; no shared Raman state is mutated.
=#

_poltensor_molecule(gas::Symbol) =
    gas === :N2 ? N₂() :
    gas === :O2 ? O₂() :
    gas === :H2 ? H₂() :
    throw(ArgumentError("no polarizability tensor for $gas"))

function Scattering._poltensor_rayleigh_coefficients(gas::Symbol,
                                                     ν_cm1::Float64,
                                                     T_K::Float64,
                                                     anisotropy::Scattering.AbstractPolarizabilityDispersion =
                                                         ConstantAnisotropy())
    # The volume mixing ratio does not enter any per-molecule quantity used
    # here; 1.0 denotes the pure gas.
    molecule = _poltensor_molecule(gas)
    mol = getMolecularConstants(molecule, 1.0)
    compute_effective_coefficents!(ν_cm1, T_K, mol, molecule, anisotropy)
    compute_σ_Rayl_coeff!(mol)
    e = mol.effCoeff
    return (α̅ = Float64(e.α̅), γ̅ = Float64(e.γ̅), γ_C = Float64(e.γ_C_Rayl),
            ρ = Float64(e.rho_depol_Rayl), σ_coeff = Float64(e.σ_Rayl_coeff))
end
