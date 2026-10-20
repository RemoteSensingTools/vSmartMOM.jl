using Test
using vSmartMOM
using vSmartMOM.InelasticScattering

# Per-J Raman line strengths obey the Placzek–Teller sum rule: the anisotropic
# scattering γ² is shared out between the O, Q and S branches with weights
# b_{J-2,J} + b_{J,J} + b_{J+2,J} = 1 (Sanghavi 2022, Eqs. 19 and 31), while the
# isotropic part α² stays in the Q branch. Line by line, Q + S + O therefore
# equals the J-independent total α²/3 + (2/27)γ² of Rayleigh scattering
# (elastic, ᾱ and γ̄) or of the vibrational band (α′ and γ′), i.e.
# Rayleigh = Cabannes + RRS (Sanghavi & Frankenberg 2023, Eq. 3).

const _IS = InelasticScattering

# Boltzmann factor exp(-hcE/kT) of level (v, J), with E in cm⁻¹.
_boltzmann(mol, v, J, T) = exp(-_IS.h * _IS.c * mol.effCoeff.E_vJ[v, J] / (_IS.k_B * T))
_spin_weight(mol, J) = mol.gₛ[2 - J % 2]

# Partition function used to normalize every per-J line coefficient.
_partition(mol, T, Jmax) = sum(_spin_weight(mol, J) * (2J + 1) *
    (_boltzmann(mol, 0, J, T) + _boltzmann(mol, 1, J, T)) for J in 0:Jmax)

@testset "Placzek–Teller sum rule of the per-J Raman lines" begin
    Jmax = 30
    for T in (250.0, 296.0)
        n2, o2 = _IS.getRamanAtmoConstants(1e7 / 395.0, T)
        for mol in (n2, o2)
            c = mol.effCoeff
            Z = _partition(mol, T, Jmax)
            elastic = 128π^5 * (c.α̅^2 / 3 + 2c.γ̅^2 / 27)
            vibrational = 128π^5 * (c.α_prime^2 / 3 + 2c.γ_prime^2 / 27)
            for J in 0:(Jmax - 2)
                if _spin_weight(mol, J) == 0
                    @test c.σ_Rayl_coeff_hires[J] == 0
                    continue
                end
                per_state = Z / (_spin_weight(mol, J) * (2J + 1))
                # The elastic Q line counts v = 0, 1, 2; the RRS lines v = 0 only.
                N_Q = sum(_boltzmann(mol, v, J, T) for v in 0:2)
                N_0 = _boltzmann(mol, 0, J, T)
                rotational = per_state * (c.σ_Rayl_coeff_hires[J] / N_Q +
                    (c.σ_RoRaman_coeff_JtoJp2[J] + c.σ_RoRaman_coeff_JtoJm2[J]) / N_0)
                @test rotational ≈ elastic rtol = 1e-10
                vib = per_state * (c.σ_VibRaman_coeff_0to1_hires[J] +
                    c.σ_RoVibRaman_coeff_0to1_JtoJp2[J] +
                    c.σ_RoVibRaman_coeff_0to1_JtoJm2[J]) / N_0
                @test vib ≈ vibrational rtol = 1e-10
                # The 1 → 0 (anti-Stokes) lines start from v = 1.
                N_1 = _boltzmann(mol, 1, J, T)
                vib_10 = per_state * (c.σ_VibRaman_coeff_1to0_hires[J] +
                    c.σ_RoVibRaman_coeff_1to0_JtoJp2[J] +
                    c.σ_RoVibRaman_coeff_1to0_JtoJm2[J]) / N_1
                @test vib_10 ≈ vibrational rtol = 1e-10
            end
        end
    end
end

@testset "Per-J Cabannes lines reproduce the Cabannes depolarization" begin
    # The J-summed elastic Q lines carry α²/3 + (2/27)⟨b_JJ⟩γ²; their depolarization
    # γ_C = 3⟨b⟩γ²/(45α² + 4⟨b⟩γ²) must agree with the energy-conservation
    # route ϖ_Cab = 1 - Σσ_RRS/σ_Rayl (Sanghavi 2026, Eqs. 1–3), up to the
    # (ν̃ - Δν̃)⁴ weighting of the RRS lines that the per-J sum leaves out.
    Jmax = 30
    T = 300.0
    for λ in (395.0, 760.0)
        n2, o2 = _IS.getRamanAtmoConstants(1e7 / λ, T)
        for mol in (n2, o2)
            c = mol.effCoeff
            Z = _partition(mol, T, Jmax)
            w = [_spin_weight(mol, J) * (2J + 1) *
                 sum(_boltzmann(mol, v, J, T) for v in 0:2) / Z for J in 0:Jmax]
            Q = sum(c.σ_Rayl_coeff_hires[J] for J in 0:Jmax) / (128π^5)
            X = (Q - c.α̅^2 / 3 * sum(w)) * 27 / 2 / sum(w)
            γ_lines = 3X / (45c.α̅^2 + 4X)
            γ_cabannes = _IS.compute_γ_mol_Cabannes!(λ, mol)[2]
            @test γ_lines ≈ γ_cabannes rtol = 2e-2
            # The anisotropy fraction is ⟨b_JJ⟩ ≈ 0.26, not ⟨b_JJ²⟩ ≈ 0.07.
            @test 0.24 < X / c.γ̅^2 < 0.28
        end
    end
end

@testset "Q-branch depolarization helper" begin
    # Polarizabilities of order 1e-24 cm³ must not underflow in Float32.
    for FT in (Float32, Float64)
        α, γ = FT(1.7406e-24), FT(0.71e-24)
        for J in 1:6
            b = FT(J * (J + 1) / ((2J - 1) * (2J + 3)))
            γ_C = _IS._q_branch_depolarization(α, γ, b)
            @test γ_C isa FT
            @test isfinite(γ_C)
            @test γ_C ≈ 3b * (γ / α)^2 / (45 + 4b * (γ / α)^2) rtol = 10eps(FT)
        end
        @test iszero(_IS._q_branch_depolarization(α, γ, zero(FT)))
    end
end
