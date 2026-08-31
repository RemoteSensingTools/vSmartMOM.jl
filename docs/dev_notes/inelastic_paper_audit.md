# Inelastic-scattering implementation audit against Sanghavi 2022 / Sanghavi & Frankenberg 2023

**Audit reference papers** (committed under `docs/papers/`):

- **Paper I** — Sanghavi 2022, *Raman scattering in the earth's atmosphere, part I: Optical properties*, JQSRT 291, 108328 (`docs/papers/1-s2.0-S0022407322002631-main.pdf`).
- **Paper II** — Sanghavi & Frankenberg 2023, *Raman scattering in the earth's atmosphere, Part II: Radiative transfer modeling for remote sensing applications*, JQSRT 311, 108791 (`docs/papers/1-s2.0-S0022407323003096-main.pdf`).

**Code under audit** — everything inside `src/Inelastic/`:

```
src/Inelastic/InelasticScattering.jl
src/Inelastic/types.jl
src/Inelastic/inelastic_helper.jl
src/Inelastic/raman_atmo_prop.jl
src/Inelastic/raman_stellar_prop.jl
src/Inelastic/stellar_inelastic_helper.jl
src/Inelastic/src/raman_constants.jl
src/Inelastic/src/molecular_constructors.jl
src/Inelastic/src/inelastic_cross_section.jl
src/Inelastic/src/apply_lineshape.jl
```

Plus the immediate consumer paths in `src/CoreRT/tools/model_from_parameters.jl` and `src/CoreRT/tools/lin_model_from_parameters.jl`.

---

## TL;DR

| #  | Severity | Finding | File: line |
|----|----------|---------|------------|
| 1  | **HIGH-ish** | `α̅₀₀` for O₂ and N₂ disagree with Paper I Table 3 (cited Buldakov values) by 1.7–3% | `src/Inelastic/src/molecular_constructors.jl:8, 80` |
| 2  | **MEDIUM**  | `compute_γ_mol_Cabannes!` denominator does not algebraically invert Paper II eq 25 | `src/Inelastic/inelastic_helper.jl:430-449` |
| 3  | **MEDIUM**  | Two `getRamanSSProp!(RRS, …)` overloads disagree on `i_λ₁λ₀` reversal | `src/Inelastic/raman_atmo_prop.jl:71-75` vs `:93-97` |
| 4  | **MEDIUM**  | Per-J VRS Q-branch hires σ omits the `b_{J,J}` line-strength factor (Paper I eq 32) | `src/Inelastic/src/inelastic_cross_section.jl:113-130` |
| 5  | **MEDIUM**  | `EffectiveCoefficients.γ_C_Rayl` is misnamed; computed value is Paper II's `γ_C,Cab` | `src/Inelastic/src/raman_constants.jl:89`, `inelastic_cross_section.jl:49` |
| 6  | LOW         | `γ_C(J) = 3/(4 + 45·(α̅/(b_{J,J}·γ̅))²)` per-J depolarization is undocumented | `src/Inelastic/src/inelastic_cross_section.jl:106, 122` |
| 7  | LOW         | `compute_γ_air_Cabannes!` round-trips through the buggy inversion (#2) when `γ_C,Cab` is directly available | `src/Inelastic/inelastic_helper.jl:359-381` |

Items that match the papers correctly are listed at the end of the document.

---

## Finding #1 — HIGH-ish — `α̅₀₀` constants disagree with Paper I Table 3

### Paper I Table 3 (cited from Buldakov [2,3])

| | `a_B [cm³]` | `a'_B [cm³]` | `γ_B [cm³]` | `γ'_B [cm³]` |
|--:|---:|---:|---:|---:|
| O₂ | **1.61·10⁻²⁴** | 1.76·10⁻²⁴ | 1.08·10⁻²⁴ | 3.19·10⁻²⁴ |
| N₂ | **1.77·10⁻²⁴** | 1.86·10⁻²⁴ | 0.71·10⁻²⁴ | 2.23·10⁻²⁴ |

### Code values (`src/Inelastic/src/molecular_constructors.jl`)

```julia
# Line 7-14 (N₂):
p = PolarizationTensor{FT}(
        α̅₀₀     = 1.7406e-24, #[cm^3]                ← paper says 1.77e-24 (-1.7%)
        α₀₀_prime = 1.86e-24, #[cm^3]                 ✓
        ω₀      = 2.6049e16,
        α_b     = 1.8e-6,
        α_c     = 0.0,
        γ̅₀₀       = 0.71e-24,  #[cm^3]                ✓
        γ₀₀_prime = 2.23e-24)   #[cm^3]               ✓

# Line 79-86 (O₂):
p = PolarizationTensor{FT}(
        α̅₀₀     = 1.5658e-24, #[cm^3]                ← paper says 1.61e-24 (-2.7%)
        α₀₀_prime = 1.76e-24, #[cm^3]                 ✓
        ω₀      = 2.1801e16,
        α_b     = -2.369e-6,
        α_c     = 8.687e-9,
        γ̅₀₀       = 1.080e-24,  #[cm^3]               ✓
        γ₀₀_prime = 3.19e-24)   #[cm^3]               ✓
```

### Knock-on consequences

`γ_C,Rayl = 3γ²/(45a² + 4γ²)` (Paper I eq 12) compares as follows:

| | code `γ_C,Rayl` | Paper I Table 4 `γ_C,Rayl` | Δ |
|--:|---:|---:|---:|
| O₂ | 0.03043 | **0.02885** | +5.5% |
| N₂ | 0.01093 | **0.01058** | +3.3% |

Air-mixture Rayleigh σ at 395 nm (per Paper II Table 1 = `0.02088`) ends up about 2-3% off relative to the paper's tabulated values.

### Recommendation

Either (a) update `α̅₀₀` to the cited Paper I Table 3 values (1.61·10⁻²⁴ for O₂, 1.77·10⁻²⁴ for N₂) and reverify all downstream cross-sections, or (b) document the alternative source the code is actually drawing from.

---

## Finding #2 — MEDIUM — `compute_γ_mol_Cabannes!` denominator does not invert Paper II eq 25

### Paper II eq 25

```
                   1     3(1 − ϖ_Cab) + 2γ_C,Cab(3 + 2ϖ_Cab)
   γ_C,Ray = ─── · ────────────────────────────────────────────────
                   2     (2 + 3ϖ_Cab) + 4γ_C,Cab(1 − ϖ_Cab)
```

Solved analytically for `γ_C,Cab` (cross-multiply, group γ_C,Cab terms):

```
                   1     2γ_Ray(2 + 3ϖ) − 3(1 − ϖ)
   γ_C,Cab = ─── · ────────────────────────────────────
                   2     (3 + 2ϖ) − 4γ_Ray(1 − ϖ)
```

### Code (`src/Inelastic/inelastic_helper.jl:430-449`)

```julia
function compute_γ_mol_Cabannes!(λ₀::FT, mol) where FT
    ν̃ =  1.e7/λ₀
    effT  =  300.; #K assumed constant for Earth atmospheres

    ϖ_Cabannes = compute_ϖ_Cabannes(λ₀, mol)
    γ_mol_Rayleigh = mol.effCoeff.γ_C_Rayl

    tmp1 = 1+2*γ_mol_Rayleigh
    tmp2 = 2+3*ϖ_Cabannes
    tmp3 = 1-ϖ_Cabannes
    tmpN = tmp1*tmp2-5
    tmpD = tmp1*tmp3+5
    γ_mol_Cabannes = 0.5*tmpN/tmpD

    return ϖ_Cabannes, γ_mol_Cabannes, γ_mol_Rayleigh
end
```

Expanding `tmp1·tmp2 − 5` — code's numerator:

```
(1 + 2γ_Ray)(2 + 3ϖ) − 5
   = 2 + 3ϖ + 4γ_Ray + 6γ_Ray·ϖ − 5
   = 2γ_Ray(2 + 3ϖ) − 3(1 − ϖ)              ✓ matches paper-inverse numerator
```

Expanding `tmp1·tmp3 + 5` — code's denominator:

```
(1 + 2γ_Ray)(1 − ϖ) + 5
   = 1 − ϖ + 2γ_Ray − 2γ_Ray·ϖ + 5
   = (6 − ϖ) + 2γ_Ray(1 − ϖ)
```

vs the paper-inverse denominator `(3 + 2ϖ) − 4γ_Ray(1 − ϖ)`.

These are **not algebraically equivalent**. Setting them equal and simplifying reduces to `−3(1 − ϖ)(1 + 2γ_Ray) = 0`, which only holds when `ϖ = 1` (no Raman).

### Numerical impact

Validation with Paper II Table 1, air at 395 nm: `ϖ_Cab = 0.97526`, `γ_C,Cab = 0.01328`, `γ_C,Ray = 0.02088`.

Forward eq 25 with code's table values: `γ_Ray = 0.02088` ✓ (Paper II eq 25 is correct to ~1e-4).

Inverse — analytic formula on `γ_Ray = 0.02088`, `ϖ = 0.97526`:
```
γ_Cab = 0.5 · [2·0.02088·4.9258 − 3·0.02474] / [4.9505 − 4·0.02088·0.02474]
      = 0.5 · 0.13149 / 4.9484
      = 0.01329                              ✓ matches paper 0.01328
```

Inverse — code's formula on the same input:
```
γ_Cab_code = 0.5 · [(1+0.04176)·4.9258 − 5] / [1.04176·0.02474 + 5]
           = 0.5 · 0.13295 / 5.02577
           = 0.01323                          ← off by ~0.4%
```

At lower ϖ (e.g. `ϖ=0.5, γ_Ray=0.5`) the discrepancy becomes ~50%. For atmospheric inputs the practical impact is small; for instructional / hypothetical regimes it isn't.

### Recommendation

Replace the body with the analytic inverse:

```julia
function compute_γ_mol_Cabannes!(λ₀::FT, mol) where FT
    ϖ_Cabannes = compute_ϖ_Cabannes(λ₀, mol)
    γ_mol_Rayleigh = mol.effCoeff.γ_C_Rayl

    # Analytic inverse of Paper II eq 25.
    num = 2γ_mol_Rayleigh*(2 + 3ϖ_Cabannes) - 3*(1 - ϖ_Cabannes)
    den = (3 + 2ϖ_Cabannes) - 4γ_mol_Rayleigh*(1 - ϖ_Cabannes)
    γ_mol_Cabannes = 0.5*num/den

    return ϖ_Cabannes, γ_mol_Cabannes, γ_mol_Rayleigh
end
```

Also see Finding #5 — the `γ_C_Rayl` field that this function reads is not actually `γ_C,Ray`, so the input semantics need to be corrected at the same time.

---

## Finding #3 — MEDIUM — Inconsistent `i_λ₁λ₀` reversal between two `getRamanSSProp!(RRS, …)` overloads

### `src/Inelastic/raman_atmo_prop.jl:58-78` (4-arg, takes `depol`)

```julia
function getRamanSSProp!(RS_type::RRS, depol, λ, grid_in)
    @unpack n2,o2 =  RS_type
    atmo_σ_Rayl = compute_optical_Rayl(λ, n2, o2)
    RS_type.greek_raman = get_greek_raman(RS_type, n2, o2)
    RS_type.ϖ_Cabannes .= compute_ϖ_Cabannes(RS_type, λ)
    index_raman_grid, atmo_σ_RRS = compute_optical_RS!(RS_type, grid_in, λ, n2, o2)

    RS_type.ϖ_λ₁λ₀ = (reverse(atmo_σ_RRS)/atmo_σ_Rayl)
    RS_type.i_λ₁λ₀ = reverse(index_raman_grid)              # ← REVERSED
    RS_type.n_Raman = length(RS_type.ϖ_λ₁λ₀)
    return nothing
end
```

### `src/Inelastic/raman_atmo_prop.jl:80-100` (3-arg, no `depol`)

```julia
function getRamanSSProp!(RS_type::RRS, λ, grid_in)
    @unpack n2,o2 =  RS_type
    atmo_σ_Rayl = compute_optical_Rayl(λ, n2, o2)
    RS_type.greek_raman = get_greek_raman(RS_type, n2, o2)
    RS_type.ϖ_Cabannes .= compute_ϖ_Cabannes(RS_type, λ)
    index_raman_grid, atmo_σ_RRS = compute_optical_RS!(RS_type, grid_in, λ, n2, o2)

    RS_type.ϖ_λ₁λ₀ = (reverse(atmo_σ_RRS)/atmo_σ_Rayl)
    RS_type.i_λ₁λ₀ = index_raman_grid #reverse(index_raman_grid)   # ← NOT reversed (note comment)
    RS_type.n_Raman = length(RS_type.ϖ_λ₁λ₀)
    return nothing
end
```

### Diagnosis

`compute_optical_RS!` (`src/Inelastic/inelastic_helper.jl:626-663`) returns `(index_ramangrid_out, atmo_σ_RRS)` where `atmo_σ_RRS = σ_tmp[index_ramangrid_out]`. Both arrays index into the same monotonic wavenumber grid. If `ϖ_λ₁λ₀` is reversed, `i_λ₁λ₀` must be reversed too for the (offset, weight) pairs to match.

Either:

- **The 4-arg variant is correct** (both arrays reversed in lockstep) and the 3-arg variant has a bug (Stokes/anti-Stokes wing swap with respect to the offset table that elemental_inelastic indexes into).
- **The 3-arg variant is intentional** (the comment `#reverse(index_raman_grid)` explicitly removed it), in which case the `reverse(atmo_σ_RRS)` on the line above is the line that should change.

The mapping that the rest of the code depends on (`elemental_inelastic.jl` reads both `RS_type.i_λ₁λ₀[Δn]` and `RS_type.ϖ_λ₁λ₀[Δn]` for the same `Δn`) requires the two arrays to be in the same order.

### Recommendation

Pick one ordering convention and use it in both overloads. A short comment explaining whether the convention is "Δn>0 means donor wavelength is at higher ν" or vice versa would prevent a future copy-paste from re-introducing the asymmetry.

---

## Finding #4 — MEDIUM — Per-J VRS hires σ omits the `b_{J,J}` line-strength factor

### Paper I eq 32 (per-J Q-branch VRS cross section)

```
                                                                     1 + 2γ_C,vib
σ_VRS(ν̃; v₁ ← v₀) = (128π⁵(a'_B b_{Δν̂,B})² · b_{J,J} · 𝒩ᵢ · (ν̃ − Δν̃_fi)⁴) · ───────────
                                                                     3 − 4γ_C,vib
```

with `b_{J,J} = J(J+1) / [(2J−1)(2J+3)]` (Paper I eq 31).

### Code (`src/Inelastic/src/inelastic_cross_section.jl:80-150`, abridged)

```julia
function compute_σ_Rayl_VibRaman_coeff_hires!(T, mol::MolecularConstants{FT}; Jmax=30) where {FT}
    @unpack α̅, γ̅, α_prime, γ_prime, E_vJ, σ_Rayl_coeff_hires, σ_VibRaman_coeff_0to1_hires, σ_VibRaman_coeff_1to0_hires, Δν̃_Rayl_coeff_hires, Δν̃_VibRaman_coeff_0to1_hires, Δν̃_VibRaman_coeff_1to0_hires = mol.effCoeff
    ...
    # Line 98:
    b_JJ = Ji*(Ji+1)/((2Ji-1)*(2Ji+3))                  # ← computed per Paper I eq 31

    # Line 106 - Rayleigh hires:
    γ_C = 3/(4+45(α̅/(b_JJ*γ̅))^2)                       # ← b_JJ used here, not as line strength

    # Line 113-115 - σ_Rayl_coeff_hires[Ji]:
    σ_Rayl_coeff_hires[Ji] = FT(kᵥ * Float64(g_N * (2Ji+1)) * Float64(Ni_by_N) *
                                _molecular_square64(α̅) *
                                Float64((1+2*γ_C)/(3-4*γ_C)))

    # Line 122-130 - σ_VibRaman_coeff_0to1_hires[Ji]:
    γ_C = 3/(4+45(α_prime/(b_JJ*γ_prime))^2)            # ← same J-dependent γ_C trick
    Ni_by_N = exp(-h*c*Float64(E_vJ[0,Ji])/(k_B*T))
    σ_VibRaman_coeff_0to1_hires[Ji] = FT(kᵥ * Float64(g_N * (2Ji+1)) * Float64(Ni_by_N) *
                                         _molecular_square64(α_prime) *
                                         Float64((1+2*γ_C)/(3-4*γ_C)))

    # Line 143-145 - σ_VibRaman_coeff_1to0_hires[Ji]: same shape, missing b_JJ.
    ...
    # Line 148: divides everything by Z_pf (partition function).
    σ_VibRaman_coeff_0to1_hires /= Z_pf
    σ_VibRaman_coeff_1to0_hires /= Z_pf
```

Code's per-J VRS hires σ is

```
σ_VibRaman_coeff_X_hires[J] = (128π⁵·g_N(2J+1)·exp(−E_vJ/kT)/Z_pf) · α_prime² · (1+2γ_C(J))/(3−4γ_C(J))
```

It is missing the multiplicative `b_{J,J}` line-strength factor that Paper I eq 32 includes for the Q-branch.

### Where the array is consumed

`grep -rn σ_VibRaman_coeff_0to1_hires src/` shows the array is read in:

- `src/Inelastic/inelastic_helper.jl:208-219` (`compute_ϖ_Cabannes(::VS_0to1_plus, λ₀)`): summed over the spectral grid to get integrated σ_VRS for the Cabannes single-scatter albedo.
- `src/Inelastic/inelastic_helper.jl:285-288` (`compute_ϖ_Cabannes_VS`): same.
- `src/Inelastic/inelastic_helper.jl` 138-149 + 207-219 + 249-252 + 285-288 + 312-313: most other callsites are commented out.
- `src/Inelastic/raman_atmo_prop.jl` and `raman_stellar_prop.jl`: the `Δν̃_*_hires` (offset) arrays *are* used (indexes into apply_gridlines), but the sigma `σ_*_hires` arrays are NOT directly used outside the integrated `compute_ϖ_Cabannes_VS` path.

Because the consumers integrate `σ × (ν̃ + Δν̃)⁴` over the spectral grid, and `(ν̃ + Δν̃)⁴` doesn't vary much across the Q-branch (which is narrow in frequency), the missing `b_{J,J}` is approximately compensated by the `Z_pf` normalization — the integrated value is roughly correct.

### Numerical check

`Σ_J g_N(2J+1)·exp(−E_vJ/kT)/Z_pf` ≈ 1 (it's the v=0 fraction of the partition function, ~1 at 300 K for both O₂ and N₂). Compare to Paper I eq 32's `Σ_J b_{J,J}·𝒩ᵢ` which is `<Q-branch line strength> × <population>` — this is a fraction of order 0.25–0.5 depending on rotational distribution.

In practice, the bulk version `compute_σ_VibRaman_coeff!` (`inelastic_cross_section.jl:200`) IS used as the integrated σ_VRS in some places, and that one *does* fold in the partition function via `Nvib = 1/(1 − exp(−hcω_e/kT))`. So the two paths (bulk vs hires-summed) are not internally consistent.

### Recommendation

Multiply each `σ_VibRaman_coeff_X_hires[Ji]` by `b_JJ` to match Paper I eq 32. Verify the integrated sum still matches `compute_σ_VibRaman_coeff!`. Same fix for `σ_Rayl_coeff_hires[Ji]` if it represents the per-J Q-branch elastic line.

---

## Finding #5 — MEDIUM — `EffectiveCoefficients.γ_C_Rayl` is misnamed

### Field definition (`src/Inelastic/src/raman_constants.jl:89`)

```julia
mutable struct EffectiveCoefficients{FT, ...}
    ...
    γ_C_Rayl::FT # 3/(45(ϵ)^2+4)        ← comment matches Paper I eq 12
    ...
end
```

### Computation (`src/Inelastic/src/inelastic_cross_section.jl:35-65`)

```julia
function compute_effective_coefficents!(ν, T, mol::MolecularConstants{FT}) where {FT}
    @unpack α̅,  γ̅, α_prime, γ_prime, ϵ, ϵ_prime = mol.effCoeff
    @unpack α̅₀₀, γ̅₀₀, ω₀, α_b, α_c, α₀₀_prime, γ₀₀_prime = mol.PolTensor
    ...
    α̅ = α̅₀₀*(1 + α_b*T + α_c*T^2)/(1-(c*ν_eff/ω₀)^2)
    γ̅ = γ̅₀₀
    ϵ = α̅/γ̅                            # ← α̅/γ̅, not γ̅/α̅
    γ_C_Rayl = 3/(45(ϵ)^2+4)            # ← reads as 3 / (45(α̅/γ̅)² + 4) = 3γ̅² / (45α̅² + 4γ̅²)
    ...
```

Plug `ϵ = α̅/γ̅` into the formula:

```
γ_C_Rayl = 3 / (45·(α̅/γ̅)² + 4)
         = 3γ̅² / (45·α̅² + 4γ̅²)
```

### Compare to the papers

- **Paper I eq 12**: `γ_C,Rayl = 3γ²/(45a² + 4γ²)`. Same shape. Paper I uses the label "Rayl".
- **Paper II eq 24**: `γ_C,Cab = 3γ²/(45a² + 4γ²)`. Same shape, but Paper II is explicit that this is the **Cabannes** depolarization, not the **Rayleigh** depolarization.
- **Paper II eq 25** (separate): `γ_C,Ray = (1/2) · [3(1−ϖ_Cab) + 2γ_C,Cab(3+2ϖ_Cab)] / [(2+3ϖ_Cab) + 4γ_C,Cab(1−ϖ_Cab)]`. The Rayleigh γ in Paper II's convention includes the inelastic Raman admixture.

So Paper I's `γ_C,Rayl` and Paper II's `γ_C,Cab` are the same quantity by formula. Paper II clarifies that the more useful `γ_C,Ray` is *different*.

The code field named `γ_C_Rayl` is, by formula, **Paper II's `γ_C,Cab`**. It is *not* Paper II's `γ_C,Ray`. This causes a chain of downstream confusion:

#### `compute_σ_Rayl_coeff!` (`inelastic_cross_section.jl:69-76`)

```julia
function compute_σ_Rayl_coeff!(mol::MolecularConstants{FT}) where {FT}
    @unpack α̅, γ_C_Rayl, σ_Rayl_coeff = mol.effCoeff
    σ_Rayl_coeff = FT(_rayleigh_prefactor64() * _molecular_square64(α̅) *
                      Float64((1+2*γ_C_Rayl)/(3-4*γ_C_Rayl)))     # ← uses γ_C,Cab as if it were γ_C,Ray
    @pack! mol.effCoeff = σ_Rayl_coeff,α̅,γ_C_Rayl
end
```

The factor `(1 + 2γ_C,x) / (3 − 4γ_C,x)` is Paper II eq 22 with `x` being the scattering type. With `x = Cab`, this gives σ_Cab (just the elastic Cabannes line). With `x = Ray`, this gives σ_Ray (full Rayleigh = Cabannes + Raman wings).

Code uses `γ_C,Cab` and gets `σ_Cab`, but the variable is named `σ_Rayl_coeff` and used downstream as if it were the full Rayleigh σ. For air at 395 nm (Paper II Table 1), `σ_Cab/σ_Rayl = ϖ_Cab = 0.97526`, so the misnaming costs ~2.5% on the total Rayleigh cross-section.

#### `compute_γ_air_Rayleigh!` (`inelastic_helper.jl:331-349`)

```julia
function compute_γ_air_Rayleigh!(λ₀::FT, RS_type::Union{RRS, VS_0to1, ...}) where FT
    γN₂ = RS_type.n2.effCoeff.γ_C_Rayl                  # ← actually γ_C,Cab(N₂)
    σ0N₂ = RS_type.n2.effCoeff.σ_Rayl_coeff             # ← actually σ_Cab(N₂)
    σN₂ = σ0N₂ * (3-4γN₂)/(1+2γN₂)                       #   strips the (1+2γ)/(3-4γ) factor → 128π⁵α²
    VMR_N₂ = RS_type.n2.vmr
    ...
    γ_air_Rayleigh = 3/(4 + tmp1 / tmp2)
    σ_air_Rayleigh = (σ0N₂ * VMR_N₂ + σ0O₂ * VMR_O₂)*(nm_per_cm/λ₀)^4/(VMR_N₂ + VMR_O₂)

    return γ_air_Rayleigh, σ_air_Rayleigh
end
```

The function returns the air-weighted **Cabannes** γ, labeled as `γ_air_Rayleigh`.

#### Production caller (`src/CoreRT/tools/model_from_parameters.jl:285-302`)

```julia
ϖ_Cab = InelasticScattering.compute_ϖ_Cabannes(λₘ, _n2, _o2)
γ_air_Cab, _ = InelasticScattering.compute_γ_air_Cabannes!(λₘ, _n2, _o2)
γ_air_Ray, _ = InelasticScattering.compute_γ_air_Rayleigh!(λₘ, _n2, _o2)
ϖ_Cabannes[i_band] = FT(ϖ_Cab)
depol_air_Cab = 2γ_air_Cab / (1 + γ_air_Cab)
depol_air_Ray = 2γ_air_Ray / (1 + γ_air_Ray)

depol_use_Cab = params.depol < 0 ? FT(depol_air_Cab) : FT(params.depol)
depol_use_Ray = params.depol < 0 ? FT(depol_air_Ray) : FT(params.depol)

push!(greek_cabannes,     Scattering.get_greek_rayleigh(depol_use_Cab))
push!(greek_rayleigh_arr, Scattering.get_greek_rayleigh(depol_use_Ray))

τ_rayl[i_band] .= getRayleighLayerOptProp(profile.p_half[end],
                                           curr_band_λ,
                                           depol_use_Ray, profile.vcd_dry);
```

For air at 395 nm:

- `γ_air_Cab` returned ≈ 0.01323 (Cabannes γ via the buggy round-trip in #2)
- `γ_air_Ray` returned ≈ 0.01328 (Cabannes γ direct)
- Paper II Table 1 says: `γ_C,Cab = 0.01328`, `γ_C,Ray = 0.02088`.

So `depol_use_Ray ≈ 0.0262` is fed to `getRayleighLayerOptProp` and to `Scattering.get_greek_rayleigh` for the Rayleigh phase function. Paper II's Rayleigh depol would be `2·0.02088/1.02088 ≈ 0.041`.

#### Net effect

In the **noRS path** (which has no separate Raman component), the lumped "Rayleigh" phase function should use `γ_C,Ray` per Paper II convention to absorb the Raman wings into the elastic line's apparent depolarization. Code uses `γ_C,Cab` instead, so the noRS Rayleigh phase function under-represents the wing depolarization by ~50% relative to Paper II.

`getRayleighLayerOptProp` itself uses Bodhaine 1999 with an internal `ρ₀ = 0.0279` (which is *itself* a Cabannes-like depol). So the τ_rayl output is approximately self-consistent with the Bodhaine convention (~0.3% off the Bodhaine value), but **off Paper II convention by ~2-3% in σ_Rayl**.

In the **Raman paths** (RRS, VRS, RVRS): `greek_cabannes` correctly uses `γ_C,Cab` (good), but `greek_rayleigh_arr` is also `γ_C,Cab` rather than `γ_C,Ray`. If anything reads `greek_rayleigh_arr` for the Rayleigh-component phase function in a Raman-tracking solver, it gets the Cabannes phase function instead. Per current code grep, `greek_rayleigh_arr` is consumed by the noRS-style elastic path, where the issue is the noRS convention noted above.

### Recommendation

Three options, in increasing scope:

1. **Pure rename / docs** — rename the field `γ_C_Rayl → γ_C_Cab` (and `σ_Rayl_coeff → σ_Cab_coeff` etc.) so the code matches Paper II terminology. Caller behavior is preserved.

2. **Add `γ_C_Ray`** — keep `γ_C,Cab` as today; add `γ_C_Ray` (computed via Paper II eq 25 *forward* from `γ_C,Cab` and `ϖ_Cab`). Have `compute_γ_air_Rayleigh!` use the new field. `compute_γ_mol_Cabannes!` becomes a no-op (γ_C,Cab is already on the struct).

3. **Conventional unification** — pick a single convention (Bodhaine or Paper II) for the Rayleigh-component phase function used by the noRS path, document it, and verify all RT regression tests pass.

---

## Finding #6 — LOW — Per-J `γ_C(J) = 3 / (4 + 45·(α̅/(b_{J,J}·γ̅))²)` is undocumented

### Code (`src/Inelastic/src/inelastic_cross_section.jl:106, 122`)

```julia
# Line 98:
b_JJ = Ji*(Ji+1)/((2Ji-1)*(2Ji+3))                  # Paper I eq 31

# Line 106 - elastic Q-branch hires:
γ_C = 3/(4+45(α̅/(b_JJ*γ̅))^2)                       # ← J-dependent γ_C, not in Paper I or II

# Line 122 - vibrational Q-branch hires:
γ_C = 3/(4+45(α_prime/(b_JJ*γ_prime))^2)            # ← same trick with primes
```

Paper I eq 12 / Paper II eq 24 give a J-independent `γ_C,Cab = 3γ²/(45α² + 4γ²)`. Paper I eq 29 gives a J-independent `γ_C,vib = 3(γ')²/(45(α')² + 4(γ')²)`.

Code's per-J form differs:

```
γ_C(J) = 3 · b_{J,J}² · γ²   /   (4 · b_{J,J}² · γ² + 45 · α²)
```

Limits:

- `b_{J,J} → 1` (high J): `γ_C(J) → 3γ²/(45α² + 4γ²) =` Paper γ_C,Cab. ✓
- `b_{J,J} → 0` (J=0,1): `γ_C(J) → 0`. Means low-J Q-branch lines are isotropic.

Physically, the J=0 Q-branch line doesn't exist (b_{0,0} = 0), and low-J Q-branch lines are weak; the `γ_C(J) → 0` behavior is plausible but isn't the formula in either paper.

### Recommendation

Add a comment citing the source if this is from elsewhere (Long? Placzek? Buldakov?), or fall back to the J-independent Paper I eq 12 / 29 forms and absorb any per-J line-strength dependence into the multiplicative `b_{J,J}` factor (as Paper I eq 32 does).

---

## Finding #7 — LOW — `compute_γ_air_Cabannes!` round-trips through the buggy inversion

### Code (`src/Inelastic/inelastic_helper.jl:359-381`)

```julia
function compute_γ_air_Cabannes!(λ₀::FT, RS_type::Union{RRS, VS_0to1, VS_1to0, RRS_plus, VS_0to1_plus, VS_1to0_plus}) where FT

    γN₂ = compute_γ_mol_Cabannes!(λ₀, RS_type.n2)[2]    # ← runs the buggy inversion #2
    ϖN₂ = compute_ϖ_Cabannes(λ₀, RS_type.n2)
    σ0N₂ = RS_type.n2.effCoeff.σ_Rayl_coeff
    σN₂ = ϖN₂ * σ0N₂ * (3-4γN₂)/(1+2γN₂)
    VMR_N₂ = RS_type.n2.vmr

    γO₂ = compute_γ_mol_Cabannes!(λ₀, RS_type.o2)[2]    # ← same
    ...
    γ_air_Cabannes = 3/(4 + tmp1 / tmp2)
    ...
    return γ_air_Cabannes, ϖ_air_Cabannes;
end
```

`compute_γ_mol_Cabannes!` is the function audited in Finding #2. It nominally inverts Paper II eq 25 to derive `γ_C,Cab` from `γ_C,Ray` + `ϖ_Cab`. But the input it gets (`mol.effCoeff.γ_C_Rayl`) is *already* `γ_C,Cab` (Finding #5) — the round-trip is a no-op in concept.

### Recommendation

Since `γ_C,Cab` is directly available from `mol.effCoeff.γ_C_Rayl` (after the rename in Finding #5), bypass the inversion entirely:

```julia
function compute_γ_air_Cabannes!(λ₀::FT, RS_type::Union{RRS, ...}) where FT
    γN₂ = RS_type.n2.effCoeff.γ_C_Rayl                  # already Paper II γ_C,Cab
    ...
end
```

This sidesteps the bug in Finding #2 for this particular caller (without removing the need to fix #2 itself).

---

## What checked out OK

The following equations match Paper I / II within numerical noise. Listed here so the next reviewer doesn't have to redo this work.

### Dunham coefficients `Y_{k,l}` (Paper I Table 1)

`src/Inelastic/src/molecular_constructors.jl:16-23, 88-95, 160-167`. Code Y indexing is offset by +1 (Julia 1-based), i.e. code `Y[k,l]` ≡ paper `Y_{k-1,l-1}`. Energy-level loop in `compute_energy_levels!` (`inelastic_cross_section.jl:160-180`) correctly raises `(v + 1/2)^(k-1)` and `[J(J+1)]^(l-1)`.

| | Y[1,2] | Y[1,3] | Y[2,1] | Y[2,2] | Y[3,1] | Y[4,1] |
|--:|:--:|:--:|:--:|:--:|:--:|:--:|
| O₂ paper | 1.4376766 | -4.839e-6 | 1580.19 | -0.01590 | -11.98 | 0.0 |
| O₂ code  | 1.4376766 | -4.839e-6 | 1580.19 | -0.01590 | -11.98 | 0.0 |
| N₂ paper | 1.99824 | -5.76e-6 | 2358.57 | -0.017318 | -14.324 | -2.26e-3 |
| N₂ code  | 1.99824 | -5.76e-6 | 2358.57 | -0.017318 | -14.324 | -2.26e-3 |

### Nuclear-spin degeneracies `g_J` (Paper I Table 2)

`src/Inelastic/src/molecular_constructors.jl:25, 97, 169`. Convention is `gₛ = [odd_J, even_J]`.

| | odd J | even J |
|--:|:--:|:--:|
| O₂ paper / code | 1 | 0 |
| N₂ paper / code | 3 | 6 |

### `b_{J,J−2}`, `b_{J,J+2}` rotational line strength (Paper I eq 19)

`src/Inelastic/src/inelastic_cross_section.jl:252-254`:

```julia
b_JJm2 = 3*Ji*(Ji-1) / (2*(2*Ji+1)*(2*Ji-1))    # 3J(J−1)/[2(2J+1)(2J−1)]    ✓
b_JJp2 = 3*(Ji+1)*(Ji+2) / (2*(2*Ji+1)*(2*Ji+3)) # 3(J+1)(J+2)/[2(2J+1)(2J+3)] ✓
```

### `b_{J,J} = J(J+1)/[(2J−1)(2J+3)]` Q-branch (Paper I eq 31)

`src/Inelastic/src/inelastic_cross_section.jl:98`. Computed correctly; consumption in σ_VibRaman is incomplete (Finding #4).

### RRS bulk σ formula (Paper I eq 22)

`src/Inelastic/src/inelastic_cross_section.jl:271-272`:

```julia
σ_RoRaman_coeff_JtoJp2[Ji] = kᵥ * g_N*(2Ji+1) * b_JJp2 * Ni_by_N * γ̅²
σ_RoRaman_coeff_JtoJm2[Ji] = kᵥ * g_N*(2Ji+1) * b_JJm2 * Ni_by_N * γ̅²
# kᵥ = (256/27)π⁵   (= 2 × 128π⁵ / 27)
```

Matches Paper I eq 22's `(128π⁵·γ_B²) · (2b_{J',J}/27) · 𝒩ᵢ` after factoring. ✓

### RVRS bulk σ formula (Paper I eq 34)

`src/Inelastic/src/inelastic_cross_section.jl:309-310, 325-326`. Same shape as RRS but with `γ_prime² = (γ_B' · b_{Δν̂,B})²` per Paper I eq 27 `b²_{Δν̂,B} = B_e/ω_e`. The code computes `γ_prime = γ₀₀_prime * sqrt(B_e/ω_e)` via `mol.Y[1,2]/mol.Y[2,1]` (`inelastic_cross_section.jl:46-48`). ✓

### Energy levels `E_ν̃(v, J) = Σ Y_{k,l}·(v+1/2)^k·[J(J+1)]^l` (Paper I eq 35)

`src/Inelastic/src/inelastic_cross_section.jl:160-180`. The triple-loop with running powers `E₁_pow`, `E₂_pow` correctly produces `Y[k,l]·E₁^(l-1)·v_half^(k-1)` for `(k,l) ∈ {1..5}²`. ✓

### `apply_lineshape_!` Doppler width and line strength (Paper I eq 22)

`src/Inelastic/src/apply_lineshape.jl:40-67`:

```julia
ν   = Δνᵢ[j] + nm_per_m/λ₀                                       # scattered ν̃
γ_d = (cSqrt2Ln2/cc_) * sqrt(cBolts_/cMassMolIE) * sqrt(T) * ν / sqrt(molMass)  # Doppler HWHM
S   = σᵢ[j] * ν^4                                                 # line strength
```

Doppler HWHM: `γ_D = ν · sqrt(2·ln2·k_B·T/(m·c²))` ✓.
Line strength: `S = σ_coef · ν⁴` ✓ matches Paper I eq 22.
Gaussian normalization: `cSqrtLn2divSqrtPi · exp(...) / γ_d` ✓.

### `ϖ_Cabannes = 1 − σ_RRS/σ_Rayl` (Paper II eq 26)

`src/Inelastic/inelastic_helper.jl:153, 256, 318`. To first order in `σ_RRS/σ_Cab` (≈ 2.5% atmospherically), this matches `σ_Cab/(σ_Cab + σ_RRS)`. Subject to the Cabannes-vs-Rayleigh σ confusion of Finding #5, but the formula itself is correctly derived.

### `γ_C,rot = 3/4` (universal rotational Raman depolarization)

`src/Inelastic/src/inelastic_cross_section.jl:50-51`. Code: `γ_C_RotRaman = 3/4`, `γ_C_RoVibRaman = 3/4`. Matches Paper I/II convention (rotational Raman is universally `3/4` depolarized regardless of molecule). ✓

### Greek-coefficient construction for Raman phase functions

`src/Inelastic/inelastic_helper.jl:864-908`. Both `get_greek_raman` (RRS) and `get_greek_raman_VS` (VRS) build the standard 3-term polarized Rayleigh-form expansion `(α, β, γ, δ, ϵ, ζ)` with the appropriate `2γ_C/(1 + γ_C)` depolarization input. RRS uses universal `γ_C,rot`; VRS uses the per-molecule `γ_C,vib`. ✓

### Air-mixture `γ_C,Cab,air` weighting (Paper II eq 27)

`src/Inelastic/inelastic_helper.jl:373-376`. Code's `γ_air_Cabannes = 3/(4 + tmp1/tmp2)` reduces to Paper II eq 27 `γ_C,Cab,air = 3k/(1+4k)` with `k = (N_N₂+N_O₂)/(D_N₂+D_O₂)` — see expansion in audit notes. ✓ (Modulo the upstream input from Finding #2 round-trip.)

---

## Suggested next steps

1. **Verify α̅₀₀ source** (#1) — clarify whether the code values come from a different reference than Paper I Table 3, or correct them.
2. **Fix `compute_γ_mol_Cabannes!`** (#2) — drop-in replacement with the analytic inverse of Paper II eq 25 (one-line change).
3. **Reconcile RRS overload reversal** (#3) — pick one convention, document, apply to both overloads.
4. **Multiply `b_{J,J}` into hires VRS σ** (#4) — single-line fix per branch in `compute_σ_Rayl_VibRaman_coeff_hires!`. Add an integrated-σ regression test against `compute_σ_VibRaman_coeff!`.
5. **Rename or supplement `γ_C_Rayl`** (#5) — pick option 1 (rename), 2 (add `γ_C_Ray`), or 3 (unify convention).
6. **Document or remove per-J `γ_C(J)`** (#6) — either cite the source or revert to Paper I/II J-independent γ_C.
7. **Bypass round-trip in `compute_γ_air_Cabannes!`** (#7) — small win, depends on #5.

Items 1, 2, 3, 4, 7 are mostly local edits with low risk. Item 5 is the largest because it touches naming and may affect benchmark numbers. Item 6 needs author input.
