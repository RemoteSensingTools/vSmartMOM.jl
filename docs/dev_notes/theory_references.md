# vSmartMOM.jl — Theory Reference & Equation Crib Sheet

Verified mapping from the canonical papers (in `docs/papers/`) to source files
and target Concepts pages. This file is the **source of truth** for paper
citations and `file.jl:LINE` references that the Concepts arc draws from.

Built by reading each PDF in targeted mode (equations + design-decision
passages only, no benchmarks/intros). Updated 2026-04-30.

## Canonical references

| Tag | Citation | PDF in `docs/papers/` |
|---|---|---|
| **S2013** | Sanghavi, Martonchik, Davis, Diner (2013). *Linearization of a scalar matrix operator method radiative transfer model with respect to aerosol and surface properties.* JQSRT **116**:1–16. <https://doi.org/10.1016/j.jqsrt.2012.10.021> | `1-s2.0-S0022407312004633-main.pdf` |
| **S2014** | Sanghavi, Davis, Eldering (2014). *vSmartMOM: A vector matrix operator method-based radiative transfer model linearized with respect to aerosol properties.* JQSRT **133**:412–433. <https://doi.org/10.1016/j.jqsrt.2013.09.004> | `1-s2.0-S0022407313003592-main.pdf` |
| **SF2014** | Sanghavi (2014). *Revisiting the Fourier expansion of Mie scattering matrices in generalized spherical functions.* JQSRT **136**:16–27. | `1-s2.0-S0022407313004962-main.pdf` |
| **SS2015** | Sanghavi & Stephens (2015). *Adaptation of the delta-m and δ-fit truncation methods to vector radiative transfer.* JQSRT **159**:53–68. | `1-s2.0-S002240731500093X-main.pdf` |
| **S2022-I** | Sanghavi (2022). *Raman scattering in the earth's atmosphere, part I: Optical properties.* JQSRT **291**:108328. | `1-s2.0-S0022407322002631-main.pdf` |
| **SF2023-II** | Sanghavi & Frankenberg (2023). *Raman scattering in the Earth's atmosphere, Part II: Radiative transfer modeling for remote sensing applications.* JQSRT **311**:108791. | `1-s2.0-S0022407323003096-main.pdf` |
| **Fell1997** | Frank Fell (1997). *Validierung eines Modells zur Simulation des Strahlungstransportes in Atmosphäre und Ozean.* PhD Thesis, Freie Universität Berlin. (German.) | `1997_Thesis_Fell.pdf` |
| **JSF2022** | Jeyaram, Sanghavi, Frankenberg (2022). *vSmartMOM.jl: An open-source Julia package for atmospheric radiative transfer and remote sensing tools.* JOSS **7**(80):4575. <https://doi.org/10.21105/joss.04575> | (not in `docs/papers/` yet) |

**Primary citation for the package:** JSF2022 (JOSS software paper).
**Primary methodological citation:** S2014 (vector RT + linearization).

## Symbol convention (S2014 §2.1, used uniformly in the Concepts arc)

- `μ` = direction cosine of light propagation (`μ > 0` downward, `μ < 0` upward).
- `μ₀` = solar direction cosine; `φ₀` = solar azimuth.
- `τ` = optical depth (TOA → BOA), with subscript `_λ` for total (= absorption + scattering) and unsubscripted for scattering-only inside the layer-construction code.
- `ϖ_0` = single-scattering albedo (the paper uses `ω̄_0`; the code uses `ϖ` and `ϖ_λ`).
- `B(T)` = Planck source.
- `Z(μ,φ; μ',φ')` = phase matrix for Stokes vector scattering, decomposed in Fourier moments `m`.
- `M = diag(μ_i)`, `C = diag(c_i)` (quadrature weights), `E` = identity supermatrix.
- `D = diag(1, 1, −1, −1)` per stream — the Stokes polarization symmetry matrix.
- `δ` = elemental optical thickness; `f_tr` = δ-M forward-truncation factor.
- `Π^m_l(μ)` = generalized spherical-function matrices (4×4); `B_l` = Greek matrix with coefficients `(α_l, β_l, γ_l, δ_l, ϵ_l, ζ_l)`.

## Equation → code map

### A. Vector RTE & discretization (S2014 §2.1) → Concepts/02

| Paper Eq. | What it says | Source file:line |
|---|---|---|
| **S2014 (2)** | Plane-parallel vector RTE for Stokes vector `L(τ, μ, φ; μ₀, φ₀)`: thermal emission + direct-beam single scatter + diffuse multiple scatter. | (theory anchor; not coded directly) |
| **S2014 (4)** | Restored RTE for diffuse Stokes vector `L` after subtracting direct-beam delta. | `src/CoreRT/rt_run.jl:53–329` (the loop that solves it) |
| **S2014 (6)** | TOA boundary condition `L(0, μ, φ; μ₀, φ₀) = S₀·δ(μ−μ₀)·δ(φ−φ₀)` — direct solar source. | `src/CoreRT/rt_run.jl:189–193` (default F₀ initialization) |
| **S2014 (8)–(10)** | Phase matrix Fourier decomposition: `Z = (1/2)C₀ + Σₘ[C̄ₘ cos m(φ−φ₀) + S̄ₘ sin m(φ−φ₀)]`. Block-diagonal/anti-diagonal. | `src/Scattering/compute_Z_matrices.jl::compute_Z_moments` |
| **S2014 (12)** | Per-Fourier-moment RTE `μ dI_m/dτ̄ = −I_m + (ϖ₀/2)∫Z_m·I_m dμ' + (1−ϖ₀)B(T)δ_0m` — the form the solver actually iterates. | `src/CoreRT/rt_run.jl:208–215` (m loop) |
| **S2014 (14)–(15)** | Supermatrix form for upwelling/downwelling streams stacked into `Iₘ⁺`, `Iₘ⁻` of length `4·Nquad`. | `src/CoreRT/types.jl::CompositeLayer/AddedLayer` (3D arrays of shape `(NquadN, NquadN, nSpec)`) |
| **S2014 (17)** | Diagonal supermatrices `M = diag(μ_i)`, `C = diag(c_i)`. | `src/CoreRT/tools/rt_set_streams.jl::QuadPoints` |

### B. Elemental layer (S2014 §2.1, App. A; SF2023-II §3.2; Fell 1997 Eqs. 1.52–1.56) → Concepts/04

| Paper Eq. | What it says | Source file:line |
|---|---|---|
| **S2014 (18)** | Layer operator form: `⟨I_Δ⁻; I_0⁺⟩ = T_Δ I_0 + R_Δ I_Δ + J_Δ`. | `src/CoreRT/rt_run.jl::rt_kernel!` call site |
| **S2014 (19)** | `T_δ,m = (E − M⁻¹[E − (ϖ₀/2)Z̄ₘ⁽⁺⁺⁾C])·δ` — elemental transmission for `m=0` thermal-only. | `src/CoreRT/CoreKernel/elemental.jl::get_elem_rt!` lines 113–139 |
| **S2014 (20)** | `R_δ,m = M⁻¹·(ϖ₀/2)·Z̄ₘ⁽⁻⁺⁾·C·δ` — elemental reflection. | `src/CoreRT/CoreKernel/elemental.jl::get_elem_rt!` lines 113–139 |
| **S2014 (21)** | `J_δ,m = δ_0m·(E − T − R)·B(T)·1` — elemental thermal source (no solar SFI here, see §2.2). | `src/CoreRT/CoreKernel/elemental.jl::get_elem_rt_SFI!` |
| **S2014 (22)** | Elemental-thickness bound `0 < δ < min(μ_i)/(1 − ϖ₀/2)` — guarantees non-negative operators. | `src/CoreRT/tools/rt_helper_functions.jl::doubling_number`, consumed in `rt_kernel.jl::get_dtau_ndoubl` (lines 245–253) |
| **S2014 §2.2** | **Design choice (echo verbatim)**: vSmartMOM does *not* include solar single-scatter in `J`. Trade-off: simpler doubling/adding, no off-quadrature SFI integration; recovered via Block-Radau (App. B). | `src/CoreRT/CoreKernel/rt_kernel.jl` (no solar term in elemental J for `m>0` thermal-only path) |
| **SF2023-II (8)–(9)** | **Constant-`N_doubl` trick** for absorbing bands: `δτ = δτ_abs + δτ_scatt`; doubling sized by `δτ_scatt` only, so `N_doubl` is held constant across the spectral grid even when `τ_abs` swings over orders of magnitude. | `src/CoreRT/tools/atmo_prof.jl::construct_atm_layer` lines ~340–376; `rt_kernel.jl::get_dtau_ndoubl` lines 245–253 |
| **SF2023-II (10)–(11)** | Explicit elastic R, T, J for elemental layer with absorption baked into `δτ = δτ_abs + δτ_scatt`. Includes `1 − exp(−δτ(1/μ_i+1/μ_j))` and `exp(−δτ/μ_i) − exp(−δτ/μ_j)` factors. | `src/CoreRT/CoreKernel/elemental.jl` lines 207–252 (the `expm1` / `expdiff_neg` numerical-stability forms) |
| **Fell 1997 Eqs. 1.52–1.56** | Single-scatter elemental r/t/J formulas in the form vSmartMOM implements (matches SF2023-II (10)–(11)). | Already cited verbatim in `src/CoreRT/CoreKernel/elemental.jl` source comments |

### C. Doubling and adding (S2014 §2.1, Eqs. 23–32; SF2023-II §3.3) → Concepts/04

| Paper Eq. | What it says | Source file:line |
|---|---|---|
| **S2014 (23)** | `T_20 = T_21·(E − R_01·R_21)⁻¹·T_10` — doubling/adding T. | `src/CoreRT/CoreKernel/rt_helpers.jl::compute_geometric_progression!` lines 88–93; `interaction.jl:14–136` |
| **S2014 (24)** | `R_20 = R_10 + T_01·(E − R_21·R_01)⁻¹·R_21·T_10` — doubling/adding R. | `src/CoreRT/CoreKernel/rt_helpers.jl::doubling_rt_update!` lines 117–122; `interaction.jl::ScatteringInterface_11` |
| **S2014 (25)–(26)** | `T_02`, `R_02` — opposite-direction recurrences. | `interaction.jl:14–136` (full case) |
| **S2014 (27)–(28)** | `J_20`, `J_02` — source updates with the same `(E − R·R)⁻¹` geometric-series factors. | `src/CoreRT/CoreKernel/rt_helpers.jl::doubling_source_update!` lines 102–108 |
| **S2014 (29)–(30)** | **D-matrix symmetry for homogeneous layer**: `T_ab = D·T_ba·D`, `R_ab = D·R_ba·D`. | `src/CoreRT/CoreKernel/doubling.jl::apply_D!` lines 85–110 |
| **S2014 (31)** | `R*_10 = D·R_10` — the *starred* quantity used inside the fast doubling loop in place of both `R_10` and `R_01`. | `src/CoreRT/CoreKernel/doubling.jl::apply_D!` step 1 (the in-place row negation for `i>2`) |
| **S2014 (32)** | After the inner loop, recover the four operators: `T_ba = T_ba`, `R_ba = D·R*_ba`, `T_ab = D·T_ba·D`, `R_ab = R*_ba·D`. | `src/CoreRT/CoreKernel/doubling.jl::apply_D!` step 2 (the `(i,j)`-parity sign table) |
| **S2014 (33)–(37)** | Surface boundary `I_s⁺ = R_0s I_s⁻ + T_s0 I_0⁺ + J_0s`; emissivity boundary `I_s⁻ = ε_g B(Ts)·E + R_g M I_s⁺`; combined to give the closed system. | `src/CoreRT/Surfaces/*.jl::create_surface_layer!`; `interaction_hdrf.jl:1–42` |
| **SF2023-II (12)** | Restates S2014 (23)–(28) for the elastic component in the Raman-coupled formulation (same equations, different name). | `interaction.jl:14–136` |

### D. δ-M and δ-fit truncation (S2014 App. A; SS2015 §2–§3) → Concepts/03c

| Paper Eq. | What it says | Source file:line |
|---|---|---|
| **S2014 (A.3)** | Truncated phase matrix `Z*(τ*; μ, μ') = (Z(...) − β_T·δ(μ−μ')·E) / (1 − β_T)`. | `src/CoreRT/LayerOpticalProperties/delta_m_truncation.jl` |
| **S2014 (A.5)–(A.10)** | Greek coefficients of the Dirac delta `δ(μ−1)·E`: `β_l^δ = δ_l^δ = ζ_l^δ = α_l^δ = (2l+1)/2`; `γ_l^δ = ϵ_l^δ = 0`. | (theory anchor) |
| **S2014 (A.11)** | Truncation factor `β_T = (2/(2L+1))·β_L` — the forward-peak fraction at truncation order `L`. | `src/Scattering/types.jl::AerosolOptics.fᵗ` |
| **S2014 (A.12)** | Modified Greek matrix `B*_l = (B_l − β_T·B^δ_l) / (1 − β_T)`. | `src/CoreRT/LayerOpticalProperties/delta_m_truncation.jl` lines 44–48; `compEffectiveLayerProperties.jl::createAero` lines 67–72 |
| **SS2015 (8)** | Truncation rescaling: `τ* = τ(1 − f_tr·ϖ_0)`, `ϖ_0* = ϖ_0(1 − f_tr) / (1 − f_tr·ϖ_0)`, `Z* = (Z − f_tr·δ·E) / (1 − f_tr)`. | `compEffectiveLayerProperties.jl::createAero` lines 67–72 (the explicit `τ_mod = (1−fᵗ·ω̃)·τAer; ϖ_mod = (1−fᵗ)·ω̃/(1−fᵗ·ω̃)` form) |
| **SS2015 (26)** | `f_tr = β_{L_tr}/(2L_tr+1)` — Wiscombe choice of truncation factor. | `src/Scattering/truncate_phase.jl` |
| **SS2015 (27a–f)** | **Modified Greek coefficients** `β_l*, δ_l*, γ_l*, ϵ_l*, α_l*, ζ_l*` at truncation order `L_tr`. The exact formulas to quote in Concepts/03c. | `src/Scattering/truncate_phase.jl`; `src/CoreRT/LayerOpticalProperties/delta_m_truncation.jl` |
| **SS2015 §2.4** | Two error sources: **DSE** (delta-separation error, from setting all `B_l = 0` for `l > L_tr`) and **PTE** (phase-truncation error, from approximating the forward peak as a Dirac delta). | (design-rationale callout for Concepts/03c) |
| **SS2015 §3, Eqs. (38a–d)** | **δBGE-fit** (δ Beta-Gamma-Epsilon fit) — vector extension of Hu et al. δ-fit; `β_l*, δ_l*, α_l*, ζ_l*` from SVD of a least-squares system on `b_1(μ), b_2(μ)`. Used as default for accurate hyperspectral retrievals near backscatter. | `src/CoreRT/LayerOpticalProperties/` (δBGE-fit selector) |
| **SS2015 §4.1** | **Critical for retrievals**: δ-m distorts O₂ A-band line shapes near exact backscatter (ϑ_view ≈ −60° in principal plane); δBGE-fit eliminates this. Worth one-paragraph callout in Concepts/03c. | (design-rationale callout) |
| **SS2015 §5** | TMS (Nakajima & Tanaka) single-scattering correction reduces remaining truncation error by ~3× near exact backscatter. | (rt_run_ss connection — Concepts/04 also) |

### E. Quadrature (S2014 App. B) → Concepts/02

| Paper Eq. | What it says | Source file:line |
|---|---|---|
| **S2014 (B.1)–(B.2)** | Block-Radau scheme: composite `(N₁+N₂)`-point quadrature where solar `μ₀` is included as a quadrature root with non-zero weight (a true Radau node), and viewing angles `μ_v` are appended as zero-weight dummy nodes. Eliminates interpolation error between quadrature points for direct ray tracing. | `src/CoreRT/tools/rt_set_streams.jl::RadauQuad` lines 24–110 |
| **S2014 App. B (text)** | Without DNI/Radau, off-quadrature solar/viewing angles need interpolation between roots → error. The Radau construction is the cost paid for the §2.2 design choice (no solar SFI in J). | (design-rationale callout for Concepts/02) |

### F. Layer optical-property assembly → Concepts/03 + 03a + 03b + 03c

| Concept | Source file:line |
|---|---|
| Build `(τ_λ, ϖ_λ, Z⁺⁺, Z⁻⁺)` per layer per Fourier moment | `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl::constructCoreOpticalProperties` lines 11–65 |
| τϖ-weighted mixing of scatterers via `+` | `src/CoreRT/types.jl::Base.:+(::CoreScatteringOpticalProperties, ...)` lines 1063–1093 |
| Vertical concatenation of layers via `*` | `src/CoreRT/types.jl::Base.:*(::CoreScatteringOpticalProperties, ...)` lines 1096+ |
| Cabannes vs full-Rayleigh greek selection | `compEffectiveLayerProperties.jl::_rayleigh_greek_source` lines 8–9 |
| Aerosol Mie → AerosolOptics(GreekCoefs, ω̃, k, fᵗ) | `src/Scattering/compute_NAI2.jl:44`, `src/Scattering/compute_PCW.jl:28` |
| Phase matrix from Greek + μ for Fourier moment m | `src/Scattering/compute_Z_matrices.jl::compute_Z_moments` |
| Gas absorption τ_abs(ν, layer) | `src/Absorption/compute_absorption_cross_section.jl` lines 32–280 |
| Adding gas absorption to scattering total | `compEffectiveLayerProperties.jl` line 58 (`combo + CoreAbsorptionOpticalProperties`) |

### G. Mie / Greek coefficients (SF2014 / S2014 §A) → Concepts/03b

| Paper Eq. | What it says | Source file:line |
|---|---|---|
| **SF2014 (1)** | Mie cross sections from `aₙ`, `bₙ` series. | `src/Scattering/mie_helper_functions.jl` |
| **SF2014 (3)–(6)** | Scattering matrix elements `f_ij` from Wigner d-functions `π_n`, `τ_n`. | (theory anchor, code uses Greek directly) |
| **SF2014 (12)–(14)** | **Phase matrix Z Fourier decomposition**: `Z = (1/2)C⁰ + Σₘ[Cᵐ cos m(φ−φ') + Sᵐ sin m(φ−φ')]`, with `Cᵐ = Aᵐ·Δ + Δ·Aᵐ·Δ`, `Aᵐ = Σ_{l=m}^∞ Π^m_l(μ)·B_l·Π^m_l(μ')`. | `compute_Z_matrices.jl::compute_Z_moments` |
| **SF2014 (15)** | Generalized spherical-harmonic matrices `Π^m_l(μ)` (4×4 with `P^m_l`, `R^m_l`, `T^m_l`). | `compute_Z_matrices.jl` |
| **SF2014 (16)** | **Greek matrix B_l** with the six coefficients `(α_l, β_l, γ_l, δ_l, ϵ_l, ζ_l)` populating its 4×4 form. | `src/Scattering/types.jl::GreekCoefs` lines 231–244 |
| **SF2014 (17)** | Six integral formulas for the Greek coefficients in terms of `f_ij` integrals — used by NAI-2 path. | `src/Scattering/compute_NAI2.jl` |
| **SF2014 (24)** | Greek matrix `B_l` from precomputed Wigner-3j combinations of `(aₙ, bₙ)` products — used by PCW path. | `src/Scattering/compute_PCW.jl` |
| **SF2014 §5** | **NAI-2 vs PCW choice**: NAI-2 ≥13× faster than NAI-1 for individual large particles (>228 size parameter); PCW 6–8× faster than NAI-2 for polydispersions over a finite size range. NAI-2 is the safer default; PCW dominates for big polydisperse populations. | `src/Scattering/Scattering.jl` (DECOMP_MAP entry) |

### H. Linearization (S2013 §3, S2014 App. C) → Concepts/06

| Paper Eq. | What it says | Source file:line |
|---|---|---|
| **S2014 (C.1)–(C.4)** | Differentiated RTE and boundary conditions; linearization is exact at the operator level. | `src/CoreRT/rt_run_lin.jl` |
| **S2014 (C.5)–(C.7)** | Differentiation rules for matrix products and inverses (used everywhere in the chain rule). | (foundation) |
| **S2014 (C.8)–(C.10)** | **Elemental derivatives**: `Ṫ_δ,m`, `Ṙ_δ,m`, `J̇_δ,m` w.r.t. the three core layer variables `(τ, ϖ, Z)`. | `src/CoreRT/CoreKernel/elemental_lin.jl` |
| **S2014 (C.11)–(C.16)** | **Derivatives propagated through doubling/adding**: same shape as Eqs. (23)–(28) but for tangent-linear operators. | `src/CoreRT/CoreKernel/{doubling,interaction}_lin.jl` |
| **S2014 (C.17)–(C.18)** | **D-matrix symmetry on derivatives**: `Ṫ_ab = D·Ṫ_ba·D`, `Ṙ_ab = D·Ṙ_ba·D` — halves the linearized doubling cost too. | `doubling_lin.jl::apply_D_lin!` (or equivalent) |
| **S2014 (C.19)** | `Ṙ*_10 = D·Ṙ_10` — the fast linearized doubling counterpart of (31). | `doubling_lin.jl` |
| **S2014 (C.21)** | Final assembled derivative form `⟨İ_Δ⁻; İ_0⁺⟩ = J̇_Δ + Ṙ_Δ I_Δ + T_Δ İ_0 + …`. | `src/CoreRT/CoreKernel/lin_added_layer_all_params.jl` lines 1–100 |
| **S2014 (C.22)–(C.24)** | Layer-averaged optics `τ̄ = τ_r + τ_g + Σ_i τ_i`, `ϖ̄_0 = (τ_r + Σ ϖ_0,i τ_i)/τ̄`, `Z̄ = (τ_r Z_r + Σ τ_i ϖ_0,i Z_i)/(ϖ̄_0 τ̄)`. | `compEffectiveLayerProperties.jl` |
| **S2014 (C.25)–(C.26)** | **Chain rule** through the elemental SS variables `x_SS = (τ_i, ϖ_0,i, Z_i, β_i)` to the microphysical parameters `x_µ = (n_r,i, n_i,i, r_m,i, σ_i)`. | `lin_added_layer_all_params.jl`; `src/Scattering/types_lin.jl` |
| **S2014 (C.27)** | `δ = τ̄ / 2^N_dbl` — elemental thickness used in derivatives. | `rt_kernel.jl::get_dtau_ndoubl` |
| **S2014 (C.28)–(C.31)** | Derivatives of `δ` w.r.t. `τ_i`, `ϖ_0,i`, `Z_i`, `β_i`. | `compEffectiveLayerProperties_lin.jl` |
| **S2014 (C.32)–(C.39)** | Derivatives of `ϖ̄_0` and `Z̄` (post-truncation forms). | same |
| **S2014 (C.40)** | `Ż_m(μ_i, μ_j)` formula via generalized spherical harmonics. | `compute_Z_matrices.jl` (linearized variant) |
| **S2014 (C.41)–(C.42)** | `β̇*` and `Ḃ*_l` for the truncated case — chain rule through the δ-M factor `f_tr`. | `delta_m_truncation_lin.jl` |
| **S2013 (35)–(46)** | Scalar predecessor of S2014 App. C — same structure for scalar `(τ, ϖ_0, P)`. Useful as the simpler-derivation reference for readers learning the chain rule. | (theory anchor for Concepts/06 history paragraph) |
| **S2013 (61)–(67)** | δ-M chain rule for scalar (`P_mod`, `τ_mod`, `ϖ_mod`). | `compEffectiveLayerProperties_lin.jl` |
| **S2013 (90)** | `N_aer = τ_λ / k_λ` — number density invariant w.r.t. wavelength. | `compEffectiveLayerProperties.jl` |
| **S2014 §C.2.2 / S2013 §3.2.2** | Surface BRDF derivatives (Lambertian, mRPV, Cox-Munk wind speed). | `Surfaces/{lambertian,rpv,coxmunk}_surface_lin.jl` |

### I. Inelastic / Raman (S2022-I, SF2023-II) → Concepts/08 (BRIEF)

User asked to keep this round light on RRS. Crib-sheet entries here are
sufficient for the brief Concepts/08 page; do not expand.

| Paper Eq. | What it says | Source file |
|---|---|---|
| **S2022-I (12), (16)** | Cabannes depolarization `γ_C,Cab = 3γ²/(45a²+4γ²)`, Rayleigh σ. | `src/Inelastic/types.jl` |
| **S2022-I (20)–(22)** | RRS phase matrix and σ_RRS, `b_{J1,J0}` coefficients. | `src/Inelastic/` |
| **S2022-I (28)–(32)** | VRS phase matrix and σ_VRS for vibrational Raman. | `src/Inelastic/` |
| **S2022-I (33)–(34)** | RVRS phase matrix and σ_RVRS for rovibrational. | `src/Inelastic/` |
| **S2022-I (35)–(36)** | Dunham expansion for ν̃_fi, energy levels. | `src/Inelastic/` |
| **S2022-I Tables 1–5** | Dunham coefficients, nuclear spin degeneracy `g_J`, polarizability invariants `a_B`, `γ_B`. | `src/Inelastic/` constants |
| **SF2023-II (5)–(6)** | **Coupled RTE** at λ and λ_r. | `src/CoreRT/CoreKernel/` (Raman dispatch) |
| **SF2023-II (8)–(9)** | **Constant-`N_doubl` trick** in absorbing bands (already covered in Concepts/03c — referenced here only to note it benefits Raman too). | `compEffectiveLayerProperties.jl`; `rt_kernel.jl` |
| **SF2023-II (14)–(15)** | Inelastic elemental R, T, J for `λ → λ_r`. | `src/CoreRT/CoreKernel/elemental_inelastic.jl` |
| **SF2023-II (16)–(21)** | Inelastic adding/doubling. **Linear-in-inelastic-scattering approximation** (one inelastic event per photon path: E-E-E-I or E-E-I-E patterns). | `interaction_inelastic.jl`, `doubling_inelastic.jl` |
| **SF2023-II (32)** | Single-scattering inelastic correction `I₁ ≈ I₀ + (I'₁ − I'₀)`. Justifies `rt_run_ss` alongside `rt_run`. | `src/CoreRT/rt_run.jl::rt_run_ss` lines 364–524 |

### J. Architecture-agnostic / GPU evidence anchors → Concepts/07

(Not paper equations — but the differentiator section needs concrete code
citations. Listed here so the Concepts/07 page can grep this file.)

| Claim | Source file:line |
|---|---|
| Architecture types `CPU`/`GPU`/`MetalGPU` | `src/Architectures.jl:33–96` |
| Weak-dep CUDA injection at load time | `ext/vSmartMOMCUDAExt.jl:21–66` |
| Weak-dep Metal injection at load time | `ext/vSmartMOMMetalExt.jl:19–101` |
| `@kernel` apply_D! (small but representative) | `src/CoreRT/CoreKernel/doubling.jl:85–110` |
| `@kernel` line shape (Doppler/Lorentz/Voigt) | `src/Absorption/compute_absorption_cross_section.jl:229–280` |
| `@kernel` Mie coefficients (DoubleSingle precision for FP32 GPU) | `src/Scattering/gpu_mie_kernels.jl:30–100` |
| `@kernel` portable LU inverse with `@localmem` | `src/CoreRT/tools/ka_batched_kernels.jl:100–178` |
| Batched matmul CPU (threaded BLAS) | `src/CoreRT/tools/cpu_batched.jl:24–73` |
| Batched matmul CUDA (CUBLAS strided) | `ext/gpu_batched_cuda.jl:122–139` |
| Batched matmul Metal (KA portable) | `ext/vSmartMOMMetalExt.jl:29–54` |
| ForwardDiff.Dual through CUDA batched_mul | `ext/gpu_batched_cuda.jl:141–177` |
| Allocation by architecture | `src/CoreRT/tools/rt_helper_functions.jl::make_added_layer` lines 91–142 |
| Spectral axis as third dimension | `src/CoreRT/types.jl::CompositeLayer` lines 150–163 (3D arrays `(NquadN, NquadN, nSpec)`) |

## "Design choice" passages to echo verbatim in the docs

These are not equations but explicit decisions made in the papers that the
docs should quote (or paraphrase tightly) so users understand *why* the code
is shaped as it is:

1. **No solar SFI in J** (S2014 §2.2). The matrix-operator Eqs. (23)–(28)
   suffice instead of (23)–(33); for thermal-only `m=0`, `J` is isotropic
   so its computations can be reused. The cost is needing Block-Radau
   (App. B) for off-quadrature angles. → Concepts/02 + Concepts/04.
2. **Constant-`N_doubl` across the spectral grid** (SF2023-II §3.2,
   Eqs. 8–9). `N_doubl` is sized by the scattering optical depth, not the
   total. Lets line-by-line on hyperspectral grids run on GPU at one
   batched call per layer per Fourier moment. → Concepts/03c +
   Concepts/04 + Concepts/07.
3. **Cabannes vs full-Rayleigh greek source** (`compEffectiveLayerProperties.jl:8–9`).
   Pure-elastic uses the higher-depol Rayleigh greek (rotational Raman is
   rolled into effective depolarization). Raman-aware modes use Cabannes
   greek and handle rotational Raman explicitly. The mismatch costs ~1%
   on Stokes I and ~3% on Q. → Concepts/03b + Concepts/08.
4. **Linear-in-inelastic-scattering approximation** (SF2023-II §3.4). Only
   one inelastic event per photon path; second-order has been shown
   negligible. → Concepts/08.
5. **Why the package even has `rt_run_ss`** (SF2023-II §5.4.2, Eq. 32).
   Not a debug helper; it's the single-scattering inelastic correction
   that makes the linear-in-Raman approximation usable in absorbing
   bands. → Concepts/04 (mention) + Concepts/08 (mention).
6. **δBGE-fit over δ-m for hyperspectral retrievals near backscatter**
   (SS2015 §4.1). δ-m distorts O₂ A-band line shapes near `ϑ_view ≈ −60°`
   in the principal plane, by an amount that biases XCO₂ retrievals.
   → Concepts/03c.
7. **Exact finite-δ elemental, not the linear approximation.** S2014
   Eqs. (19)–(20) are written in the **infinitesimal-δ limit** (first-order
   in δ, equivalent to `1−exp(−x) ≈ x` and `exp(−x) ≈ 1−x`). Many MOM codes
   stop there — they need a very thin elemental layer (large `N_doubl`)
   for that linear form to be accurate. vSmartMOM's `elemental.jl` instead
   implements the **exact finite-δ single-scatter formulas** (Fell 1997
   Eqs. 1.52–1.56, restated as SF2023-II Eqs. 10–11):
   ```
   r⁻⁺[i,j] = ϖ_λ·Z⁻⁺[i,j]·(μⱼ/(μᵢ+μⱼ))·wⱼ·(1 − exp(−δτ(1/μᵢ + 1/μⱼ)))
   t⁺⁺[i,j] = ϖ_λ·Z⁺⁺[i,j]·(μⱼ/(μᵢ−μⱼ))·wⱼ·(exp(−δτ/μᵢ) − exp(−δτ/μⱼ))   (i ≠ j)
   t⁺⁺[i,i] = exp(−δτ/μᵢ)·(1 + ϖ_λ·Z⁺⁺[i,i]·(δτ/μᵢ)·wᵢ)                  (i = j L'Hôpital limit)
   ```
   coded with `-expm1(-x)` and `expdiff_neg(a, b)` for numerical stability
   in optically thin layers and especially in `Float32`. This means the
   elemental layer can be **thicker at the same single-scatter accuracy**
   → fewer doublings, less round-off accumulation. It compounds with the
   constant-`N_doubl` trick: `N_doubl` is sized by single-scatter accuracy
   on `τ_scat`, not by Taylor-series truncation on the exponential.
   → Concepts/04 § Elemental (load-bearing); Concepts/01 (differentiator
   list).

## CITATION.bib entries to add (for `pages/vSmartMOM/References.md`)

```bibtex
@article{Sanghavi2014vSmartMOM,
  author = {Sanghavi, S. and Davis, A. B. and Eldering, A.},
  title = {vSmartMOM: A vector matrix operator method-based radiative transfer model linearized with respect to aerosol properties},
  journal = {JQSRT},
  volume = {133},
  pages = {412--433},
  year = {2014},
  doi = {10.1016/j.jqsrt.2013.09.004}
}

@article{Sanghavi2013smartMOM,
  author = {Sanghavi, S. and Martonchik, J. V. and Davis, A. B. and Diner, D. J.},
  title = {Linearization of a scalar matrix operator method radiative transfer model with respect to aerosol and surface properties},
  journal = {JQSRT},
  volume = {116},
  pages = {1--16},
  year = {2013},
  doi = {10.1016/j.jqsrt.2012.10.021}
}

@article{Sanghavi2014MieFourier,
  author = {Sanghavi, S.},
  title = {Revisiting the Fourier expansion of Mie scattering matrices in generalized spherical functions},
  journal = {JQSRT},
  volume = {136},
  pages = {16--27},
  year = {2014}
}

@article{SanghaviStephens2015,
  author = {Sanghavi, S. and Stephens, G.},
  title = {Adaptation of the delta-m and delta-fit truncation methods to vector radiative transfer: Effect of truncation on radiative transfer accuracy},
  journal = {JQSRT},
  volume = {159},
  pages = {53--68},
  year = {2015}
}

@article{Sanghavi2022Raman1,
  author = {Sanghavi, S.},
  title = {Raman scattering in the earth's atmosphere, part I: Optical properties},
  journal = {JQSRT},
  volume = {291},
  pages = {108328},
  year = {2022}
}

@article{SanghaviFrankenberg2023Raman2,
  author = {Sanghavi, S. and Frankenberg, C.},
  title = {Raman scattering in the Earth's atmosphere, Part II: Radiative transfer modeling for remote sensing applications},
  journal = {JQSRT},
  volume = {311},
  pages = {108791},
  year = {2023}
}

@article{Jeyaram2022vSmartMOMjl,
  author = {Jeyaram, R. and Sanghavi, S. and Frankenberg, C.},
  title = {vSmartMOM.jl: An open-source Julia package for atmospheric radiative transfer and remote sensing tools},
  journal = {Journal of Open Source Software},
  volume = {7},
  number = {80},
  pages = {4575},
  year = {2022},
  doi = {10.21105/joss.04575}
}

@phdthesis{Fell1997,
  author = {Fell, F.},
  title = {Validierung eines Modells zur Simulation des Strahlungstransportes in Atmosph\"are und Ozean},
  school = {Freie Universit\"at Berlin},
  year = {1997}
}
```
