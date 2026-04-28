# vSmartMOM.jl — Theory Reference Map

Crib sheet mapping the canonical Sanghavi papers to source files and to where
they should be cited in the documentation. Use this when writing the
"Core RT Theory" page (`docs/src/pages/vSmartMOM/CoreRTTheory.md`), the
inelastic page, the references page, and `CITATION.bib`.

## Canonical references

| Tag | Citation |
|---|---|
| **S2014** | Sanghavi, S., Davis, A.B., Eldering, A. (2014). "vSmartMOM: A vector matrix operator method-based radiative transfer model linearized with respect to aerosol properties." *JQSRT* 133:412–433. <https://doi.org/10.1016/j.jqsrt.2013.09.004> |
| **S2015** | Sanghavi, S., Stephens, G. (2015). "Adaptation of the delta-m and δ-fit truncation methods to vector radiative transfer: Effect of truncation on radiative transfer accuracy." *JQSRT* 159:53–68. |
| **S2022-I** | Sanghavi, S. (2022). "Raman scattering in the earth's atmosphere, part I: Optical properties." *JQSRT* 291:108328. <https://doi.org/10.1016/j.jqsrt.2022.108328> |
| **SF2023-II** | Sanghavi, S., Frankenberg, C. (2023). "Raman scattering in the Earth's atmosphere, Part II: Radiative transfer modeling for remote sensing applications." *JQSRT* 311:108791. <https://doi.org/10.1016/j.jqsrt.2023.108791> |
| **JSF2022** | Jeyaram, R., Sanghavi, S., Frankenberg, C. (2022). "vSmartMOM.jl: An open-source Julia package for atmospheric radiative transfer and remote sensing tools." *JOSS* 7(80):4575. <https://doi.org/10.21105/joss.04575> *(software paper — primary citation for the package itself; not in `~/papers/` yet, fetch and add)* |
| **S2013** (older) | Sanghavi, S., Martonchik, J.V., Davis, A.B., Diner, D.J. (2013). "Linearization of a scalar matrix operator method radiative transfer model with respect to aerosol and surface properties." *JQSRT* 116:1–16. *(scalar predecessor, largely superseded by S2014)* |

## Equation → code map

### Elastic vector RT (S2014)

| Paper | Equation | Source file | What it does |
|---|---|---|---|
| S2014 | (4) | — | RTE for diffuse + direct combined intensity |
| S2014 | (18) | [src/CoreRT/CoreKernel/rt_kernel.jl](src/CoreRT/CoreKernel/rt_kernel.jl) | Layer interaction supermatrix form |
| S2014 | (19)–(20) | [elemental.jl](src/CoreRT/CoreKernel/elemental.jl) | Elemental T, R for thin layer |
| S2014 | (21) | [elemental.jl](src/CoreRT/CoreKernel/elemental.jl) | Elemental source J (thermal-only by S2014's design choice) |
| S2014 | (22) | [LayerOpticalProperties/](src/CoreRT/LayerOpticalProperties/) | Elemental-layer optical-thickness bound `δ < min(μ)/(1−ω/2)` |
| S2014 | (23)–(28) | [doubling.jl](src/CoreRT/CoreKernel/doubling.jl), [interaction.jl](src/CoreRT/CoreKernel/interaction.jl) | Doubling and adding (T20, R20, T02, R02, J20, J02) |
| S2014 | (29)–(32) | doubling/interaction kernels | D-matrix symmetry that halves doubling cost |
| S2014 | (33)–(37) | [src/CoreRT/Surfaces/](src/CoreRT/Surfaces/) | Surface boundary, TOA boundary |
| S2014 | §2.2 | [rt_kernel.jl](src/CoreRT/CoreKernel/rt_kernel.jl) | **Design rationale** — no solar single-scattering in J term, full direct+diffuse handled together. Worth one paragraph in CoreRTTheory.md because it's a deliberate departure from VLIDORT/SCIATRAN and explains why `J` is thermal-only here |

### Truncation (S2014 App. A + S2015)

| Paper | Equation | Source file |
|---|---|---|
| S2014 | (A.1)–(A.12) | [src/CoreRT/LayerOpticalProperties/](src/CoreRT/LayerOpticalProperties/) δ-m vector truncation |
| S2015 | full | δ-fit refinement, polarization-aware truncation discussion |

### Quadrature (S2014 App. B)

| Paper | Equation | Source file |
|---|---|---|
| S2014 | (B.1)–(B.2) | `QuadPoints` construction, QUAD_MAP entries in [src/IO/Parameters.jl](src/IO/Parameters.jl) — block-Radau direct-raytracing avoids interpolation between quadrature points. **Non-obvious technique that should cite S2014 App. B at the call site** |

### Linearization (S2014 App. C)

| Paper | Equation | Source file |
|---|---|---|
| S2014 | (C.8)–(C.10) | [src/CoreRT/CoreKernel/elemental_lin.jl](src/CoreRT/CoreKernel/elemental_lin.jl) — derivatives of T, R, J for elemental layer |
| S2014 | (C.11)–(C.16) | [doubling_lin.jl](src/CoreRT/CoreKernel/doubling_lin.jl), [interaction_lin.jl](src/CoreRT/CoreKernel/interaction_lin.jl) — derivatives propagated through doubling/adding |
| S2014 | (C.17)–(C.20) | D-matrix symmetry on derivatives (halves linearized doubling cost) |
| S2014 | (C.21) | [lin_added_layer_all_params.jl](src/CoreRT/CoreKernel/lin_added_layer_all_params.jl) — assembled derivative form |
| S2014 | (C.22)–(C.45) | [src/CoreRT/parameter_layout.jl](src/CoreRT/parameter_layout.jl) + [src/Scattering/types_lin.jl](src/Scattering/types_lin.jl) — chain rule for τ, ω, Z, β with respect to aerosol microphysical parameters; Mie linearization |

### Inelastic optical properties (S2022-I)

| Paper | Equation | Source file |
|---|---|---|
| S2022-I | (12), (16) | [src/Inelastic/types.jl](src/Inelastic/types.jl) | Cabannes depolarization γ_C,Cab, Rayleigh σ_Rayl |
| S2022-I | (20)–(22) | [src/Inelastic/](src/Inelastic/) | RRS phase matrix and σ_RRS, b_{J1,J0} coefficients (Eq. 19) |
| S2022-I | (28)–(32) | [src/Inelastic/](src/Inelastic/) | VRS phase matrix and σ_VRS |
| S2022-I | (33)–(34) | [src/Inelastic/](src/Inelastic/) | RVRS phase matrix and σ_RVRS |
| S2022-I | (35)–(36) | [src/Inelastic/](src/Inelastic/) | Dunham expansion for ν̃_fi, energy levels |
| S2022-I | Tables 1–5 | [src/Inelastic/](src/Inelastic/) constants | Dunham coefficients, nuclear spin degeneracy gJ, polarizability invariants a_B, γ_B from Buldakov |

### Inelastic RT (SF2023-II)

| Paper | Equation | Source file |
|---|---|---|
| SF2023-II | (5)–(6) | [src/CoreRT/CoreKernel/](src/CoreRT/CoreKernel/) | Coupled RTE at λ and λ_r |
| SF2023-II | (8)–(9) | [src/CoreRT/LayerOpticalProperties/](src/CoreRT/LayerOpticalProperties/) | **Constant-N_doubl elemental layer trick**: δ_scatt held fixed across spectral grid in absorbing bands. Currently undocumented anywhere user-visible — worth a paragraph in the CoreRTTheory page since it's the rationale behind the layer-construction code |
| SF2023-II | (10)–(11) | [elemental.jl](src/CoreRT/CoreKernel/elemental.jl) | Elastic elemental R, T, J (matches S2014 with explicit form) |
| SF2023-II | (12) | [doubling.jl](src/CoreRT/CoreKernel/doubling.jl), [interaction.jl](src/CoreRT/CoreKernel/interaction.jl) | Elastic doubling/adding |
| SF2023-II | (14)–(15) | [elemental_inelastic.jl](src/CoreRT/CoreKernel/elemental_inelastic.jl) | Inelastic elemental R, T, J for λ → λ_r |
| SF2023-II | (16)–(21) | [interaction_inelastic.jl](src/CoreRT/CoreKernel/interaction_inelastic.jl), [doubling_inelastic.jl](src/CoreRT/CoreKernel/doubling_inelastic.jl) | Inelastic adding/doubling: linear-in-inelastic-scattering assumption (one inelastic event per photon path) |
| SF2023-II | (32) | [src/CoreRT/rt_run.jl](src/CoreRT/rt_run.jl) `rt_run_ss` | **Single-scattering inelastic correction** `I₁ ≈ I₀ + (I'₁ − I'₀)`. Justifies why `rt_run_ss` exists alongside `rt_run` — it is not just a debug helper. The API reference should make this connection explicit |

## Suggested CoreRTTheory.md outline

1. **Overview** — what the matrix operator method does, why MOM (S2014 §2 intro). Cite S2014, JSF2022.
2. **Vector RTE** — Eq. (4) of S2014, sign convention for μ, polarization basis.
3. **Elemental layer** — S2014 Eqs. (19)–(22). Show file:line refs to [elemental.jl](src/CoreRT/CoreKernel/elemental.jl).
4. **Doubling and adding** — S2014 Eqs. (23)–(32). file:line refs to [doubling.jl](src/CoreRT/CoreKernel/doubling.jl), [interaction.jl](src/CoreRT/CoreKernel/interaction.jl). Mention D-matrix symmetry.
5. **Why no solar SFI in J** — S2014 §2.2. One-paragraph design-decision callout.
6. **Truncation** — S2014 Appendix A + S2015. Pointer to [LayerOpticalProperties/](src/CoreRT/LayerOpticalProperties/).
7. **Quadrature** — S2014 Appendix B. Pointer to QuadPoints + IO/Parameters.jl QUAD_MAP.
8. **Linearization (Jacobians)** — S2014 Appendix C. File refs to `_lin.jl` kernels and `parameter_layout.jl`.
9. **Inelastic extension** — SF2023-II §3. File refs to `*_inelastic.jl` kernels. Constant-N_doubl trick (§3.2) gets its own subsection.
10. **Single-scattering inelastic correction** — SF2023-II §5.4.2 + Eq. (32). Explicit pointer to `rt_run_ss`.
11. **Other dispatch arms** — name the `_ss`, `_canopy`, `_multisensor`, `_hdrf` variants and where to look for each.

Each equation block should show `file.jl:LINE` next to it so readers can grep directly into the source.

## CITATION.bib

The package's primary citation is **JSF2022** (the JOSS software paper). The
methodological citations are **S2014** (vector RT + linearization),
**S2022-I** (Raman optical properties), and **SF2023-II** (Raman RT in
vSmartMOM). **S2015** for truncation. Add all five to `CITATION.bib`.
