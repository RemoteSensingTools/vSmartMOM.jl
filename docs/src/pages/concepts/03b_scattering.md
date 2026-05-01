# 3b · Mie, Rayleigh, and the Phase Matrices

> **For:** users adding aerosols to their scene; method developers extending the scattering codebase; retrieval developers reasoning about polarization sensitivity.
>
> **Prev:** [3a · Gas Absorption](03a_absorption.md) · **Next:** [3c · Mixing & δ-M Truncation](03c_mixing.md)

This page covers three of the four layer-optics arrays plus the aerosol
contribution to the optical depth: the Fourier-moment phase matrices
``\mathbf{Z}^{++}`` and ``\mathbf{Z}^{-+}``, the aerosol single-scattering
albedo ``\tilde{\varpi}``, and the extinction cross-section ``k`` (which
maps microphysics to ``\tau_\mathrm{aer}``).

## What this page produces

For each aerosol type and each Fourier moment ``m``:

| Symbol | Meaning | Used for |
|---|---|---|
| ``\mathbf{Z}^{++}_m(\mu', \mu)``, ``\mathbf{Z}^{-+}_m(\mu', \mu)`` | phase matrix Fourier moments | layer scattering |
| ``\tilde{\varpi}`` | aerosol single-scattering albedo | layer SSA |
| ``k`` | extinction cross-section | aerosol optical depth from microphysics |
| ``f^t`` | δ-M forward-truncation factor | consumed in [Concepts/03c](03c_mixing.md) |

For Rayleigh (and Cabannes), the same ``\mathbf{Z}^{++}_m, \mathbf{Z}^{-+}_m`` are
computed from the Rayleigh / Cabannes Greek source.

## The Mie pipeline

```
   Aerosol(size_dist, n_r, n_i)
        │
        ▼
    MieModel
        │
        ▼
   compute_aerosol_optical_properties
        │
        ▼
   AerosolOptics(GreekCoefs, ω̃, k, fᵗ)
        │
        ▼
   compute_Z_moments  (per Fourier moment m)
        │
        ▼
   (Z⁺⁺, Z⁻⁺) per (μ', μ, m)
```

The user-facing entry is `compute_aerosol_optical_properties(model::MieModel)`
in `src/Scattering/`. The dispatcher branches on the decomposition algorithm
(`NAI2()` or `PCW()`) carried in `MieModel`.

## The Greek matrix `B_l`

Polarized Mie scattering is fully captured by the Greek expansion of the
``4 \times 4`` phase matrix (Sanghavi 2014 App. A; Sanghavi 2014 Mie-Fourier
Eq. 16):

```math
\mathbf{B}_l =
\begin{bmatrix}
\beta_l & \gamma_l & 0 & 0 \\
\gamma_l & \alpha_l & 0 & 0 \\
0 & 0 & \zeta_l & -\epsilon_l \\
0 & 0 & \epsilon_l & \delta_l
\end{bmatrix}
```

The six coefficients ``(\alpha_l, \beta_l, \gamma_l, \delta_l, \epsilon_l, \zeta_l)``
fully determine the polarized scattering pattern — ``\beta_l`` is the
*phase function* expansion (the only one a scalar code uses); the other
five carry polarization coupling. They live in the `GreekCoefs` struct at
[`src/Scattering/types.jl:231–244`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/types.jl#L231-L244):

```julia
struct GreekCoefs{FT}
    α::Vector{FT}
    β::Vector{FT}
    γ::Vector{FT}
    δ::Vector{FT}
    ϵ::Vector{FT}
    ζ::Vector{FT}
end
```

The output of `compute_aerosol_optical_properties` packages this with the
scalar bookkeeping ([`src/Scattering/types.jl:281–292`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/types.jl#L281-L292)):

```julia
struct AerosolOptics{FT}
    greek_coefs::GreekCoefs
    ω̃::Union{FT, Array{FT}}      # single-scattering albedo
    k::Union{FT, Array{FT}}      # extinction cross-section
    fᵗ::Union{FT, Array{FT}}     # δ-M truncation factor
end
```

## Phase matrix `Z` from Greek + μ

Once the Greek coefficients are in hand, the per-Fourier-moment phase
matrices ``\mathbf{Z}^{\pm\pm}_m(\mu_i, \mu_j)`` come from
``\Pi^m_l(\mu)``, the generalized spherical-function matrices (Sanghavi
2014 Mie-Fourier Eq. 12–14):

```math
\mathbf{Z}_m(\mu, \mu') = \mathbf{C}^m(\mu, \mu') \pm \mathbf{D}\mathbf{S}^m(\mu, \mu')\mathbf{D},
\qquad
\mathbf{A}^m(\mu, \mu') = \sum_{l=m}^\infty \boldsymbol{\Pi}^m_l(\mu)\,\mathbf{B}_l\,\boldsymbol{\Pi}^m_l(\mu').
```

The coded entry point is
`compute_Z_moments(pol_type, μ, greek, m, arr_type=arr_type)` in
[`src/Scattering/compute_Z_matrices.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/compute_Z_matrices.jl). It returns
``(\mathbf{Z}^{++}_m, \mathbf{Z}^{-+}_m)`` of shape `(NquadN, NquadN)` (or
the pre-broadcast scalar-spectral version that gets expanded across `nSpec`
in `expandOpticalProperties`).

## NAI-2 vs PCW — two ways to compute the Greek expansion

Sanghavi 2014 (the Mie-Fourier paper) compares two decomposition methods:

- **NAI-2** (Numerical Angular Integration, version 2; Siewert's formalism):
  computes the scattering matrix on a Gauss-Legendre quadrature in scattering
  cosine ``\xi``, then projects to Greek coefficients via integrals against
  associated Legendre functions. The ``2N_\mathrm{max} - 1`` quadrature
  roots needed for accurate integration are determined upfront. Default
  path; entry at [`src/Scattering/compute_NAI2.jl:44`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/compute_NAI2.jl#L44).
- **PCW** (Pre-computed Wigner; Domke's formalism corrected): computes
  Greek coefficients *directly* from products of size-averaged Mie
  coefficients ``\langle a_n^* a_m \rangle`` etc. weighted by precomputed
  Wigner-3j symbols. No reconstruction of the angular scattering matrix.
  Entry at [`src/Scattering/compute_PCW.jl:28`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/compute_PCW.jl#L28).

Tradeoffs (Sanghavi 2014 Mie-Fourier §5):

| Regime | Faster method |
|---|---|
| Single particle, ``r_\mathrm{max} > 20\,\mu\mathrm{m}`` (size parameter > 228) | NAI-2 (no angular reconstruction needed) |
| Polydispersion over a finite size range up to ~80 μm | PCW (6–8× faster than NAI-2) |
| Single small particle | PCW (precomputed Wigner-3j makes individual particles cheap) |

NAI-2 is the safer default. PCW pays off for big polydisperse populations
that dominate aerosol retrievals over ocean. Selection is via the
`DECOMP_MAP` in [`src/IO/Parameters.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/IO/Parameters.jl) (`"NAI2"` or `"PCW"`).

## Rayleigh and the Cabannes choice

For the molecular contribution (no aerosols, just air molecules):

- **Pure-elastic runs** (`noRS`) use the *full Rayleigh* Greek expansion. In
  this case rotational Raman is rolled into the *effective* depolarization
  ratio ``\rho_n``; the phase matrix is the standard Rayleigh form with
  depolarization of ~0.028 for air at 532 nm.
- **Raman-aware runs** (`RRS`, `VS_*`) use the *Cabannes* Greek expansion.
  Rotational Raman is then handled *explicitly* by the inelastic kernel
  ([Concepts/08](08_inelastic.md)), and the elastic Cabannes path uses a
  lower depolarization (~0.007).

The selector at [`src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:8–9`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl#L8-L9):

```julia
_rayleigh_greek_source(::Union{noRS, noRS_plus}, greek_rayleigh, greek_cabannes) = greek_rayleigh
_rayleigh_greek_source(::AbstractRamanType, greek_rayleigh, greek_cabannes) = greek_cabannes
```

Mismatch — using Cabannes greek with `noRS`, or full Rayleigh greek with
`RRS` — produces a ~1% bias on Stokes ``I`` and a ~3% bias on ``Q``,
because the polarization-sensitive coefficients (``\beta``, ``\delta``)
shift by the depolarization difference.

## Aerosol microphysics → optics

The size distribution and refractive index are inputs to the Mie
calculation. vSmartMOM ships with two microphysics schemes (in the
`Aerosols` module):

- `TOMAS15Scheme` — TOMAS-15-bin sectional scheme used in CESM/GEOS-Chem.
- `TwoMomentScheme` — two-moment lognormal scheme (median radius + width).

Plus a refractive-index database for common aerosol species (sulfate, dust,
black carbon, organic carbon, sea salt). See [Aerosols (microphysics)](../Aerosols/Overview.md)
for the schemes (status: API still being stabilized — WIP).

## Code anchors

| Concept | Source |
|---|---|
| Mie entry (NAI-2) | [`src/Scattering/compute_NAI2.jl:44`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/compute_NAI2.jl#L44) |
| Mie entry (PCW) | [`src/Scattering/compute_PCW.jl:28`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/compute_PCW.jl#L28) |
| Mie helper (a_n, b_n) | [`src/Scattering/mie_helper_functions.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/mie_helper_functions.jl) |
| GreekCoefs type | [`src/Scattering/types.jl:231–244`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/types.jl#L231-L244) |
| AerosolOptics type | [`src/Scattering/types.jl:281–292`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/types.jl#L281-L292) |
| Polarization types | [`src/Scattering/types.jl:92–143`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/types.jl#L92-L143) |
| Z from Greek + μ | `src/Scattering/compute_Z_matrices.jl::compute_Z_moments` |
| Truncation entry | [`src/Scattering/truncate_phase.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/truncate_phase.jl) (used by `createAero`) |
| Rayleigh / Cabannes selector | [`src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:8–9`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl#L8-L9) |
| GPU Mie kernel | [`src/Scattering/gpu_mie_kernels.jl:30–100`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/gpu_mie_kernels.jl#L30-L100) |
| YAML decomposition map | `src/IO/Parameters.jl::DECOMP_MAP` |

## Hands-on tutorials

Runnable examples with Plotly figures:

- [Scattering & Greek coefficients](../tutorials/Tutorial_Scattering.md)
- [Mie deep dive](../tutorials/Tutorial_MieDeepDive.md)

## References

- **Sanghavi (2014)**, *Revisiting the Fourier expansion of Mie scattering matrices in generalized spherical functions*, JQSRT **136**:16–27, [doi:10.1016/j.jqsrt.2013.12.015](https://doi.org/10.1016/j.jqsrt.2013.12.015). **NAI-2 vs PCW comparison; corrected Domke formalism.**
- Sanghavi et al. (2014), JQSRT **133**:412–433, [doi:10.1016/j.jqsrt.2013.09.004](https://doi.org/10.1016/j.jqsrt.2013.09.004), App. A.
- Hansen, J. E. & Travis, L. D. (1974), *Light scattering in planetary atmospheres*, Space Sci. Rev. **16**:527. (Background on Mie + Greek matrix.)
- Siewert, C. E. (1981, 1982), *On the equation of transfer relevant to the scattering of polarized light*, *On the phase matrix basic to the scattering of polarized light*. (NAI formalism.)
- Domke, H. (1974), *The expansion of scattering matrices for an isotropic medium in generalized spherical functions*, Astrophys. Space Sci. **29**:379. (PCW basis.)
- de Rooij, W. A. & van der Stap, C. C. A. H. (1984), *Expansion of Mie scattering matrices in generalized spherical functions*, A&A **131**:237. (Domke corrections that PCW further fixes.)
- Crib sheet: `docs/dev_notes/theory_references.md` §G.
