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

### Effective aerosol Fourier support

The maximum particle size parameter determines a conservative Mie-series
length, but the extreme-radius tail of an integrated size distribution can
carry negligible weight. The optional `radiative_transfer.greek_beta_cutoff`
therefore applies a second test after size-distribution integration:

```math
l_{\beta,a} = \max\{l : |\beta_{a,l}(\nu)| \geq \epsilon_\beta
\text{ for at least one } \nu \text{ in the band}\}.
```

This is evaluated separately for each aerosol `a`; their maximum contributes
to the band-level `m_max` alongside Rayleigh, surface, and source support. The
absolute value is intentional because β can be signed. Only β is tested: the
other five Greek families can be zero or change sign and are not used as the
scalar phase-function support benchmark. The last passing degree is retained,
so the implementation does not assume monotonic decay.

The maximum-size result remains the allocation ceiling, and stream/user caps
are still applied after component aggregation. This criterion shortens the
Fourier loop but does not itself reduce `nstreams` or angular-operator size.

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

### Pure-gas no-redistribution Rayleigh data

New science workflows select the calculation explicitly. Two methods are
available, selected by dispatch on an `AbstractMolecularRayleighMethod`:

| Method | Physics | Species |
|---|---|---|
| `UsingRefInd()` | refractive index plus (effective) King factor | He, Ar, N₂, O₂, CO₂ |
| `UsingPolTen()` | molecular polarizability tensor (mean polarizability ``\bar\alpha``, anisotropy ``\gamma``) from the Raman-module constants | N₂, O₂, H₂ |

```julia
method = UsingRefInd()            # RefractiveIndexRayleigh() is an alias
p = molecular_rayleigh_properties(method, :CO2, wavelength_nm)
ratio = molecular_rayleigh_cross_section_ratio(method, :CO2, λ, λ_ref)
source = molecular_rayleigh_provenance(method, :CO2)

q = molecular_rayleigh_properties(UsingPolTen(), :N2, wavelength_nm;
                                  temperature_K = 288.15)
```

Both return the same named-tuple fields: `species`, `wavelength_nm`,
`refractive_index`, `king_factor`, `depolarization`, and `cross_section_cm2`.
The two-argument forms `molecular_rayleigh_properties(species, λ)` and
`molecular_rayleigh_cross_section_ratio(species, λ, λ_ref)` remain the
refractive-index method. `RefractiveIndexRayleigh` is kept as an alias of
`UsingRefInd`, so existing calling code is unaffected.

#### Refractive index plus King factor (`UsingRefInd`)

The RI helper returns a refractive index, King factor, solver depolarization
ratio, and scattering cross section in cm² molecule⁻¹ for `He`, `Ar`, `N2`,
`O2`, and `CO2`. Its Ar and CO₂ coefficients were reviewed for the Habitable
Worlds Observatory (HWO) exo-Earth studies and carry a machine-readable
provenance record; He, N₂ and O₂ have no such record yet. The method argument is intentionally visible so an RI calculation cannot be
confused with the polarizability-tensor implementation.

The cross-section-ratio helper is useful for a matched-reference-optical-depth
experiment,

```math
\tau_R(\lambda) = \tau_R(\lambda_{ref})
                  \frac{\sigma_R(\lambda)}{\sigma_R(\lambda_{ref})}.
```

The cross section follows the conventional refractive-index/King-factor
expression used by the cited coefficient compilations:

```math
\sigma_R = \frac{24\pi^3\nu^4}{N^2}
           \left(\frac{n^2-1}{n^2+2}\right)^2 F_K,
\qquad
\rho = \frac{6(F_K-1)}{7F_K+3}.
```

Here ``\rho`` is exactly the depolarization convention consumed by
`get_greek_rayleigh`. He and Ar are monatomic, so ``F_K=1`` and ``\rho=0``;
they still have different cross-section dispersion. These helpers deliberately
do not infer pressure, mean molecular mass, or layer column. They also do not
turn on Raman redistribution. A caller comparing gases at fixed
``\tau_R(\lambda_{ref})`` must set the optical-depth profile explicitly and
build the Greek phase source from the returned depolarization at each
wavelength (or place wavelengths with distinct depolarization in separate
bands).

For Ar and CO₂, the coefficient provenance and distinct evidence intervals are machine-readable through
`molecular_rayleigh_provenance`:

- The [Peck--Fisher Ar measurement](https://doi.org/10.1364/JOSA.54.001362)
  spans 17 vacuum wavelengths from 467.9465 to 2058.6909 nm. Sneep and Ubachs
  adopt its dispersion formula over 5000--33000 cm⁻¹
  (303.0303--2000 nm). Ar has the monatomic identity ``F_K=1``; it has no
  empirical King-factor interval.
- The [Bideau-Mehu CO₂ measurement](https://doi.org/10.1016/0030-4018%2873%2990289-7)
  spans 180.7--253.7 nm. The implementation is the corrected/scaled formula
  recorded by [Sneep and Ubachs](https://doi.org/10.1016/j.jqsrt.2004.07.025).
  Its broader 180.7--1694.5 nm envelope is a secondary database-declared
  formula range, not the primary measurement range. The source record is
  frozen at refractiveindex.info database commit
  `37773dce8128ccabcd9d5108d9958d6041244661`, path
  `database/data/main/CO2/nk/Bideau-Mehu.yml`, SHA-256
  `d3b52d2a0b667a2e95aec110f35a7fced0c2ab3f02583d2e604aeb58c204f55c`.
- The CO₂ King fit in Sneep and Ubachs is based on
  [Alms--Burnham--Flygare measurements](https://doi.org/10.1063/1.431821) at
  457.9, 488.0, 514.5, 568.2, and 647.1 nm. Therefore 457.9--647.1 nm is its
  empirical interpolation span; use outside that span is King-factor
  extrapolation, independently of whether the dispersion formula is in range.
- [Wilmouth and Sayres (2019)](https://doi.org/10.1016/j.jqsrt.2019.05.031)
  directly check both gases continuously over 333--363 nm.
  [Thalman et al. (2014)](https://doi.org/10.1016/j.jqsrt.2014.05.030), with
  [its erratum](https://doi.org/10.1016/j.jqsrt.2016.12.014), check Ar in
  selected regions between 350 and 660 nm, with corrected tabulated points at
  370.0, 405.8, 532.2, and 660.0 nm. Those selected regions are not a
  continuous 350--660 nm validation interval.
  [He et al. (2021)](https://doi.org/10.5194/acp-21-14927-2021) provide a
  continuous 307--725 nm direct check for CO₂.

The primary-measurement range, adopted dispersion envelope, King-factor
support, and direct-scattering checks are not interchangeable; dispersion and
King-factor extrapolation should be reported separately.

The N₂ and O₂ RI coefficient sets are those listed in
[He et al. (2021)](https://doi.org/10.5194/acp-21-14927-2021), Table 1, which
gives their dispersion-fit domains as 21 360--39 370 cm⁻¹ (254--468 nm) for N₂
and 18 315--34 722 cm⁻¹ (288--546 nm) for O₂. At longer wavelengths both are
dispersion extrapolations.

This API performs no explicit frequency redistribution. Its conventional
Rayleigh cross section and effective King factor are not a separately computed
Cabannes-only cross section. The RI method neither calls nor modifies the
tensor/Raman machinery in `InelasticScattering`.

#### Polarizability tensor (`UsingPolTen`)

`UsingPolTen()` uses the polarizability tensors that the Raman module already
carries for N₂, O₂, and H₂ (`getMolecularConstants` in
`src/Inelastic/src/molecular_constructors.jl`). It calls the same routines as
the RT model's molecular-depolarization path (`compute_effective_coefficents!`,
`compute_σ_Rayl_coeff!`). No new spectroscopic constants are introduced.

- Mean polarizability: the Buldakov et al. single-pole dispersion and
  temperature dependence,
  ``\bar\alpha(\nu,T) = \bar\alpha_{00}(1+\alpha_b T+\alpha_c T^2)
  /(1-(2\pi c\nu/\omega_0)^2)``, evaluated at ``\nu = 10^7/\lambda``
  cm⁻¹ and `temperature_K` (default 288.15 K, the reference temperature of
  the He, Ar, N₂ and CO₂ RI fits).
- Anisotropy: ``\gamma = \bar\gamma_{00}``, which is **independent of
  wavelength and temperature** in the Raman-module constants.
- Derived quantities:
  ``\rho = 6\gamma^2/(45\bar\alpha^2+7\gamma^2)`` (the same convention as
  the RI method), ``F_K = 1 + 2\gamma^2/(9\bar\alpha^2)``, and
  ``\sigma = \tfrac{128}{3}\pi^5\nu^4\bar\alpha^2F_K``.
- `refractive_index` is the Lorentz--Lorenz equivalent,
  ``(n^2-1)/(n^2+2) = 4\pi N\bar\alpha/3``, at `number_density_cm3`
  (default 2.546899×10¹⁹ cm⁻³). Pass `number_density_cm3 = nothing` to get
  `refractive_index = nothing`.

Other species raise an `ArgumentError` that suggests `UsingRefInd()`.
`molecular_rayleigh_provenance` is currently defined only for the RI method.

At 288.15 K over 350--800 nm the two methods differ, as PolTen/RI − 1, by
the amounts below. Refractivities are compared at the density of each RI fit:
2.546899×10¹⁹ cm⁻³ for N₂ and 2.68678×10¹⁹ cm⁻³ (273.15 K) for O₂. The O₂ RI
`refractive_index` refers to that denser state, so comparing it with the
PolTen default density would show a spurious ≈5 % difference in ``n-1``.
The cross sections, King factors and depolarizations are per molecule and do
not depend on density.

| Quantity | N₂ | O₂ |
|---|---|---|
| ``\sigma`` | −0.51 % (350 nm) to +0.46 % (800 nm) | −1.54 % (350 nm) to +0.05 % |
| ``n-1`` | within ±0.14 % | within ±0.46 % |
| ``F_K`` | −0.26 % to +0.18 % | −1.60 % to +0.98 % |
| ``\rho`` | −7.3 % to +5.2 % | −15 % to +11 % |

The depolarization differences arise mainly because the tensor anisotropy has
no dispersion: ``\rho_\mathrm{PolTen}`` rises with wavelength, while the RI
King factors fall. ``\rho \propto F_K - 1`` turns small King-factor
differences into percent-level ``\rho`` differences.

#### Model defaults

These helpers do **not** change any model default. The automatic
`model_from_parameters` elastic Rayleigh optical depth remains the terrestrial
Bodhaine et al. (1999) air formula, scaled by surface pressure and rescaled by
King factor to the configured depolarization. With the default `depol: -1`,
that depolarization and the Rayleigh/Cabannes Greek coefficients come from
the N₂/O₂ polarizability constants (`getRamanAtmoConstants`), evaluated at the
band-centre wavenumber and the column-weighted mean temperature. Using these
pure-gas helpers inside an RT run requires setting the optical-depth profile
and the Greek source explicitly. Atmospheric composition is not yet connected
to molecular columns and scattering-weighted phase-matrix mixing.

### Elastic Rayleigh vs Cabannes in RT runs

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
- `TwoMomentScheme` — two-moment bulk scheme (AOD + effective radius, with fixed geometric width).

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
