# PRIVATE TODO — Species-dependent Rayleigh and Raman expansion

Status: private implementation handoff for the Exo-Earth pressure/composition
study. Do not present these items as existing production capability.

## Current capability audit

- Molecular Raman constant constructors currently exist for N2, O2, and H2 in
  `src/Inelastic/src/molecular_constructors.jl`.
- Atmospheric Raman mixture calculations in `raman_atmo_prop.jl` and related
  helpers are substantially specialized to N2/O2.
- H2 therefore has partial molecular/stellar machinery but is not yet a
  general member of an arbitrary atmospheric mixture.
- CO2, H2O, and CH4 molecular constructors/constants are absent.
- Existing effective Cabannes/Rayleigh mixture behavior must be preserved and
  validated while generalizing dispatch.

## 1. Species-dependent elastic Rayleigh scattering

- [ ] Define a common molecule interface for refractive index or molecular
  polarizability, King/depolarization factor, Rayleigh cross section, and
  wavelength-validity metadata.
- [ ] Retain and validate N2 and O2 implementations.
- [ ] Promote H2 from its partial inelastic implementation into the general
  atmospheric Rayleigh mixture path.
- [ ] Add reviewed molecular data and constructors for CO2, H2O, and CH4.
- [ ] Implement wavelength-dependent sigma_Rayl(lambda) rather than assuming a
  universal lambda^-4 scaling when dispersion is material over 250--450 nm.
- [ ] Compute mixture extinction and polarization quantities by number-column
  weighting, not by averaging species VMRs or depolarization ratios directly.
- [ ] Return both total Rayleigh and elastic Cabannes cross sections so the
  elastic source and Raman redistribution conserve scattering consistently.
- [ ] Expose per-species contributions for ExoOptics provenance and Jacobians.

## 2. Rotational Raman scattering (RRS)

- [ ] Generalize RRS mixture loops from hard-coded N2/O2 arguments to a typed
  collection of participating molecules.
- [ ] Validate existing N2, O2, and H2 rotational constants, nuclear-spin
  statistical weights, energy levels, Raman shifts, polarizability tensors,
  and temperature-dependent populations against primary references.
- [ ] Add CO2 rotational Raman constants and selection rules.
- [ ] Add H2O rotational Raman support appropriate to an asymmetric top; do
  not force it through the linear-molecule J-to-J+/-2 implementation.
- [ ] Add CH4 rotational Raman support with spherical-top degeneracies and
  selection rules; do not approximate it as N2/O2.
- [ ] Decide whether water-vapor and CH4 RRS are required for the first
  350--450 nm science run based on quantitative cross-section tests, but keep
  the mixture API capable of representing them.
- [ ] Test Raman-grid padding so all transitions that feed 350--450 nm are
  present even when their incident wavelengths lie outside the reported band.

## 3. Water vibrational Raman scattering (VRS)

- [ ] Add H2O vibrational modes with documented frequencies, degeneracies,
  temperature populations, Raman activities/polarizability derivatives, and
  polarization properties.
- [ ] Include at minimum the bending and OH-stretch manifolds needed to assess
  redistribution into/out of 350--450 nm; document any mode aggregation.
- [ ] Support rotational-vibrational structure only where justified by the
  spectral resolution and available molecular data.
- [ ] Treat both Stokes and anti-Stokes transitions and enforce detailed-
  balance/population consistency.
- [ ] Verify photon/energy bookkeeping and the wavelength-versus-wavenumber
  Jacobian conventions used by the inelastic source terms.

## 4. Pressure/composition integration

- [ ] Provide an API that accepts layer T, pressure boundaries, and arbitrary
  species VMRs and returns species/mixture Rayleigh, Cabannes, RRS, and VRS
  optical properties.
- [ ] Keep spectral arrays in the package's internal wavenumber convention and
  preserve the spectral third dimension used by batched RT.
- [ ] Ensure changing surface pressure/composition updates both elastic and
  inelastic optical properties through `BatchContext`/`update_model!` or the
  appropriate current model-update boundary.
- [ ] Add analytic or upstream-AD derivatives with respect to surface pressure
  and composition after the forward implementation is validated.
- [ ] Extend gas-abundance Jacobians through CIA and MT_CKD, including the H2O
  self-continuum and H2O line self-broadening dependence on abundance.
- [ ] Extend the analytic surface-pressure Jacobian through pressure-dependent
  line shapes; its CIA/MT_CKD midpoint-pressure and column factors are already
  included, but ordinary line cross sections are currently held fixed.

## 5. 350--450 nm validation campaign

- [ ] Use native 0.01 nm absorption sampling over the 250--450 nm shortwave
      windows; retain 0.01 cm^-1 for the O2 A-band and NIR line windows before
      convolution.
- [ ] Compare integrated observables over 360--400 and 400--440 nm from the
  same simultaneous RT call.
- [ ] Establish elastic-only, absorption-only, RRS, and RRS+H2O-VRS cases.
- [ ] Validate N2/O2 against current Raman regression cases before adding one
  species at a time: H2, CO2, H2O, CH4.
- [ ] Check optical-depth closure:
  total molecular scattering = Cabannes + sum(all Raman channels).
- [ ] Check spectral redistribution closure on a grid padded beyond the output
  window and verify convergence versus padding and spectral step.
- [ ] Test Stokes sign/rotation conventions using
  `docs/src/pages/conventions.md` before any external-code comparison.
- [ ] Add CPU reference tests first, then CUDA/Metal coverage supported by the
  inelastic implementation.

## Likely source locations

- `src/Inelastic/src/raman_constants.jl`
- `src/Inelastic/src/molecular_constructors.jl`
- `src/Inelastic/src/inelastic_cross_section.jl`
- `src/Inelastic/raman_atmo_prop.jl`
- `src/Inelastic/raman_stellar_prop.jl`
- `src/Inelastic/inelastic_helper.jl`
- `src/Inelastic/stellar_inelastic_helper.jl`
- `src/Inelastic/types.jl`

## ExoOptics consumer

The experiment-level pressure, composition, spectral-window, LUT, and disk-
integration work is tracked in ExoOptics:

`TODO/PRIVATE_EXO_EARTH_PRESSURE_RAMAN_TODO.md`
