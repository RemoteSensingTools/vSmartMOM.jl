# vSmartMOM principles

This page summarizes the core ideas behind vSmartMOM (vector Smart Matrix-Operator Method), consolidated from the primary references.

## What problem vSmartMOM solves

- Vector (polarized) radiative transfer in layered media with molecular (Rayleigh) and aerosol scattering, plus gaseous absorption.
- Efficient, accurate multiple scattering using the matrix-operator method (doubling/adding) with discrete-ordinate quadrature.
- Modular treatment of polarization states (I, IQ, IQU, IQUV) and Fourier decomposition in azimuth.
- Optional linearization with respect to aerosol optical properties to facilitate Jacobians and sensitivity studies.

## Method at a glance

- Discretize angle with a quadrature (e.g., Gauss full-sphere, hemisphere, Radau) to get a finite set of streams.
- Expand the vector phase matrix into Fourier series over azimuth and into Legendre polynomials over polar angle; use a truncation strategy.
- Build per-layer single-scattering properties (optical thickness τ, single-scattering albedo ϖ, and polarized phase matrices Z⁺⁺, Z⁻⁺).
- Compose layers via matrix-operator doubling/adding into a composite medium, preserving polarization.
- Couple to surfaces via BRDF models (Lambertian, RPV, Ross-Li, Cox-Munk, Canopy) in a vector-consistent way.  The Cox-Munk ocean surface provides full polarization support via Mueller-matrix Fresnel reflection from wind-roughened wave facets; the Canopy surface couples a multi-layer vegetation canopy with a soil BRDF via internal adding-doubling.

For equation-level details and code mapping of `elemental`, `doubling`, and `interaction`, see `Core RT Theory (Doubling/Adding)`.

## Truncation strategies and accuracy

- δ-m and δ-fit approaches are used to remove/approximate strong forward peaks in scattering, improving convergence at low order.
- vSmartMOM applies these truncations consistently to the vector phase matrix, following Sanghavi & Stephens (2015).
- Practical guidance:
  - Choose L_trunc consistent with aerosol phase function complexity and desired accuracy.
  - Use Radau or Gauss full-sphere quadrature depending on hemispheric symmetry and source configuration.

## Polarization handling

- Supports Stokes vector subsets via fixed types: I; IQ; IQU; IQUV.
- Phase matrix moments and the discrete-ordinate system are built per polarization type; this ensures type-stable kernels.

## Linearization (aerosol Jacobians)

- The matrix-operator formalism enables analytic derivatives with respect to aerosol properties (e.g., τ_ref, phase function parameters) by differentiating the composed operators.
- vSmartMOM exposes aerosol optical property construction (Mie-based) and truncation, enabling efficient sensitivity calculations.

## Layer optics and composition

- Per-layer quantities include:
  - τ_λ (with gaseous absorption), ϖ_λ, τ, ϖ, and Fourier-expanded Z matrices.
  - Rayleigh properties are computed from depolarization factor and wavelength; aerosol properties from Mie models at reference λ.
- Composition:
  - Added layer (A) and composite layer (C) operators are combined according to the standard matrix-operator algebra.
  - Doubling logic computes appropriate optical thickness sub-steps for numerical stability.

## Surface interaction and HDRF

- Surfaces are modeled via BRDFs; interaction terms are included in the boundary condition and in HDRF post-processing.
- **Scalar surfaces** (Lambertian, RPV, Ross-Li) populate only the (1,1) Stokes block of the surface reflectance matrix.
- **Polarized surfaces** (Cox-Munk) fill the full Mueller-matrix blocks (I-Q, U-V coupling), requiring numerical azimuthal integration with polarized Fourier kernels and a TMS single-scattering correction for the specular sun-glint peak.
- **Canopy surfaces** internally solve canopy sub-layers via adding-doubling before presenting an effective reflectance to the atmospheric RT.
- The HDRF and VZA postprocessing utilities provide azimuthal averaging and RAMI-style outputs.

## Practical setup tips

- Accuracy knobs: quadrature choice, L_trunc, max Fourier moment m, and depolarization factor for Rayleigh.
- Stability: small per-step optical depths via doubling; ensure consistent units (wavenumber bands in cm⁻¹).
- Performance: choose architecture CPU/GPU; GPU accelerates batched operations and kernel calls.

## References

- Sanghavi, S., Davis, A.B., Eldering, A. (2014). vSmartMOM: A vector matrix-operator method-based radiative transfer model linearized with respect to aerosol properties. JQSRT 133: 412–433.
- Sanghavi, S., Stephens, G. (2015). Adaptation of the δ-m and δ-fit truncation methods to vector radiative transfer: Effect of truncation on radiative transfer accuracy. JQSRT 159: 53–68.
- Grant, I.P., Hunt, G.E. (1969). Discrete space theory of radiative transfer I. Proc. Roy. Soc. A 313(1513): 183–197.
- Plass, G.N., Kattawar, G.W., Catchings, F.E. (1973). Matrix operator theory of radiative transfer. 1: Rayleigh scattering. Applied Optics 12(2): 314–329.
