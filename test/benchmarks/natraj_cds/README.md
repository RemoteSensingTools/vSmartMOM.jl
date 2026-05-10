# Natraj CDS Rayleigh Tables

This directory stores the Coulson-Dave-Sekera Rayleigh reference tables used by
Natraj, Li, and Yung (2009):

V. Natraj, K.-F. Li, and Y. L. Yung, "Rayleigh Scattering in Planetary
Atmospheres: Corrected Tables Through Accurate Computation of X and Y
Functions", ApJ 691, 1909-1920. DOI: 10.1088/0004-637X/691/2/1909.

Source archive:
`https://web.gps.caltech.edu/~vijay/Rayleigh_Scattering_Tables/CDS/CDS.tar.gz`

The extracted `CDS/` files cover:

- `tau`: 0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1
- Lambertian surface albedo: 0.0, 0.25, 0.8
- `mu0`: 0.1, 0.2, 0.4, 0.6, 0.8, 0.92, 1.0
- `mu`: 0.02, 0.06, 0.1, 0.16, 0.2, 0.28, 0.32, 0.4, 0.52,
  0.64, 0.72, 0.84, 0.92, 0.96, 0.98, 1.0
- relative azimuth: 0, 30, 60, 90, 120, 150, 180 degrees
- Stokes `I`, `Q`, `U`, for both upwelling TOA and downwelling BOA radiation

Manual benchmark:

```bash
cd test
julia --project=. benchmarks/benchmark_natraj_cds.jl --quick
julia --project=. benchmarks/benchmark_natraj_cds.jl
```
