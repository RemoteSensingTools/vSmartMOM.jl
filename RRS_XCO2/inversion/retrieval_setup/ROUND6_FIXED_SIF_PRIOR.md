# Round 6: fixed SIF, original UTLS prior

Approved definition (2026-09-10): fix both SIF parameters, but do not tighten
or fix the UTLS aerosol parameters. This supersedes the earlier prospective
use of the name "round 6" for a CO2-covariance experiment.

Round 6 has the same 28-coordinate state and fixed SIF boundary as round 5:
SIF-off has `Lnu759=mSIF=0`; SIF-on has `Lnu759=0.004818031987713776`
and `mSIF=1.2291230681458325e-5`, in native radiance/wavenumber units.
`SIF760 = Lnu759 + mSIF*(nu760-nu759)` is injected into the core model.
Both SIF columns are excluded from the active solve.

| Prior coordinate | Round 5 standard deviation | Round 6 standard deviation |
|---|---:|---:|
| UTLS `ln(AOD760)` | 0.075 | 0.75 |
| UTLS `ln(z0/km)` | 0.01 | 0.10 |

The source is the approved round-3/round-4 tapered prior, SHA-256
`34c0e81b7a853af157b68bb879db0771342a9579460835675f92df8aab7f9375`.
All non-SIF means `xa[1:32,:]` and all non-SIF covariances
`Sa[1:32,1:32,:]` must be exactly identical to that source, including the
entire CO2 block, aerosol cross-covariances, pressure and surface priors.
Only the two SIF rows/columns are zeroed; the four upper CO2 layers remain
fixed as before. The active/full map is `[1; 6:32]`, and the active/core map
is `1:28`. No CO2 decorrelation or performance optimization is introduced.

The forward model, OE implementation, analytical Jacobians, spectroscopy,
solar input, nine streams, Fourier stopping policy, grids, detector operator,
truth spectra, and noise realizations are reused from the accelerated round-5
checkpoint `7acab57000dae259207a6760faae156cfa1734f6`. New round-6 workflow
files are deployed as a separately checksummed overlay; the base checkout is
read-only to this workflow. NetCDF output provenance identifies both the base
checkpoint and the combined base/overlay/input campaign identity.

Local no-SIF partitions: Curry physical GPU 1 owns the 20 no-aerosol states
`1:5,21:25,41:45,61:65`; Wurst physical GPU 1 owns the 20 aerosol states
`11:15,31:35,51:55,71:75`. Each runs 11 noise members and both measurement
classes (440 files per worker). The first noiseless corrected/uncorrected
pair must pass convergence and fit checks. Device visibility is restricted
to physical GPU 1, which appears as CUDA device 0 inside each process.

Gattaca2 owns all 40 SIF-on states. Its submitter schedules a two-retrieval
noiseless smoke job for state 018, then a dependent full array, at most two
GPU tasks concurrently. Production is released only after smoke success.
The existing corrected-SIF truth release barrier is retained.

See [GATTACA_ROUND6_FIXED_SIF.md](../GATTACA_ROUND6_FIXED_SIF.md) for transfer,
submission, and monitoring commands. Round-3/4/5 products are not overwritten.
