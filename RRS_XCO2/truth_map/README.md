# 64-scene high-resolution truth map

Retrieval-facing acceptance criteria and the complete truth/forward/
linearized consistency record are maintained in
[`../inversion/retrieval_setup/TRUTH_FORWARD_LINEARIZED_CONSISTENCY.md`](../inversion/retrieval_setup/TRUTH_FORWARD_LINEARIZED_CONSISTENCY.md).

`true_states.dat` is the authoritative whitespace-separated state table. Its
nesting order is surface, aerosol, SIF, then XCO2 (XCO2 changes fastest).
For a compact description of the repeated surface and aerosol definitions,
see `scene_components.dat`; this avoids decoding those definitions from all
64 rows of the state table.

The SIF-on truth case uses the version-2 campaign definition
`2pi * L_lambda(760 nm) = 0.5 mW m^-2 nm^-1`. Thus every isotropic upwelling
BOA stream has
`L_lambda(760 nm) = 0.5/(2pi) mW m^-2 sr^-1 nm^-1`. This is an unweighted
upward-solid-angle integral at one reference wavelength, not a
wavelength-integrated SIF area. Files and state tables carrying the former
`total_0p5` label or `sif_total` field are superseded inputs and must not be
mixed with this campaign.

The aerosol-on case has total AOD760 0.28. At 550 nm, its species AODs retain
the requested 8:1.8:0.2 ratio: 0.3665052845 sulfate, 0.0824636890 organic
carbon, and 0.0091626321 stratospheric sulfate (total AOD550 0.4581316056).
The values were normalized using the species-specific vSmartMOM Mie optics
and the same endpoint interpolation used by the production O2 solve grid.
The aerosol-off case retains the
same microphysics and vertical profiles with all optical depths set to zero.

All current scenes use a 1000 hPa surface, a 30 degree solar zenith angle,
nadir viewing, and 16 reduced atmospheric layers. The p/T/q source is the current
`RRS_XCO2/config/oco_grass_3aerosol.yaml` profile. Its warm lower troposphere
and relatively moist near-surface specific humidity are closer to a
midlatitude-summer than a midlatitude-winter atmosphere. CO2 is vertically
uniform at the tabulated VMR.

O2 absorption in the A band uses the ABSCO v5.2 `o2_v52_v2.jld2` table and
O2 VMR 0.21. Because ABSCO v5.2 has no A-band H2O table, the A-band truth also
uses the rebuilt HITRAN H2O LUT driven by the common q profile. The retrieval
forward model deliberately uses that same combination. The weak and strong
bands use their corresponding CO2 and H2O ABSCO v5.2 tables. Every active
scene records the exact table paths and profile-preparation provenance.

Earlier instrument and retrieval products were moved to
`../obsolete/pre_absco_closure_20260829/` and must not be used as active
inputs. The archived high-resolution A-band truth does use ABSCO v5.2 O2 with
O2 VMR 0.21, but its aerosol calculation predates the exact retained-core
grid construction described below. Closure testing found the resulting
one-ULP Float32 node displacement to be spectrally significant. The currently
copied A-band arrays were retained only for diagnosis and are being replaced
by `scripts/regenerate_o2_preserve_co2.jl`. That producer writes only the
three A-band variables and refuses to run unless the regenerated weak and
strong CO2 arrays are already finite and complete; it therefore cannot erase
or silently replace the accepted CO2 products.

The planned no-SIF land-glint extensions are isolated under
[`land_glint/`](land_glint/README.md). They use SZA=VZA at the five native
nine-stream directions from 10.2373 through 60 degrees, with relative azimuth
0 degrees, preserving the same atmospheric, surface, aerosol, spectral, and
numerical configuration as this nadir map.

### Chunked aerosol/Raman production

For the current A-band-only correction, run
`scripts/regenerate_o2_preserve_co2.jl`. Aerosol-on O2 A-band RRS scenes cannot
be solved over the complete Raman-shouldered spectrum in one allocation. The
producer splits
the O2 output band into core intervals, adds the full Raman shoulder to both
sides of every interval, solves the expanded grid, and writes only the core
into the existing NetCDF variables. It groups the four XCO2 values because
CO2 does not absorb in this A-band configuration. The general full-map
producer remains `scripts/generate_truth_map_aerosol_chunked.jl`.

Only Cabannes and RRS use the shoulder-expanded solve. The separately stored
Rayleigh/noRS spectrum is monochromatic and is therefore solved directly on
the retained core, with the same canonical full-band aerosol and surface
anchors. This avoids repeating all Fourier/layer operations for thousands of
discarded shoulder samples. The shoulder-invariance closure test bounds the
core-only versus shoulder-carried noRS difference at relative L2 `1.4e-6`.

The retained O2 core is inserted verbatim between independently constructed
shoulders. Do not reconstruct the core as one Float32 range beginning at the
left shoulder: the changed range origin displaced narrow-band nodes by up to
about `9.8e-4 cm-1`. The shared `RRSXCO2Common.raman_solve_grid` constructor
now owns this rule. Retrieval noRS calculations use the same canonical core
directly and therefore do not evaluate the solve-only Raman shoulders. A
regression asserts `rrs_solve_grid[keep] == noRS_grid` bit-for-bit for both
Float32 and Float64. Checkpoints encode this grid-construction version, the
ABSCO version, and the requested shoulder width, and reject older runs.

The runner reuses results across physically identical states: O2 is keyed by
surface and SIF (CO2 has no A-band absorption), while each CO2 band is keyed by
surface and XCO2 (SIF is absent). Results are written incrementally and a JLD2
checkpoint is updated atomically after every physical-state/chunk unit.

```bash
CUDA_VISIBLE_DEVICES=1 CUDA_DEVICE=0 \
AEROSOL_CASE_FILTER=aerosol O2_CHUNK_POINTS=64 FORCE=1 \
TRUTH_OUT="$PWD/RRS_XCO2/truth_map/aerosol_chunked" \
julia --project=. RRS_XCO2/scripts/regenerate_o2_preserve_co2.jl
```

The A-band-only producer never recreates scene files and never writes either
CO2 variable. `FORCE=1` clears only its own atomic JLD2 checkpoint. A scene is
marked incomplete before its first A-band write and complete only after all
three A-band components have finite values.

The corrected SIF-on campaign is built outside the canonical tree and must
pass its all-or-nothing release gate before use:

```bash
julia --project=. RRS_XCO2/scripts/validate_publish_sif_truth_restart.jl
```

That default command does not publish data. It validates all 32 SIF-on scenes,
the version-2 normalization metadata, exact four-XCO2 A-band reuse, and
bit-identical preservation of both CO2 bands against the immutable legacy
archive and matching no-SIF scenes. It writes a checksum receipt inside
`.sif_v2_restart/` only after every check passes. Publication additionally
requires `SIF_RELEASE_ACTION=publish` and the explicit confirmation token
documented at the top of the script. The 32 scene files are promoted from
verified same-directory sidecars before `true_states.dat`; the corrected table
is the final commit point. A caught error rolls all files back from the
archive. Because a POSIX filesystem cannot atomically replace 33 independent
files, a process kill or power loss during those renames can leave a mixed set
of scene files under the still-legacy table. The persistent
`.sif_v2_publication_in_progress` marker makes that state fail-visible; do not
run downstream processing while it exists, and use only the script's explicit
audited resume path.

Superseded and invalid simulation products are isolated under
`../obsolete/`; see `../obsolete/README.md` for their provenance. In
particular, the original 12-layer dataset with the altitude-profile coordinate
bug is under `../obsolete/truth_map/sims_12layer_uncorrected/`. The
`sims_12layer/` directory here is reserved for a future corrected 12-layer
dataset and is intentionally empty until that production run is completed.

`sim_wavelength.nc` stores the reported 0.1 cm-1 grids for O2 A, weak CO2, and
strong CO2. O2 A is solved with an additional 234 cm-1 Raman shoulder on both
sides; those solve-only points are discarded from `hiressim_NNN.nc`. The
strong-CO2 grid contains 995 points: its unchanged 987-point base spectrum plus
eight short-wavelength samples needed for the synthetic Gaussian instrument
convolution. These appended samples span approximately 2041.7062--2041.9980
nm and are integrated directly into every scene and both wavelength files by
`scripts/merge_convolution_shoulders.jl`.

Each scene file stores the TOA I,Q,U Rayleigh/noRS,
Cabannes/RRS-elastic, and rotational-Raman radiances for O2 A, and only the
Rayleigh/noRS radiance for both CO2 bands. No redundant Cabannes+RRS sum is
stored.

### O2 component plots with synthetic OCO overlays

`scripts/plot_truth_state_o2_components.py` generates the four-panel O2
component figures and overlays the matching Mueller-analyzer-processed,
Gaussian-convolved, resampled `I_OCO` values from `OCO_radiances/` as a
thicker black line. For visualization only, the black curve is divided by
the analyzer's unpolarized throughput `M11 = 0.5`. This normalization is not
applied to the stored measurement vector and must not be applied to its
measurement covariance. The high-resolution panels use radiance per cm-1,
whereas the instrument products are stored per nm; the plotting script also
applies the inverse spectral-density Jacobian at the OCO sample centers before
overlaying them, so both curves share valid units.

The regenerated figures also show the matching OCO-2 detector-noise standard
deviation as a gray ±1σ cloud. Absolute Rayleigh panels center the cloud on
the corrected OCO radiance. Component and effect-minus-reference panels center
it on zero: Rayleigh uses the corrected covariance, while Cabannes/RRS panels
use the uncorrected covariance. For aerosol-minus-clear and SIF-on-minus-off,
the effect-present truth spectrum is the noisy measurement and the reference
simulation is deterministic, so the cloud uses only the effect-present scene's
noise. Each component panel reports the maximum `|signal|/σ` and the fraction
of OCO samples exceeding 1σ. The M11 and spectral-density conversions are
applied to both the plotted radiance and plotted σ only; stored measurements
and covariance products remain unchanged. These annotations are per-sample
diagnostics; a spectrally coherent signal can accumulate significance across
many samples even when every individual point lies within the cloud.

```bash
python3 RRS_XCO2/scripts/plot_truth_state_o2_components.py \
  --scene-root RRS_XCO2/truth_map/aerosol_chunked
python3 RRS_XCO2/scripts/plot_truth_state_o2_components.py \
  --scene-root RRS_XCO2/truth_map
```

The first command covers the 32 aerosol scenes; the second covers the 32
no-aerosol scenes stored at the top level. `sims_12layer/` remains empty and
is not populated with these 16-layer products.

For the Rayleigh/noRS spectrum across all three bands, generate one combined
figure per state with:

```bash
python3 RRS_XCO2/scripts/plot_truth_scene_rayleigh_three_bands.py
```

This discovers both clear and aerosol scene roots and writes all 64 plots to
`truth_map/rayleigh_three_bands/`. The same plot-only `M11 = 0.5`
normalization and per-nm to per-cm-1 conversion are applied to the black OCO
overlay; the stored measurement vectors are unchanged.

### Synthetic OCO noise covariances

`OCO_radiances/noise_covariances/OCO2noise_NNN.nc` contains the corrected and
uncorrected diagonal measurement-noise covariance for every state. The
covariances use OCO-2 L1B ATBD Eq. (3-8), Table 3-5 MaxMS values, Table 3-6
MinMS range checks, and representative wavelength-dependent L1B `snr_coef`
arrays. They use the raw stored analyzer measurements and are not divided by
M11. See `OCO_radiances/noise_covariances/README.md` for the exact units,
matrix representation, coefficient provenance, validation, and reproduction
commands.

Matched XCO2 sensitivity plots for the weak and strong CO2 bands can be made
directly from the processed OCO radiances and noise covariances with:

```bash
python3 RRS_XCO2/scripts/plot_xco2_scene_differences.py \
  --surface urban --aerosol none --sif off \
  --xco2-low 400 --xco2-high 420
```

The plotted difference is the high-XCO2 scene minus the otherwise identical
low-XCO2 scene. Following the other truth-map effect plots, the high-XCO2
scene is treated as one noisy measurement and the low-XCO2 simulation as a
deterministic reference; the zero-centered cloud therefore uses only the
high-XCO2 scene's 1-sigma detector noise. Outputs are saved under
`OCO_radiances/xco2_differences/` by default.

Raman is evaluated with the package's supported embedded-solar production
path. Cabannes and noRS use the same quadrature/source geometry so their
difference is physically comparable. The source solar spectrum is the
high-resolution `solar.out` used by the OCORaman workflows; its Fraunhofer
transmission multiplies the 5777 K Planck irradiance.

The driver is resumable and skips completed scene files by default:

```bash
CUDA_DEVICE=1 julia --project=. RRS_XCO2/scripts/generate_truth_map.jl
```

Set `FIRST_STATE`, `LAST_STATE`, or `FORCE=1` to select or overwrite states.
O2 results are cached over XCO2 because CO2 is absent from that band; CO2-band
results are cached over SIF because SIF is absent from those bands.
