# Bottom-layer retrieval visualizations

The authoritative cross-campaign physical-ensemble entry point now lives at
`RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py`.  Its required
index selects full-column (1), tight bottom-layer (2), loose bottom-layer (3),
or loose bottom-layer with the reduced round-4 SIF state (4), and sends plots
to the corresponding campaign directory automatically.  See
`RRS_XCO2/visualization/README.md` for exact commands.

The compatibility entry points in this directory adapt every Python figure
producer in `RRS_XCO2/inversion/` to the isolated bottom-layer campaign.  They
delegate to the shared implementation rather than copying it, while supplying
the correct bottom-layer retrieval root, truth table, component catalog,
atmospheric profile, and output directories.

The scripts live here, outside `retrievals/`, so code remains trackable while
the private generated retrieval products stay ignored by Git.  Default figures
and tables are written beneath `../retrievals/` and therefore cannot be
published accidentally by an ordinary `git add`.

## Round-4 state compatibility

The state-oriented entry points support both the original 30-coordinate
retrievals and the reduced round-4 state:

- SIF-off round 4 saves 28 active coordinates.  SIF at 759 nm and its slope
  are fixed to zero.
- SIF-on round 4 saves 29 active coordinates.  SIF at 759 nm is known exactly,
  `mSIF` is retrieved, and `SIF760` is derived from those two quantities.

Plots and exported tables reconstruct a common physical-state view and state
clearly which SIF coordinates were retrieved, fixed, or derived.  The
gain/noise overlay applies the tangent map only: the fixed 759-nm intercept is
never added to `G delta_y`.

For SIF-on round 4, the plotted comparison curve is explicitly the
**retrieval-space reference truth**: the exact nonlinear truth-template value
at 759 nm plus the truth slope under the retrieval's affine spectral model.
The template's independently evaluated 760-nm value is a forward-model
mismatch diagnostic, not the `x_truth` coordinate used in state-error plots.

Round-4 SIF-on products must be plotted against the exact corrected SIF truth
table transferred with that campaign.  The local historical bottom-layer
table still contains the obsolete `total_0p5` convention; the plotting code
rejects that pairing instead of silently drawing a physically incorrect truth
curve.

## CO2 convention

Bottom-layer plots distinguish two quantities:

- **bottom-layer CO2** is the injected/retrieved layer-16 VMR and is the
  primary coordinate of this experiment;
- **XCO2** is the corresponding dry-air-column mean.

Truth profiles are reconstructed as 400 ppm in layers 1--15 and the tabulated
`bottom_co2_ppm` in layer 16.  They are never reconstructed as 16 copies of
the much smaller column perturbation in `xco2_ppm`.

The five bottom-layer values occur in this order within every state block:

| Within-block index | Bottom-layer CO2 | Nominal XCO2 |
|---:|---:|---:|
| 1 | 360 ppm | 397.508088 ppm |
| 2 | 380 ppm | 398.754044 ppm |
| 3 | 400 ppm | 400.000000 ppm |
| 4 | 420 ppm | 401.245956 ppm |
| 5 | 440 ppm | 402.491912 ppm |

Every physical-ensemble figure prints only its selected relationship, for
example `bottom-layer CO2=360 ppm -> column XCO2=397.5081 ppm`; it no longer
repeats the complete five-case mapping above every plot.  Each CO2 card also
reports the corrected and uncorrected XCO2 mean and sample standard deviation
from its matched noisy perturbations 01--10.  Perturbation 11 remains excluded
from these statistics even when its noiseless solution is displayed.  The
truth annotation is derived from `truth/true_states.dat`, while retrieved XCO2
is reconstructed from each terminal pressure/profile state and checked against
the saved diagnostic.

CO2 molecular-concentration cards use logarithmic concentration and altitude
axes to resolve the lower atmosphere.  Their 16th--84th percentile clouds are
supplemented by horizontal one-sample-standard-deviation bars at alternating
layer centers; alternating corrected/uncorrected markers reduce overlap.

## Entry points

Run these from this directory with Python 3 or directly as executables:

- `plot_retrieval_fit.py RETRIEVAL.nc`: terminal fit and residual/noise panels.
- `plot_aerosol_trajectory.py RETRIEVAL.nc`: AOD and height LM trajectory.
- `plot_surface_legendre_moments.py RETRIEVAL.nc`: truth/prior/retrieved
  Legendre contributions.
- `plot_retrieval_state_convergence.py RETRIEVAL.nc`: full compact state path,
  including bottom-layer CO2 and XCO2, plus OE diagnostics.
- `plot_compact_state_trajectory.py RETRIEVAL1.nc [RETRIEVAL2.nc ...]`: overlay
  several trajectories for one state and retrieval class.
- `plot_state001_corrected_uncorrected_compact.py [STATE]`: historical entry
  point, now generalized to any state; state 001 remains its default.
- `plot_corrected_vs_uncorrected_errors.py STATE`: one-state corrected versus
  uncorrected displacement plot, including the optional gain/noise prediction.
- `plot_all_corrected_vs_uncorrected_errors.py`: unified plot over every
  currently completed matched pair.
- `plot_physical_retrieval_ensembles.py`: aerosol profiles, surface spectra,
  CO2 molecular-concentration profiles, and SIF spectra.  This bottom-layer
  wrapper has explicit `--campaign round3-nosif` and `--campaign
  round4-nosif` presets (round 4 is the default).  It defaults to the 400 ppm
  bottom-layer control; use `--bottom-co2-ppm 360`, `380`, `420`, or `440` for
  another experiment. It defaults to the SIF-off truth category;
  use `--sif-case angular_integral760_0p5` for the corresponding SIF-on
  retrievals. Here `0.5` means
  `2pi * L_lambda(760 nm) = 0.5 mW m^-2 nm^-1`, not a wavelength-integrated
  SIF area. SIF-on plots use a distinct
  `sif_angular_integral760_0p5` filename component and therefore cannot
  overwrite established `nosif` plots. The plot loader also requires the
  version-2 SIF provenance embedded in every SIF-on retrieval and matches its
  truth provenance against the selected state-table row.  A legacy
  `total_0p5` table is rejected with an explicit error rather than being mixed
  silently with the corrected campaign. `--aerosol-category no-aerosol` or
  `with-aerosol` restricts one invocation to one physical category;
  `--product aerosol-surface` or `co2-sif` restricts it to one of the two
  output figures. The default values `both` and `all` retain both categories
  and products, respectively.
- `plot_spectral_fit_visualization_trials.py`: spectral residual atlas, fit
  card, RRS budget, and fit-quality summary.  It requires all perturbations
  01--11 for at least one corrected/uncorrected state pair.
- `instrument/fit_oco2_eofs_to_retrieval_residuals.py`: EOF fits for all
  completed matched retrieval pairs.
- `instrument/compare_oco2_observed_radiances.py`: 80-scene synthetic-versus-
  OCO-2 radiance-range check.

`export_retrieval_state_trajectory.py RETRIEVAL.nc` is also provided because
it writes the physical tables consumed alongside several figures.

Use `--help` on any entry point to see optional output paths.  Explicit command
line paths override all wrapper defaults.

## Examples

```bash
cd RRS_XCO2/bottom_layer_XCO2_retrievals/visualization

./plot_retrieval_fit.py \
  ../retrievals/corrected/retrieval_state013_perturbation11.nc

./plot_corrected_vs_uncorrected_errors.py 13 --no-gain-prediction

./plot_physical_retrieval_ensembles.py \
  --bottom-co2-ppm 400 --show-noiseless

# Imported SIF-on round-4 data require both explicit transferred paths:
./plot_physical_retrieval_ensembles.py \
  --inversion-root /path/to/round4_sif_on/retrievals \
  --truth-table /path/to/round4_sif_on/truth/true_states.dat \
  --bottom-co2-ppm 400 \
  --sif-case angular_integral760_0p5 \
  --show-noiseless
```

## Individual round-4 no-SIF commands

The following block selects the completed noiseless corrected retrieval for
state 001.  Change `corrected` to `uncorrected`, `001` to another zero-padded
state index, or `11` to a noise perturbation from `01` through `10`.

```bash
cd RRS_XCO2/bottom_layer_XCO2_retrievals/visualization

R4=../round4_known_sif759/retrievals_nosif
OUT=../round4_known_sif759/plots_nosif/individual
FILE="$R4/corrected/retrieval_state001_perturbation11.nc"
mkdir -p "$OUT"
```

One terminal spectral fit with residual/noise panels:

```bash
./plot_retrieval_fit.py "$FILE" \
  --output "$OUT/state001_p11_corrected_fit.png"
```

One complete state trajectory and convergence diagnostic:

```bash
./plot_retrieval_state_convergence.py "$FILE" \
  --output "$OUT/state001_p11_corrected_convergence.png" \
  --table-output "$OUT/state001_p11_corrected_trajectory.dat"
```

One aerosol-only trajectory:

```bash
./plot_aerosol_trajectory.py "$FILE" \
  --output "$OUT/state001_p11_corrected_aerosol.png"
```

One surface-Legendre diagnostic:

```bash
./plot_surface_legendre_moments.py "$FILE" \
  --output "$OUT/state001_p11_corrected_surface.png"
```

One human-readable state export (including fixed/derived SIF coordinates):

```bash
./export_retrieval_state_trajectory.py "$FILE" \
  --output "$OUT/state001_p11_corrected_states.dat" \
  --surface-output "$OUT/state001_p11_corrected_surface.dat"
```

Two trajectories from the same retrieval class (useful for comparing a noisy
realization with its noiseless counterpart):

```bash
./plot_compact_state_trajectory.py \
  "$R4/corrected/retrieval_state001_perturbation01.nc" \
  "$R4/corrected/retrieval_state001_perturbation11.nc" \
  --output "$OUT/state001_p01_p11_corrected_trajectories.png" \
  --table-output "$OUT/state001_p01_p11_corrected_trajectories.dat"
```

All completed perturbations for one truth state:

```bash
./plot_state001_corrected_uncorrected_compact.py 1 \
  --inversion-root "$R4" \
  --output "$OUT/state001_corr_uncorr_summary.png" \
  --table-output "$OUT/state001_corr_uncorr_summary.dat"
```

The 3-by-7 corrected-versus-uncorrected truth-displacement plot, including the
terminal gain prediction by default:

```bash
./plot_corrected_vs_uncorrected_errors.py 1 \
  --inversion-root "$R4" \
  --output "$OUT/state001_corr_uncorr_displacements.png" \
  --table-output "$OUT/state001_corr_uncorr_displacements.dat"
```

Add `--no-gain-prediction` to the preceding command to suppress its open-marker
gain overlay.

The four-surface physical ensemble for the 400-ppm bottom-layer control:

```bash
./plot_physical_retrieval_ensembles.py \
  --inversion-root "$R4" \
  --bottom-co2-ppm 400 \
  --sif-case off \
  --show-noiseless \
  --output-dir ../round4_known_sif759/plots_nosif/physical_ensembles
```

The spectral-fit visualization suite for state 001 (requires a complete
matched perturbation ensemble):

```bash
./plot_spectral_fit_visualization_trials.py \
  --inversion-dir "$R4" \
  --state 1 \
  --output-dir ../round4_known_sif759/plots_nosif/spectral_state001
```

Finally, the unified truth-displacement figure over every completed state is:

```bash
./plot_all_corrected_vs_uncorrected_errors.py \
  --inversion-root "$R4" \
  --output ../round4_known_sif759/plots_nosif/all_completed_displacements.png \
  --table-output ../round4_known_sif759/plots_nosif/all_completed_displacements.dat
```

For imported SIF-on results, replace `R4` with their local retrieval root and
pass their transferred corrected truth table explicitly as
`--truth-table /path/to/corrected/true_states.dat` to every state- or
truth-oriented command.  The fit-only command does not need a truth table.

Ensemble scripts intentionally retain pending placeholders or report that a
complete 01--11 ensemble is unavailable while production retrievals are still
running.  They do not reinterpret partial data as a finished ensemble.
