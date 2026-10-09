# Central RRS/XCO2 visualizations

This directory is the user-facing home for visualization entry points that
span more than one retrieval campaign.  Generated PNG/DAT products remain in
the data-bearing campaign directories; they are not written beside these
scripts.

For a durable, ready-to-copy reference that explains every option and the
scientific meaning of all six indices, see
[`CHEATSHEET.md`](CHEATSHEET.md).

## Physical retrieval ensembles

For synchronized browsing of rounds 3–6, use the local
[RRS Ensemble Viewer for VS Code](vscode-ensemble-viewer/README.md).
Run **RRS: Open Ensemble Comparison** from the Command Palette in the
Curry-connected window. Its top row shows aerosol/surface figures and its
bottom row shows CO2/SIF figures, each ordered round 3 through round 6.
Left/right arrows cycle the 20 four-surface groups together; click a figure
to enlarge it and press Escape to restore the grid.

Run [`plot_physical_retrieval_ensembles.py`](plot_physical_retrieval_ensembles.py)
with one required retrieval-type index:

| Index | Retrieval interpretation | Local input root | Default plot root |
|---:|---|---|---|
| 1 | Full-column CO2 | `RRS_XCO2/inversion/{corrected,uncorrected}` | `RRS_XCO2/inversion/physical_ensemble_visualizations/full_column/` |
| 2 | Tightly constrained bottom-layer CO2 (`acos_mapped`) | `bottom_layer_XCO2_retrievals/retrievals/` | that root's `plots/physical_ensembles/` |
| 3 | Loosely constrained bottom-layer CO2 (`acos_mapped_tapered_vertical_correlation`); both SIF parameters are retrieved | Selected by `--SIF=0` or `--SIF=1` | the selected campaign's `plots/physical_ensembles/` |
| 4 | The same loose CO2 covariance with round-4 known SIF at 759 nm; no SIF intercept co-retrieval | Selected by `--SIF=0` or `--SIF=1` | the selected campaign's `plots/physical_ensembles/` |
| 5 | Round 5: both SIF coefficients fixed, tighter UTLS aerosol prior, unchanged loose CO2 covariance | Selected by `--SIF=0` or `--SIF=1` | `round5_fixed_sif/plots_nosif/physical_ensembles/` or `round5_fixed_sif/plots_sif/physical_ensembles/` |
| 6 | Round 6: both SIF coefficients fixed, original UTLS aerosol prior, unchanged loose CO2 covariance | Selected by `--SIF=0` or `--SIF=1` | `round6_fixed_sif/plots_nosif/physical_ensembles/` or `round6_fixed_sif/plots_sif/physical_ensembles/` |

For type 4 no-SIF truth, both SIF coordinates are fixed to zero and absent
from the active state.  In the later SIF-on type-4 campaign, the 759-nm value
is known, `SIF760` is derived, and only the spectral slope remains active.
Use `--SIF=0` for no-SIF truth (the default) and `--SIF=1` for corrected
SIF-on truth.  Types 3 through 6 have protected private SIF-on campaign routes;
types 1 and 2 reject `--SIF=1` because no corrected-SIF route is registered.
The private SIF-on NetCDF inputs stay outside the repository, while their PNG
outputs are made visible under these local campaign directories:

- type 3: `bottom_layer_XCO2_retrievals/retrievals_acos_mapped_tapered_vertical_correlation_sif/plots/physical_ensembles/`;
- type 4: `bottom_layer_XCO2_retrievals/round4_known_sif759/plots_sif/physical_ensembles/`;
- type 5: `bottom_layer_XCO2_retrievals/round5_fixed_sif/plots_sif/physical_ensembles/`;
- type 6: `bottom_layer_XCO2_retrievals/round6_fixed_sif/plots_sif/physical_ensembles/`.

These locations are protected by the repository's generated-plot ignore rules.

Round 5 excludes both SIF coefficients from the active state. Its SIF panels
reconstruct the fixed 759-nm anchor and slope, derive the 760-nm radiance, and
label these as fixed quantities without retrieval spread. The plotted SIF
curve is the affine retrieval-state representation, not the nonlinear truth
template. The source priors must be the round-5 `*_tight_utls_*` priors.

Round 6 uses the same fixed-SIF boundary but the original UTLS aerosol prior.
Use index `6` and the `*_standard_utls_*` priors. Its PNGs and reproducible
commands are linked in the [Round 6 plot index](../bottom_layer_XCO2_retrievals/round6_fixed_sif/ENSEMBLE_PLOTS.md).
The [Round 3 versus Round 6 comparison](plot_round3_vs_round6_xco2_psurf.py)
accepts `--SIF=0`, `--SIF=1`, or `--SIF=both` and writes separately named
four-panel PNGs and DAT tables under `comparisons/round3_vs_round6_*`.

The index locks the input root, truth table, aerosol catalog, atmospheric
profile, and figure title.  This prevents a bottom-layer file from being
combined accidentally with a full-column truth table.  The dispatcher also
checks a completed retrieval's campaign geometry, source prior, and saved
state schema before plotting.

### Category controls

For type 1 select the uniform column with `--xco2-ppm`; for types 2--6 select
the injected layer-16 VMR with `--bottom-co2-ppm`.

Use one of:

```text
--aerosol-category no-aerosol
--aerosol-category with-aerosol
--aerosol-category both
```

Use `--product aerosol-surface`, `--product co2-sif`, or `--product all`.
`--product all` writes two figures per selected aerosol category, so it writes
four when `--aerosol-category both` is used.  Each figure retains the 2-by-2
urban/rural/desert/forest comparison.

### Examples

From the repository root, type 1 (full column):

```bash
python3 RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 1 \
  --xco2-ppm 400 \
  --aerosol-category no-aerosol \
  --product all \
  --show-noiseless
```

Type 2 (tightly constrained bottom layer):

```bash
python3 RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 2 \
  --bottom-co2-ppm 400 \
  --aerosol-category with-aerosol \
  --product all \
  --show-noiseless
```

Type 3 (loosely constrained bottom layer, SIF coordinates co-retrieved even
though this local truth category has no SIF):

```bash
python3 RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 3 \
  --SIF=0 \
  --bottom-co2-ppm 400 \
  --aerosol-category with-aerosol \
  --product all \
  --show-noiseless
```

Type 4 (same loose CO2 prior, round-4 reduced SIF state):

```bash
python3 RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 4 \
  --SIF=0 \
  --bottom-co2-ppm 400 \
  --aerosol-category with-aerosol \
  --product all \
  --show-noiseless
```

All five bottom-layer CO2 cases for one type and aerosol category:

```bash
for bottom_ppm in 360 380 400 420 440; do
  python3 RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 3 \
    --bottom-co2-ppm "$bottom_ppm" \
    --aerosol-category with-aerosol \
    --product all \
    --show-noiseless
done
```

For the imported round-3 SIF-on campaign, change only the flag:

```bash
python3 RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 3 \
  --SIF=1 \
  --bottom-co2-ppm 420 \
  --aerosol-category both \
  --product all \
  --show-noiseless
```

Use `INDEX --backend-help` to list styling and diagnostic options such as
`--include-fit-failures`, `--hide-individual`, and `--scale-mode`.

Round 5, all 40 ensemble figures (both SIF cases, five bottom-layer CO2
levels, both aerosol categories, and both products):

```bash
for sif in 0 1; do
  for bottom_ppm in 360 380 400 420 440; do
    python3 RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 5 \
      --SIF="$sif" --bottom-co2-ppm "$bottom_ppm" \
      --aerosol-category both --product all --show-noiseless
  done
done
```

## Physical retrieval ensemble tables

[`tabulate_physical_retrieval_ensembles.py`](tabulate_physical_retrieval_ensembles.py)
writes the same ensembles as Markdown tables.  It takes the same six
retrieval-type indices, is routed and validated by the same registry, and
reuses the plotting backend's loading, pairing, and physical-reconstruction
functions, so a table cannot disagree with the figure drawn from the same
files.

Tables go to a `tables*` directory beside the corresponding `plots*` figure
directory:

| Index | Default table root |
|---:|---|
| 1 | `RRS_XCO2/inversion/physical_ensemble_tables/full_column/` |
| 2 | `bottom_layer_XCO2_retrievals/retrievals/tables/physical_ensembles/` |
| 3 | that campaign's `tables/physical_ensembles/` |
| 4 | `round4_known_sif759/tables_nosif/physical_ensembles/` |
| 5 | `round5_fixed_sif/tables_nosif/physical_ensembles/` or `round5_fixed_sif/tables_sif/physical_ensembles/` |
| 6 | `round6_fixed_sif/tables_nosif/physical_ensembles/` or `round6_fixed_sif/tables_sif/physical_ensembles/` |

Unlike the plotter, the default is a sweep over every truth CO2 case of the
selected campaign; naming `--bottom-co2-ppm` or `--xco2-ppm` narrows it to one
case.  Cases with no paired retrieval yet are reported as pending and skipped,
so this is safe to re-run as a campaign fills in:

```bash
python3 RRS_XCO2/visualization/tabulate_physical_retrieval_ensembles.py 3
```

Each file is named after the corresponding figure base name and contains, for
the four surfaces of one CO2 case and aerosol category:

1. the XCO2 ensemble summary with truth, both classes, and their biases;
2. the dry-air-weighted partial-column contributions of layers 1--4, 5--15, and
   16, which sum to XCO2 and therefore show where a column bias is produced;
3. light-path diagnostics (total AOD, surface pressure, CO2 column, XCO2 bias);
4. retrieved aerosol AOD and median height per species;
5. surface Legendre coefficients per band;
6. the O2 A-band SIF state;
7. the noiseless perturbation-11 retrieval, excluded from all statistics; and
8. ensemble availability, including missing pairs and fit-quality failures.

An `index.md` in the same directory is rebuilt from whichever tables are
present.  `INDEX --backend-help` lists the remaining options, including
`--include-fit-failures` and `--include-empty`.

## Round-3 versus round-4 XCO2 and pressure errors

[`plot_round3_vs_round4_xco2_psurf.py`](plot_round3_vs_round4_xco2_psurf.py)
creates the dedicated 2-by-2 no-SIF comparison.  Round 3 occupies the left
column and round 4 the right; surface pressure is on top and column-averaged
XCO2 below.  Every matched complete pair is retained.  Color identifies the
surface, filled/open markers identify no-aerosol/aerosol scenes, circles show
noisy perturbations, squares show perturbation 11, and an x overlay marks a
pair whose uncorrected retrieval did not converge.  Each populated panel also
reports the corrected and uncorrected mean bias and sample standard deviation
over every plotted pair, including marked failures.

```bash
python3 RRS_XCO2/visualization/plot_round3_vs_round4_xco2_psurf.py
```

Use `--SIF=1` for the corresponding SIF-on comparison:

```bash
python3 RRS_XCO2/visualization/plot_round3_vs_round4_xco2_psurf.py --SIF=1
```

Use `--SIF=both` to combine both truth categories in the same four panels:

```bash
python3 RRS_XCO2/visualization/plot_round3_vs_round4_xco2_psurf.py --SIF=both
```

The complete combined plot has **880 matched pairs per round**: 800 noisy
pairs and 80 noiseless pairs across all surfaces, aerosol cases, and five
bottom-layer CO2 levels. SIF-on markers have a small central `+`; no-SIF
markers are plain. All existing surface colors and circle/square/open/filled
conventions remain unchanged. Summaries pool both SIF categories and include
noiseless pairs and marked failures, as in the separate-category figures.
The PNG and companion DAT are saved under
`bottom_layer_XCO2_retrievals/comparisons/round3_vs_round4_combined_sif/`;
the DAT includes a final `sif_category` column. Separate-category outputs
are not overwritten.

The SIF-on round-4 panels remain empty until all 440 matched pairs are locally
available.  Their placeholders report the current complete-pair count, so a
partial transfer cannot be mistaken for the final round-4 ensemble.

To inspect an explicitly partial round-4 transfer, opt in with:

```bash
python3 RRS_XCO2/visualization/plot_round3_vs_round4_xco2_psurf.py \
    --SIF=1 --allow-partial-round4
```

## Round-3 versus round-5 XCO2 and pressure errors

The round-5 entry point retains the same four-panel layout and marker
conventions, using the 28-coordinate fixed-SIF retrieval schema. Both SIF759
and mSIF are prescribed in round 5: the loader validates the saved parameter
mapping and fixed coefficients, reconstructs absolute SIF coordinates, and
keeps their tangent increments zero. No retrieval products are modified.

```bash
python3 RRS_XCO2/visualization/plot_round3_vs_round5_xco2_psurf.py --SIF=0
python3 RRS_XCO2/visualization/plot_round3_vs_round5_xco2_psurf.py --SIF=1
python3 RRS_XCO2/visualization/plot_round3_vs_round5_xco2_psurf.py --SIF=both
```

Round-5 no-SIF inputs come from `round5_fixed_sif/retrievals_nosif` under
`bottom_layer_XCO2_retrievals`. SIF-on inputs come from
`$RRS_XCO2_PRIVATE_RESULTS_ROOT/bottom_layer_round5_fixed_sif_on_tight_utls_acos_mapped_tapered_vertical_correlation_v1/retrievals`
(the default private-results root is `$HOME/RRS_XCO2_private/results`).
Round-3 inputs are unchanged. Complete figures use 440 matched pairs per
round per SIF category, or 880 when combined, including noiseless members
and marked convergence failures.

PNG figures and companion DAT tables are saved in
`bottom_layer_XCO2_retrievals/comparisons/round3_vs_round5_{nosif,sif,combined_sif}/`.
Round-3/4 outputs are not overwritten. For an explicitly partial round-5
SIF-on snapshot, add `--allow-partial-round5`; otherwise incomplete SIF-on
panels remain placeholders. The original entry point also accepts
`--compare-round=5`, while keeping round 4 as its default.
