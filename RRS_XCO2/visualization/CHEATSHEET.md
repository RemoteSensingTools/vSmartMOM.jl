# Physical retrieval ensemble plotting cheatsheet

This is the practical reference for
[`plot_physical_retrieval_ensembles.py`](plot_physical_retrieval_ensembles.py).
It is intentionally self-contained: begin here even if the details of the
retrieval campaigns are no longer fresh in memory.

## The shortest useful recipe

Run commands from the repository root:

```bash
./RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py INDEX \
  --SIF=0_OR_1 \
  CO2_SELECTION \
  --aerosol-category AEROSOL_SELECTION \
  --product PRODUCT \
  --show-noiseless
```

The launcher is executable.  This equivalent form is useful on a copy whose
executable permission was not preserved:

```bash
python3 RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py \
  INDEX [OPTIONS]
```

Choose the four capitalized items as follows:

1. Choose `INDEX` from the retrieval-type table below.
2. Choose `--SIF=0` for no-SIF truth or `--SIF=1` for corrected SIF-on
   truth.  SIF-on is currently registered for indices 3 through 6.
3. For index 1, use `--xco2-ppm VALUE`.  For indices 2--6, use
   `--bottom-co2-ppm VALUE`.
4. Choose `AEROSOL_SELECTION` as `no-aerosol`, `with-aerosol`, or `both`.
5. Choose `PRODUCT` as `aerosol-surface`, `co2-sif`, or `all`.

The script prints the selected input root, a validation file, the output
directory, and a state/availability inventory before printing each PNG path.

If only the index is supplied, the defaults are: 400 ppm CO2, no-SIF truth,
both aerosol categories, both physical products, local panel scales,
individual perturbations 01--10 shown, perturbation 11 hidden, fit failures
excluded, no paired-difference diagnostic, and 180 DPI.  This normally writes
four PNGs.

## Step 1: choose the retrieval-type index

The index is not merely a label.  It selects the retrieval files, truth table,
prior family, state-vector interpretation, plot title, and default output
directory as one protected unit.

| Index | Meaning | Important distinction |
|---:|---|---|
| `1` | Full-column CO2 retrievals | The truth CO2 VMR changes uniformly through the column.  This is the original full-column campaign. |
| `2` | Tightly constrained bottom-layer CO2 retrievals | Only truth layer 16 changes, but the original ACOS-mapped prior correlates the lower layers relatively tightly.  The legacy 30-coordinate state retrieves both SIF parameters. |
| `3` | Loosely constrained bottom-layer CO2 retrievals | Uses the tapered vertical CO2 correlation: lower layers can move more independently, while correlation tightens aloft.  The legacy 30-coordinate state still retrieves `SIF760` and `mSIF`, even for no-SIF truth. |
| `4` | Loose bottom-layer CO2 with the reduced round-4 SIF state | Uses the same tapered CO2 prior as type 3.  For the currently registered no-SIF campaign, SIF at 759 nm and its slope are fixed to zero and excluded from the 28-coordinate active state. |
| `5` | Round-5 fixed SIF with tighter UTLS aerosol priors | Both SIF coefficients are fixed (zero for no-SIF); UTLS log-AOD and log-height prior standard deviations are reduced tenfold. The CO2 prior is unchanged from types 3 and 4. |
| `6` | Round-6 fixed SIF with the original UTLS aerosol priors | Same fixed-SIF state as round 5, but all non-SIF priors, including UTLS and CO2, match the round-3/4 source prior. |

Round-4 SIF-on retrievals have a slightly different 29-coordinate state:
SIF at 759 nm is known, `mSIF` remains active, and `SIF760` is derived.  That
campaign is registered separately with its corrected-v2 SIF truth table via
`--SIF=1`. It is not inferred from the no-SIF type-4 root.

Round 5 uses 28 active coordinates for both truth cases. Its SIF plots show
the fixed affine state anchored at 759 nm, not an independently retrieved
760-nm radiance or slope. No SIF ensemble spread is expected.

### Default input and output locations

| Index | Retrieval input root | Default figure directory |
|---:|---|---|
| `1` | `RRS_XCO2/inversion/{corrected,uncorrected}/` | `RRS_XCO2/inversion/physical_ensemble_visualizations/full_column/` |
| `2` | `RRS_XCO2/bottom_layer_XCO2_retrievals/retrievals/` | `RRS_XCO2/bottom_layer_XCO2_retrievals/retrievals/plots/physical_ensembles/` |
| `3`, `--SIF=0` | `RRS_XCO2/bottom_layer_XCO2_retrievals/retrievals_acos_mapped_tapered_vertical_correlation_nosif/` | that retrieval root's `plots/physical_ensembles/` |
| `3`, `--SIF=1` | private round-3 SIF retrieval results | `RRS_XCO2/bottom_layer_XCO2_retrievals/retrievals_acos_mapped_tapered_vertical_correlation_sif/plots/physical_ensembles/` |
| `4`, `--SIF=0` | `RRS_XCO2/bottom_layer_XCO2_retrievals/round4_known_sif759/retrievals_nosif/` | `RRS_XCO2/bottom_layer_XCO2_retrievals/round4_known_sif759/plots_nosif/physical_ensembles/` |
| `4`, `--SIF=1` | private round-4 SIF retrieval results | `RRS_XCO2/bottom_layer_XCO2_retrievals/round4_known_sif759/plots_sif/physical_ensembles/` |
| `5`, `--SIF=0` | `RRS_XCO2/bottom_layer_XCO2_retrievals/round5_fixed_sif/retrievals_nosif/` | `RRS_XCO2/bottom_layer_XCO2_retrievals/round5_fixed_sif/plots_nosif/physical_ensembles/` |
| `5`, `--SIF=1` | private round-5 fixed-SIF/tight-UTLS retrieval results | `RRS_XCO2/bottom_layer_XCO2_retrievals/round5_fixed_sif/plots_sif/physical_ensembles/` |
| `6`, `--SIF=0` | `RRS_XCO2/bottom_layer_XCO2_retrievals/round6_fixed_sif/retrievals_nosif/` | `RRS_XCO2/bottom_layer_XCO2_retrievals/round6_fixed_sif/plots_nosif/physical_ensembles/` |
| `6`, `--SIF=1` | private round-6 fixed-SIF/standard-UTLS retrieval results | `RRS_XCO2/bottom_layer_XCO2_retrievals/round6_fixed_sif/plots_sif/physical_ensembles/` |

The directories are separate, so identical-looking filenames from different
retrieval types cannot overwrite one another.

SIF-on retrieval NetCDF files remain in the private results tree.  Only the
generated PNG products are written to the repository-facing campaign folders
listed above, where they are covered by the repository's PNG/`plots/` ignore
rules and therefore cannot be added to Git accidentally.

## Step 2: choose the CO2 case

### Index 1: full-column truth

Use:

```text
--xco2-ppm 380
--xco2-ppm 400
--xco2-ppm 420
--xco2-ppm 440
```

This value is the vertically uniform truth VMR and therefore the nominal
column XCO2.  If omitted, it defaults to 400 ppm.

Example:

```bash
./RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 1 \
  --xco2-ppm 400 \
  --aerosol-category both \
  --product all \
  --show-noiseless
```

### Indices 2--6: bottom-layer truth

Use:

```text
--bottom-co2-ppm 360
--bottom-co2-ppm 380
--bottom-co2-ppm 400
--bottom-co2-ppm 420
--bottom-co2-ppm 440
```

This is the truth VMR in layer 16 only.  Layers 1--15 remain at 400 ppm, so
the number is not the resulting column XCO2.  The approximate mappings are:

| Bottom-layer VMR | Column XCO2 |
|---:|---:|
| 360 ppm | 397.5081 ppm |
| 380 ppm | 398.7540 ppm |
| 400 ppm | 400.0000 ppm |
| 420 ppm | 401.2460 ppm |
| 440 ppm | 402.4919 ppm |

If omitted, the bottom-layer selector defaults to 400 ppm.

Do not mix the selectors: index 1 rejects `--bottom-co2-ppm`, and indices 2--4
reject `--xco2-ppm`.  This guard prevents a plausible-looking plot of the
wrong physical experiment.

## Step 3: choose the aerosol category

Use exactly one of:

| Option | Meaning | Number of physical categories plotted |
|---|---|---:|
| `--aerosol-category no-aerosol` | Select truth scenes with zero aerosol optical depth. | 1 |
| `--aerosol-category with-aerosol` | Select truth scenes containing the three-aerosol mixture with total AOD at 760 nm of 0.28. | 1 |
| `--aerosol-category both` | Generate the no-aerosol and with-aerosol categories separately. | 2 |

The default is `both`.  “Both” does not combine the two aerosol cases in one
panel; it creates separate PNG files.

## Step 4: choose the physical product

| Option | Figure contents |
|---|---|
| `--product aerosol-surface` | Aerosol vertical AOD profiles above surface-reflectance spectra.  Each figure contains the urban, rural, desert, and forest surfaces in a 2-by-2 layout. |
| `--product co2-sif` | CO2 molecular-number-density profiles above O2 A-band SIF spectra, again in a 2-by-2 surface layout. |
| `--product all` | Generate both figures for every selected aerosol category. |

The default is `all`.  Therefore:

- one aerosol category plus `all` produces two PNGs;
- `both` aerosol categories plus one product produces two PNGs;
- `both` plus `all` produces four PNGs.

`--co2-only` is the older spelling of `--product co2-sif`.  Prefer the latter
because it is clearer.  Do not combine `--co2-only` with an explicit
non-default `--product`; the script rejects that ambiguous combination.

## Optional display and ensemble controls

### `--show-noiseless`

Overlay perturbation 11, the retrieval of the unperturbed measurement, using
the noiseless line/marker style.  Perturbation 11 is displayed only: it is
never included in the ensemble mean, standard deviation, or percentile cloud.
Both its corrected and uncorrected files must be valid for the paired overlay
to appear.  Because it is visible data, it can expand automatically chosen
axis limits.

This option is recommended for most diagnostic figures.

### `--include-fit-failures`

By default, a noisy corrected/uncorrected pair contributes to the plotted
ensemble only when both retrieval files are:

- complete;
- state-step converged; and
- accepted by the saved spectral fit-quality gate.

`--include-fit-failures` keeps complete, state-step-converged members even if
the fit-quality gate failed.  Affected cards are marked.  Use this when
studying compensating or pathological retrieval behavior; do not interpret it
as changing a failed spectral fit into a successful retrieval.

### `--hide-individual`

Suppress the faint curves for perturbations 01--10.  The truth, corrected and
uncorrected means, and 16th--84th percentile uncertainty clouds remain.  This
is useful for a cleaner manuscript figure; omit it when diagnosing individual
retrieval behavior.

### `--scale-mode local` or `--scale-mode shared`

- `local` is the default.  Each surface card is zoomed independently, making
  small corrected/uncorrected differences easier to see.
- `shared` uses common limits for comparable cards.  This makes absolute
  cross-surface differences honest and immediate, but close curves can become
  visually indistinguishable.  Limits are shared within one aerosol-category
  figure, not between the separately written aerosol and no-aerosol figures.

### `--include-paired-correction`

In addition to the standard aerosol/surface plot, create a diagnostic showing
the paired corrected-minus-uncorrected physical curves and their 16th--84th
percentile range.  It affects only runs that include the `aerosol-surface`
product.  It does not add a second CO2/SIF figure.

### `--dpi INTEGER`

Set PNG resolution in dots per inch.  The default is `180`.  Higher values
increase pixel dimensions and file size but do not change any scientific
calculation.

### `--output-dir PATH`

Override the campaign-specific output directory.  Usually omit this option so
the index keeps figures beside the appropriate retrieval campaign.  It is
useful for temporary tests, for example:

```bash
./RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 3 \
  --bottom-co2-ppm 400 \
  --aerosol-category with-aerosol \
  --product co2-sif \
  --output-dir /tmp/round3_plot_test
```

The output directory is created automatically when necessary.  A relative
path is interpreted relative to the shell's current directory.  An existing
same-named PNG is overwritten without prompting, so do not send several
retrieval types to one custom directory unless their names are separated.

## How the figure layout avoids obscuring data

Legend and annotation placement is automatic; there is no command-line
position setting to remember.  The plotter reserves space according to the
kind of information being shown:

- curve/statistic keys and aerosol-species line styles occupy two separate
  footer rows below the aerosol/surface panels;
- CO2/SIF state keys and ensemble-diagnostic keys likewise occupy separate
  footer rows;
- the campaign identifier has its own small header line instead of being
  appended to an already long scientific title;
- each CO2 card reports its XCO2 summary above the plotting rectangle; and
- each SIF card places its numerical summary in a dedicated strip beneath the
  wavelength axis.

The small completion badge and the pressure/column summary remain inside a
CO2 panel because their corners are consistently clear of the profiles.  If a
new campaign puts data in those corners, adjust those two card-level
annotations in the backend rather than moving the shared figure legends back
onto the data axes.

Typography follows named `FONT_*` constants near the top of the backend.  The
smallest routine metadata is kept above the former fine-print size, while
compact spectral panels, main axes, headers, and legends use progressively
larger sizes.  Adjust those constants together if a future output medium
requires another global increase.

## SIF selection

Use one simple campaign flag:

```text
--SIF=0
--SIF=1
```

`--SIF=0` is the default and selects no-SIF truth. `--SIF=1` selects the
corrected-v2 SIF-on campaign and truth table.  It is registered for types 3
and 4; types 1 and 2 reject it because they do not have a protected
corrected-SIF campaign route.  The lower-level `--sif-case` spelling is locked
out so SIF units cannot be mixed accidentally.

## Ready-to-copy commands

### Full-column, 400 ppm, both aerosol categories, all products

```bash
./RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 1 \
  --xco2-ppm 400 \
  --aerosol-category both \
  --product all \
  --show-noiseless
```

### Tight bottom layer, 400 ppm, aerosol scenes only

```bash
./RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 2 \
  --bottom-co2-ppm 400 \
  --aerosol-category with-aerosol \
  --product all \
  --show-noiseless
```

### Loose bottom layer, 380 ppm, no-aerosol CO2/SIF plot

```bash
./RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 3 \
  --SIF=0 \
  --bottom-co2-ppm 380 \
  --aerosol-category no-aerosol \
  --product co2-sif \
  --show-noiseless
```

### Loose bottom layer, 420 ppm, corrected SIF-on truth

```bash
./RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 3 \
  --SIF=1 \
  --bottom-co2-ppm 420 \
  --aerosol-category both \
  --product all \
  --show-noiseless
```

### Round 4, 420 ppm, aerosol/surface plot including fit failures

```bash
./RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 4 \
  --bottom-co2-ppm 420 \
  --aerosol-category with-aerosol \
  --product aerosol-surface \
  --show-noiseless \
  --include-fit-failures
```

### All five bottom-layer CO2 cases for types 3 and 4

```bash
for retrieval_type in 3 4; do
  for bottom_ppm in 360 380 400 420 440; do
    ./RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py \
      "$retrieval_type" \
      --bottom-co2-ppm "$bottom_ppm" \
      --aerosol-category both \
      --product all \
      --show-noiseless
  done
done
```

This last recipe can produce up to 40 PNGs: two retrieval types times five CO2
cases times two aerosol categories times two physical products.

## How to interpret the figures

- Black represents truth.
- Red/orange represents the uncorrected retrieval ensemble.
- Green/blue represents the corrected retrieval ensemble.  CO2 uses the more
  strongly separated orange/blue pair because those curves often overlap.
- Faint individual curves are matched noisy perturbations 01--10.
- Shaded bounds are the 16th--84th percentile of the reconstructed physical
  curves, not independently varied one-parameter error bars.
- Perturbation 11 is the noiseless retrieval and is excluded from statistics.
- Empty or explicitly provisional cards mean files are missing or too few
  valid paired retrievals are currently available; the script does not fill
  them using a different state.

“Corrected” and “uncorrected” refer to the two retrieval measurement classes
used throughout this experiment.  Both retrievals use the linearized
`RS_type::noRS` forward model.  The corrected class treats the Rayleigh truth
simulation as its measurement, whereas the uncorrected class fits the
Cabannes+RRS truth simulation without representing RRS in its retrieval
forward model.  Their perturbations are matched by index before ensemble
statistics are formed.

Every PNG contains all four surfaces in the fixed 2-by-2
urban/rural/desert/forest layout.  There is currently no one-surface-only
option.

## Output filenames

For full-column type 1, filenames begin with:

```text
xco2_<ppm>_nosif_<aerosol-category>_
```

For bottom-layer types 2--4, they begin with:

```text
bottom_co2_<ppm>_nosif_<aerosol-category>_
```

The standard suffixes are:

| Suffix | Product |
|---|---|
| `profiles_surface.png` | Aerosol vertical profiles and surface spectra |
| `co2_profiles.png` | CO2 concentration profiles and SIF spectra |
| `paired_correction_effect.png` | Optional corrected-minus-uncorrected aerosol/surface diagnostic |

The campaign output directories, rather than the filename alone, distinguish
retrieval types 2, 3, and 4.

## Help and troubleshooting

Central routing help:

```bash
./RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py --help
```

All plotting-backend options:

```bash
./RRS_XCO2/visualization/plot_physical_retrieval_ensembles.py 3 \
  --backend-help
```

Common messages:

- **Use `--xco2-ppm`, not `--bottom-co2-ppm`**: index 1 was selected with the
  bottom-layer selector.
- **Use `--bottom-co2-ppm`, not `--xco2-ppm`**: index 2, 3, or 4 was selected
  with the full-column selector.
- **No completed retrieval was found**: the selected campaign has not produced
  any usable NetCDF retrieval yet, or the files are not in the registered
  campaign root.
- **Missing required input / wrong prior / wrong state model**: do not bypass
  the guard with manual truth or retrieval paths.  First determine whether the
  intended campaign was copied into the correct directory.
- **A card says missing or shows fewer than 10/10 pairs**: inspect the inventory
  printed by the command.  A pair is omitted when either its corrected or
  uncorrected member is unavailable or fails the active validity gates.
- **No unique state is found**: the requested CO2 value is not one of the
  discrete truth cases, or it is not yet available for the selected campaign.
- **SIF-on input is missing**: verify that the private corrected-v2 campaign
  and its matching truth table have been transferred to Curry.

The central launcher intentionally forbids overriding `--inversion-root`,
`--truth-table`, `--scene-components`, `--vertical-profile-table`,
`--atmospheric-profile`, and `--campaign-label`.  Those values define the
scientific identity of the numeric index.  Allowing them to be changed
independently would defeat the main protection provided by this interface.

## Markdown tables of the same ensembles

[`tabulate_physical_retrieval_ensembles.py`](tabulate_physical_retrieval_ensembles.py)
is the tabular companion.  It takes the same retrieval-type index and the same
campaign guards, and reuses this plotter's backend for loading, pairing, and
physical reconstruction, so the tables and the figures cannot disagree.

```bash
./RRS_XCO2/visualization/tabulate_physical_retrieval_ensembles.py INDEX
```

That default sweeps every truth CO2 case of the campaign and writes one
Markdown file per CO2 case and aerosol category, plus an `index.md`, into the
`tables*` directory beside the index's `plots*` directory.  Cases without a
paired retrieval yet are listed as pending and skipped, so re-running as the
campaign fills in only refreshes what has data.

| Index | Default table root |
|---:|---|
| `1` | `RRS_XCO2/inversion/physical_ensemble_tables/full_column/` |
| `2` | `RRS_XCO2/bottom_layer_XCO2_retrievals/retrievals/tables/physical_ensembles/` |
| `3` | that retrieval root's `tables/physical_ensembles/` |
| `4` | `RRS_XCO2/bottom_layer_XCO2_retrievals/round4_known_sif759/tables_nosif/physical_ensembles/` |

Narrow the sweep to one case with the index's own CO2 selector, exactly as for
the figures:

```bash
./RRS_XCO2/visualization/tabulate_physical_retrieval_ensembles.py 3 \
  --bottom-co2-ppm 420 \
  --aerosol-category with-aerosol
```

Every file reports, for all four surfaces: the XCO2 ensemble summary; the
dry-air-weighted partial-column contributions of layers 1--4, 5--15, and 16,
which sum to XCO2 and locate where a column bias is produced; light-path
diagnostics (total AOD, surface pressure, CO2 column, XCO2 bias); aerosol AOD
and median height per species; surface Legendre coefficients; the SIF state;
the noiseless perturbation-11 retrieval; and ensemble availability.  Ensemble
statistics use the same paired, converged, fit-quality-accepted members as the
figures, and perturbation 11 is never included in them.

`INDEX --backend-help` lists the remaining options, including
`--include-fit-failures`, `--include-empty`, and `--no-index`.

## Four-panel Round-3 versus Round-4 bias comparison

From the repository root, combine no-SIF and SIF-on retrievals with:

```bash
python3 RRS_XCO2/visualization/plot_round3_vs_round4_xco2_psurf.py --SIF=both
```

This separate comparison script uses `--SIF=0` for no-SIF only, `--SIF=1`
for SIF-on only, and `--SIF=both` for their union. (`both` is not an option
of the physical-ensemble script.) Round 3 is left, Round 4 right; pressure
is above column-averaged XCO2. Each point has corrected-minus-truth on x
and uncorrected-minus-truth on y.

- Surface colors: urban blue, rural mauve, desert orange, forest green.
- Open markers: aerosols; filled markers: no aerosol.
- Circles: perturbations 01–10; smaller translucent squares: noiseless 11.
- In the combined figure only, a small internal `+` means SIF-on; plain
  markers mean no-SIF.
- Red `x`: the uncorrected retrieval did not converge.
- Summaries: pooled mean bias and sample standard deviation over all
  displayed pairs, including noiseless pairs and marked failures.

With complete data there are 880 pairs in each round (440 per SIF category).
The combined PNG and its exact-value DAT are written to
`RRS_XCO2/bottom_layer_XCO2_retrievals/comparisons/round3_vs_round4_combined_sif/`.
The separate SIF-on/no-SIF figures are preserved. `--allow-partial-round4`
explicitly permits plotting an incomplete Round-4 SIF-on transfer;
otherwise the Round-4 panels await the complete dataset.
