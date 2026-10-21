# Gattaca2 round-4 SIF-on retrievals

This path retrieves all 40 bottom-layer SIF-on states with SIF radiance known
exactly at 759 nm and only the SIF spectral slope active. It reads the
published corrected-v2 truth, OCO-like measurements, and frozen noise. It
never regenerates or modifies truth.

Round 4 may be queued while round 3 is running, but **must use a second,
non-nested Git checkout**. The existing checkout remains the immutable source
of data-bearing round-3 and full-column inputs. Never switch, pull, or detach
the checkout from which the active round-3 Slurm array was launched.

The private output namespace is

```text
$HOME/RRS_XCO2_private/results/
  bottom_layer_round4_known_sif759_sif_on_acos_mapped_tapered_vertical_correlation_v1/
```

The smoke test is state 018, perturbation 11. Production contains one state
per Slurm array task (`0-39%2`), runs perturbation 11 before 1-10, and keeps
corrected/uncorrected pairs on the same task.

## 1. Create the isolated source checkout

Run these commands on the Gattaca2 login node only after an approved round-4
commit has been pushed. Replace the angle-bracket value with that exact
40-character commit SHA.

```bash
set -euo pipefail

export ROUND3_REPO="$HOME/code/uni_vSmartMOM"
export RRS_REPO="$HOME/code/uni_vSmartMOM_round4_known_sif759"
export ROUND4_CODE_CHECKPOINT_SHA='<approved-40-character-commit-sha>'

test -d "$ROUND3_REPO/.git"
test ! -e "$RRS_REPO"

git clone --single-branch --branch suniti_multi_sensor \
    https://github.com/RemoteSensingTools/vSmartMOM.jl.git \
    "$RRS_REPO"
git -C "$RRS_REPO" switch --detach "$ROUND4_CODE_CHECKPOINT_SHA"

# Manifest.toml is deliberately ignored by Git.  Reuse the exact manifest
# that ran round 3; the round-4 readiness identity hashes this copy.
test -f "$ROUND3_REPO/Manifest.toml"
cp -p "$ROUND3_REPO/Manifest.toml" "$RRS_REPO/Manifest.toml"
test "$(sha256sum "$RRS_REPO/Manifest.toml" | awk '{print $1}')" = \
     "$(sha256sum "$ROUND3_REPO/Manifest.toml" | awk '{print $1}')"

test "$(git -C "$RRS_REPO" rev-parse HEAD)" = \
    "$ROUND4_CODE_CHECKPOINT_SHA"
test -z "$(git -C "$RRS_REPO" status --porcelain --untracked-files=all)"
```

The two repositories must be siblings (or otherwise disjoint), not nested.

## 2. Bind round 4 to the existing data checkout

Keep these exports in the shell used for submission:

```bash
export ROUND3_REPO="$HOME/code/uni_vSmartMOM"
export RRS_REPO="$HOME/code/uni_vSmartMOM_round4_known_sif759"
export RRS_PRIVATE_ROOT="$HOME/RRS_XCO2_private"
export BOTTOM_XCO2_CAMPAIGN_ROOT="$ROUND3_REPO/RRS_XCO2/bottom_layer_XCO2_retrievals"
export FULL_COLUMN_TRUTH_ROOT="$ROUND3_REPO/RRS_XCO2/truth_map"
```

The launcher derives and validates three additional ignored/private runtime
inputs from that same data checkout:

```text
$ROUND3_REPO/RRS_XCO2/inversion/instrument/representative_stokes_coefficients.nc
$BOTTOM_XCO2_CAMPAIGN_ROOT/truth/scene_components.dat
$ROUND3_REPO/src/SIF_emission/sif-spectra.csv
```

Their byte hashes are included in the immutable round-4 input-set identity.
The fresh source checkout therefore needs no copied truth or ignored data.

## 3. Small private prior inputs

After the round-4 code is committed, transfer these three files into the new
private campaign's `retrieval_setup/` directory:

```text
apriori_states_round4_known_sif759_on_acos_mapped_tapered_vertical_correlation.nc
apriori_states_round4_known_sif759_on_acos_mapped_tapered_vertical_correlation.dat
source_apriori_states_acos_mapped_tapered_vertical_correlation.nc
```

The last file is a byte-exact copy of the Curry file named
`apriori_states_acos_mapped_tapered_vertical_correlation.nc`, from which the
reduced round-4 prior was constructed. Rename it with the `source_` prefix
when installing it. It lives in the dedicated round-4 namespace rather than
being read from or written over the completed round-3 campaign. Raw hashes of
independently generated NetCDF priors can differ because their metadata
contains a creation timestamp and an absolute source path, so round 4 pins
this exact copy. Its expected SHA-256 is:

```text
34c0e81b7a853af157b68bb879db0771342a9579460835675f92df8aab7f9375
```

Record the SHA-256 of all three files. No truth, measurement, noise, ABSCO,
solar, corrected-v2 release,
Stokes-coefficient, component-catalog, or SIF-template file needs another
transfer when the existing round-3 installation is retained.

## 4. Submit smoke plus production

Export the site settings and approved hashes:

```bash
export GATTACA_SLURM_ACCOUNT=cmml
export GATTACA_SLURM_QOS=normal
export ROUND4_PRIOR_SHA256='<64-character-NC-sha256>'
export ROUND4_PRIOR_SUMMARY_SHA256='<64-character-DAT-sha256>'
export ROUND4_SOURCE_PRIOR_SHA256='<64-character-source-prior-sha256>'
```

Then submit from the new checkout:

```bash
bash "$RRS_REPO/RRS_XCO2/inversion/submit_gattaca_round4_known_sif759_retrievals.sh"
```

The submitter explicitly copies `RRS_REPO`, `RRS_PRIVATE_ROOT`,
`BOTTOM_XCO2_CAMPAIGN_ROOT`, and `FULL_COLUMN_TRUTH_ROOT` into both Slurm jobs
and records them in `retrieval_jobs.env`. Every task independently verifies:

- the clean detached round-4 checkout and approved checkpoint;
- a distinct canonical data-bearing checkout;
- exact canonical bottom-layer and full-column input roots;
- one scheduler-assigned CUDA device with at least 48 GiB;
- Julia 1.12.5 and external-input checksums;
- all corrected-v2 publication receipts and all 40 SIF-on products;
- the exact 29-dimensional prior and its three static runtime inputs;
- a private output path outside both repositories and both truth roots.

Expected behavior: round 3 continues from `$ROUND3_REPO`; round 4 executes code
from `$RRS_REPO`, reads truth through the two explicit legacy roots, and writes
only beneath the dedicated private round-4 result namespace.
