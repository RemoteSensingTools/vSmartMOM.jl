# Bottom-layer CO2 campaign

This campaign keeps CO2 at 400 ppm in atmospheric layers 1--15 and changes
only layer 16, the lowest layer in vSmartMOM's TOA-to-BOA ordering. Its five
bottom-layer VMRs are 360, 380, 400, 420, and 440 ppm. On the fixed 1000 hPa,
16-layer production profile, layer 16 contains 0.06229780735737046 of the dry
column, so the corresponding XCO2 values are approximately 397.508088,
398.754044, 400.000000, 401.245956, and 402.491912 ppm.

The state nesting is unchanged from the full-column campaign:

```text
surface -> aerosol -> SIF -> CO2 case
```

The SIF-on case is normalized by
`2pi * L_lambda(760 nm) = 0.5 mW m^-2 nm^-1`, so every isotropic upwelling
BOA stream has
`L_lambda(760 nm) = 0.5/(2pi) mW m^-2 sr^-1 nm^-1`. This is an unweighted
upward-solid-angle integral at 760 nm, not a wavelength-integrated SIF area.
Its state-table label is `angular_integral760_0p5`; legacy `total_0p5` inputs
are invalid for this campaign.

The only change is that the innermost dimension now has five cases rather
than four, producing 80 states. `truth/true_states.dat` is the authoritative
table after it is generated. `truth/control_reuse_map.dat` records how the 16
accepted uniform-400 full-column controls map to their new indices.

All shared aerosol, SIF, and surface values in the new table are cloned from
the corresponding accepted full-column `true_states.dat` row. This detail is
intentional: the older table is the record actually parsed by the production
solver, while its human-readable component catalog carries a few additional
nominal digits. Reconstructing inputs from those extra digits changes two
aerosol AODs by one or two `Float32` ULPs and would make copied O2 components
slightly inconsistent with newly solved CO2 components.

## Why only the CO2 bands are recomputed

The O2 A-band configuration contains no CO2 absorber. The accepted 400 ppm
full-column spectrum with the same surface, aerosol, and SIF state therefore
supplies bit-identical Rayleigh, Cabannes, and RRS arrays for all five new CO2
profiles. Conversely, SIF is absent from the weak and strong CO2 bands. The
producer consequently does the following:

1. For no-SIF states, copy the matching accepted 400 ppm scene, retain all
   three O2 arrays, and recompute weak and strong CO2 only for bottom-layer
   values 360, 380, 420, and 440 ppm.
2. Copy the 400 ppm control without a new RT calculation.
3. Only after every no-SIF state is complete, assemble the SIF-on states from
   accepted SIF-on O2 arrays and the corresponding new no-SIF CO2 arrays.

This leaves 16 non-control atmospheric/surface calculations per machine, or
32 band solves each. It avoids repeating any aerosol/Raman A-band solve.

## Machine partition and state ranges

| Machine | Partition | no-SIF states | SIF-on states |
|---|---|---|---|
| Curry | no aerosol | 001--005, 021--025, 041--045, 061--065 | 006--010, 026--030, 046--050, 066--070 |
| Wurst | aerosol | 011--015, 031--035, 051--055, 071--075 | 016--020, 036--040, 056--060, 076--080 |

Clear scenes are written directly under `truth/`; aerosol scenes are written
under `truth/aerosol_chunked/`, matching the established truth-map layout.
The partitions and indices are disjoint.

The downstream retrieval code now accepts this 80-state table and isolated
campaign paths without changing the full-column defaults. Synthetic OCO
radiances and frozen diagonal noise covariances are generated within
`truth/OCO_radiances/`; retrievals are written under `retrievals/`. No existing
retrieval product is copied into this campaign.

Campaign-aware plotting entry points and their usage are documented in
[`visualization/README.md`](visualization/README.md).  They display the
retrieved layer-16 VMR separately from the resulting column XCO2 and rebuild
bottom-layer truth profiles from the exact 16-layer state rather than treating
the column mean as a vertically uniform VMR.

## Launch sequence

Do not start either truth producer while a retrieval worker is active on the
same host. The launcher waits by default and enforces Curry/clear and
Wurst/aerosol placement. Run the first two commands independently on their
assigned hosts:

```bash
# Curry, after its retrieval worker exits
./RRS_XCO2/scripts/run_bottom_layer_truth_partition.sh none nosif

# Wurst, after its retrieval worker exits
./RRS_XCO2/scripts/run_bottom_layer_truth_partition.sh aerosol nosif
```

The launcher uses physical GPU 1 by default. Select another physical device
with `BOTTOM_TRUTH_PHYSICAL_GPU`; it is hidden behind
`CUDA_VISIBLE_DEVICES`, so Julia still addresses the selected device as
logical device 0. For example, Curry device 0 can run alongside a retrieval
on device 1 with:

```bash
BOTTOM_TRUTH_PHYSICAL_GPU=0 \
BOTTOM_TRUTH_ALLOW_RETRIEVAL_OVERLAP=1 \
./RRS_XCO2/scripts/run_bottom_layer_truth_partition.sh none nosif
```

The overlap override bypasses only the host-wide retrieval check. The
launcher still identifies the selected GPU by UUID and waits until that GPU
has no non-infrastructure compute process; it never assumes that free memory
alone makes a shared device safe. Curry's explicitly named, permanent
`postgres: GPU0 memory keeper` context is ignored, but every other CUDA
process remains blocking.

Wait until **both** no-SIF partitions have passed their terminal validation.
Then assemble the SIF-on states:

```bash
# Curry
./RRS_XCO2/scripts/run_bottom_layer_truth_partition.sh none sif

# Wurst
./RRS_XCO2/scripts/run_bottom_layer_truth_partition.sh aerosol sif
```

For a non-waiting check, set `WAIT_FOR_RETRIEVALS=0`. To validate the launch
without starting Julia, add `BOTTOM_TRUTH_LAUNCH_DRY_RUN=1`. The Julia
producer itself supports `BOTTOM_TRUTH_PLAN_ONLY=1`, which verifies every
source scene and prints all planned actions without writing data.

If the accepted full-column truth is archived before this campaign runs, set
`FULL_COLUMN_TRUTH_ROOT` to its new `truth/` directory. The producer validates
the source state, physics metadata, dimensions, and ABSCO version and will
fail rather than silently use an incompatible archive.

## Publication and failure safety

The new producer never edits a final NetCDF file. It claims each state with an
atomic directory under `truth/.claims/`, builds a private copy under
`.staging/`, validates its metadata and every finite radiance array, and then
publishes it by a same-filesystem rename. Existing complete products are
validated and skipped; existing invalid products are never overwritten.

If a process is interrupted, its claim, `owner.txt`, and staging directory
are intentionally retained. Do not delete a claim until the recorded host and
PID have been checked and the partial file has been inspected. There is no
automatic stale-lock stealing. This prevents Curry and Wurst, or two manual
launches, from silently overwriting one another.

Each completed scene records both `bottom_co2_ppm` and the nominal
dry-column-averaged `xco2_ppm`, the complete 16-layer CO2 profile, source-file
and state-table hashes, and component-reuse provenance. The nominal XCO2 is
computed from the intended ppm values and the production dry-air weights;
the `Float32` realization differs by less than `1e-5` ppm. The 400 ppm control
arrays and every reused O2/SIF-paired CO2 array are checked bit-for-bit before
publication.

The physical motivation, retrieval prior, and future retrieval validation
gates are documented in
[`../inversion/retrieval_setup/BOTTOM_LAYER_XCO2_RETRIEVAL_PLAN.md`](../inversion/retrieval_setup/BOTTOM_LAYER_XCO2_RETRIEVAL_PLAN.md).

## Clear/no-SIF retrieval campaign

The completed clear/no-SIF truth subset contains states
`001--005,021--025,041--045,061--065`. Its instrument processing produces 20
OCO-grid measurement files and 20 frozen-noise files, each with 2,742 samples.
Corrected and uncorrected retrievals use the same definitions as the
full-column campaign, use the `OCO_RRS_synth` analytical Jacobian with nine
streams and `noRS`, and run perturbation 11 (noiseless) before perturbations
1--10.

The campaign-local prior is in `retrieval_setup/apriori_states.nc`. It retains
the mapped ACOS CO2 covariance and the 400 ppm mean, but uses the agreed
surface slope and curvature uncertainties `sigma(P1)=sigma(P2)=2e-3` equally
in the O2 A, weak-CO2, and strong-CO2 bands. The shared full-column prior
remains unchanged. Fixed CO2 layers 1--4 use the tabulated 400 ppm background,
never the column-averaged `xco2_ppm` value.

For the no-SIF-first ensemble, the absolute wavelength-space SIF-slope prior
is centered at zero. Its 1-sigma width is three times the earlier value,
`2.625e-3 mW m-2 sr-1 nm-2`. This admits the zero-SIF truth without retaining
the former four-sigma pull toward a nonzero slope and deliberately loosens
the slope dimension for the first bottom-layer experiment.

The two Curry A100s have disjoint ownership:

| Partition | Physical GPU | Clear/no-SIF states | Products |
|---|---:|---|---:|
| `curry0` | 0 | 001--005, 021--025 | 220 |
| `curry1` | 1 | 041--045, 061--065 | 220 |

Each product count is ten states times eleven perturbations times two
measurement classes. Run or resume a partition with:

```bash
./RRS_XCO2/inversion/run_bottom_layer_retrieval_partition.sh curry0
./RRS_XCO2/inversion/run_bottom_layer_retrieval_partition.sh curry1
```

The second command waits if its assigned GPU is occupied. Atomic
partition-claim directories prevent duplicate launchers, complete files are
validated and skipped on restart, and each corrected/uncorrected pair remains
on one GPU. Device 0 first performs a paired noiseless desert state-043 smoke
solve. Production starts only after `validate_bottom_layer_retrievals.jl`
confirms scientific closure for the corrected result: `retrieval_complete=1`,
`converged=1`, `fit_quality_ok=1`, and OE `outcome=1`, with the stored
per-band chi-squared values consistent with that status. The uncorrected
state-043 result is deliberately allowed to retain scientifically meaningful
RRS model mismatch, but it must be complete, finite, schema-valid, and have a
legal, internally consistent OE outcome. Device 1 waits for this gate and
then skips the two state-043 products already computed on device 0, avoiding
duplicate work.

## Aerosol retrieval campaign, with and without SIF

The Wurst launcher covers all 40 aerosol states and writes 880 products
(40 states times 11 perturbations times two measurement classes). Before it
claims a GPU, it requires the complete 80-state OCO-radiance and frozen-noise
datasets to pass their production validators. A retrieval-specific preflight
then checks the exact aerosol state set, campaign-local paths, 30-parameter
positive-definite prior, paired corrected/uncorrected order, and perturbation
order `11,1,...,10`. It writes the authoritative aerosol-only manifest to
`retrievals/retrieval_manifest_wurst0_aerosol_all.dat`; neither production
worker overwrites the clear-scene manifest.

The static ownership is deliberately slightly asymmetric because physical
GPU 1 remains occupied when GPU 0 starts:

| Partition | Physical GPU | Aerosol/no-SIF states | Aerosol/SIF-on states |
|---|---:|---|---|
| `wurst0` | 0 | 011--015, 031--035, 051 | 016--020, 036--040 |
| `wurst1` | 1 | 052--055, 071--075 | 056--060, 076--080 |

State 013, perturbation 11 is the initial aerosol/no-SIF smoke retrieval.
The corrected and uncorrected members must both be structurally valid, and
the corrected member must converge with acceptable fit quality. Both workers
then compute only their statically assigned no-SIF states. A shared validator
barrier over all 20 no-SIF states prevents either device from entering the
SIF phase early. Only after that barrier passes does GPU 0 run the analogous
state-018 SIF smoke retrieval; GPU 1 waits for that gate before starting its
own SIF states. This ordering holds even when GPU 1 is launched hours later.

Run the workers on Wurst with:

```bash
./RRS_XCO2/inversion/run_bottom_layer_aerosol_retrievals.sh wurst0
./RRS_XCO2/inversion/run_bottom_layer_aerosol_retrievals.sh wurst1
```

The `wurst1` worker waits until physical GPU 1 becomes free. For a read-only
release check that performs every input validation without claiming a GPU or
starting a retrieval, set `BOTTOM_RETRIEVAL_PREFLIGHT_ONLY=1`. Partition
claims under `retrievals/.partition_claims/` record the exact owned states,
physical GPU UUID, host, PID, and phase barriers, preventing two launchers
from owning the same static partition.
