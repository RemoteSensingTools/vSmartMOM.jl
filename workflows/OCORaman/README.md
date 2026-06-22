# OCORaman Workflow

Research workflow for regenerating O2 A-band Raman/Cabannes outputs, sanity
checks, old OCO-style plots, and later XCO2 bias experiments against the
current `vSmartMOM` source tree.

This directory is intentionally separate from `src/` and `test/`: generated
LuT files, figures, and retrieval experiments can move fast here without
turning research products into package internals.

## Setup

From the repository root:

```bash
julia --project=workflows/OCORaman -e 'using Pkg; Pkg.instantiate()'
```

The workflow project uses the local checkout through:

```toml
[sources]
vSmartMOM = {path = "../.."}
```

Optional packages from the old `OCORaman/src/OCOPlots` layer, such as
`InstrumentOperator`, `NCDatasets`, `Interpolations`, `HDF5`, and retrieval
plotting helpers, should be added only when the corresponding recreated script
needs them. The base environment is kept small so Cabannes sanity checks and
GPU output generation are not blocked by old plotting dependencies.

## 1. Cabannes/Rayleigh Sanity Numbers

```bash
julia --project=workflows/OCORaman \
  workflows/OCORaman/scripts/01_cabannes_sanity.jl
```

Useful overrides:

```bash
WAVELENGTH_NM_LIST=760,765 TEMPERATURE_K=296 \
  julia --project=workflows/OCORaman \
  workflows/OCORaman/scripts/01_cabannes_sanity.jl
```

Outputs go under `workflows/OCORaman/output/` by default, or under
`OCORAMAN_OUTPUT_ROOT` when set.

## 2. One O2 A-Band GPU Smoke Output

This wraps the current unified driver:

`test/benchmarks/creategrid_O2Aband_RamanSIF.jl`

```bash
GRID_YAML=test/test_parameters/O2_parameters2_SIF_grid.yaml \
SZA=19 ALBEDO=0.0 SIF_ON=1 \
julia --project=workflows/OCORaman \
  workflows/OCORaman/scripts/02_o2a_smoke_output.jl
```

The default `O2_parameters2_SIF_grid.yaml` is configured for
`Architectures.GPU()` and `Float64`.

## 3. Plan Or Run A Small LuT Grid

The grid runner defaults to `DRY_RUN=1` so it prints cases without launching
large jobs.

```bash
OCORAMAN_SZA_LIST=0,30,50,70 \
OCORAMAN_ALBEDO_LIST=0.0,0.1,0.3,0.5 \
julia --project=workflows/OCORaman \
  workflows/OCORaman/scripts/03_run_o2a_grid.jl
```

Launch the cases explicitly:

```bash
DRY_RUN=0 \
OCORAMAN_SZA_LIST=0,30,50,70 \
OCORAMAN_ALBEDO_LIST=0.0,0.1,0.3,0.5 \
julia --project=workflows/OCORaman \
  workflows/OCORaman/scripts/03_run_o2a_grid.jl
```

## 4. Quicklook Plot

```bash
GRID_FILE=/path/to/rrs_grid_sza19p0_alb0p0_sifon.jld2 \
julia --project=workflows/OCORaman \
  workflows/OCORaman/scripts/04_plot_grid_quicklook.jl
```

or:

```bash
GRID_OUTDIR=/path/to/grid/output \
julia --project=workflows/OCORaman \
  workflows/OCORaman/scripts/04_plot_grid_quicklook.jl
```

## 5. Full O2 A-Band Float32 NetCDF LuT

The current replacement for hand-edited `O2_parameters2_SIF_grid.yaml` runs is:

```bash
CUDA_DEVICE=1 DRY_RUN=1 julia --project=. \
  workflows/OCORaman/scripts/createRamanLUT_O2A.jl
```

It generates two broad O2 A-band chunks with 250 cm^-1 Raman shoulders, writes
only the unshouldered 740-785 nm output, uses `Float32`, and stores
Rayleigh/Cabannes/RRS Stokes `{I,Q,U}`, `F0`, zero `SIF0`, `tau_rayl`,
`tau_abs`, and the generated `p/T/q` profiles in NetCDF. The dedicated input
template is `workflows/OCORaman/config/paramsRamanLUT_O2A.yaml`, whose `p/T/q`
values were copied from `test/test_parameters/ParamsEMIT_MODTRANcomp_newLUT.yaml`
and trace to the AFGL 1986 Midlatitude Winter profile used by the EMIT/MODTRAN
comparison workflow. O2 absorption uses ABSCO; q-driven H2O absorption uses the
rebuilt HITRAN H2O LUT.

Suggested distributed pressure slices:

```bash
CUDA_DEVICE=1 PSURFS=1000 OUT_NC=/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_psurf1000.nc \
  julia --project=. workflows/OCORaman/scripts/createRamanLUT_O2A.jl

CUDA_DEVICE=1 PSURFS=750 OUT_NC=/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_psurf750.nc \
  julia --project=. workflows/OCORaman/scripts/createRamanLUT_O2A.jl

CUDA_DEVICE=1 PSURFS=500 OUT_NC=/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_psurf0500.nc \
  julia --project=. workflows/OCORaman/scripts/createRamanLUT_O2A.jl
```

Run `DRY_RUN=1` first: the full default grid is very large
(`21 albedos x 14 SZA x 1 SIF state x 3 psurf`) and the script prints a
runtime estimate based on the latest full-band timing baseline.

## Workflow Notes

- Keep generated `.jld2`, NetCDF, PNG, and PDF files out of git by default.
- Each script writes `run_manifest.json` with the git SHA and dirty-tree status.
- The current unified `creategrid_O2Aband_RamanSIF.jl` driver runs `iBand = 1`.
  Recreating the old full LuT generator means restoring a controlled loop over
  bands, surface/albedo, aerosol, SZA, and output naming on top of this scaffold.
- For paper figures, promote final figures deliberately by moving/copying them
  into a tracked docs or manuscript directory later.
- For OCO XCO2 bias work, use RRS output as truth first and noRS/elastic
  linearized retrieval as the retrieval model:
  `Δx ≈ (K' S_e^-1 K + S_a^-1)^-1 K' S_e^-1 Δy`.
