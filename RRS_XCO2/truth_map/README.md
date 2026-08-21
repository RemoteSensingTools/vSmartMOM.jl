# 64-scene high-resolution truth map

`true_states.dat` is the authoritative whitespace-separated state table. Its
nesting order is surface, aerosol, SIF, then XCO2 (XCO2 changes fastest).
For a compact description of the repeated surface and aerosol definitions,
see `scene_components.dat`; this avoids decoding those definitions from all
64 rows of the state table.

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

### Chunked aerosol/Raman production

Aerosol-on O2 A-band RRS scenes cannot be solved over the complete
Raman-shouldered spectrum in one allocation. Run
`scripts/generate_truth_map_aerosol_chunked.jl` for those 32 states. It splits
the O2 output band into core intervals, adds the full Raman shoulder to both
sides of every interval, solves the expanded grid, and writes only the core
into the final NetCDF variable. Weak and strong CO2 spectra are also chunked,
but do not receive Raman shoulders.

The runner reuses results across physically identical states: O2 is keyed by
surface and SIF (CO2 has no A-band absorption), while each CO2 band is keyed by
surface and XCO2 (SIF is absent). Results are written incrementally and a JLD2
checkpoint is updated atomically after every physical-state/chunk unit.

```bash
CUDA_DEVICE=1 \
O2_CHUNK_POINTS=256 CO2_CHUNK_POINTS=512 \
julia --project=. RRS_XCO2/scripts/generate_truth_map_aerosol_chunked.jl
```

The default output is `truth_map/aerosol_chunked/`. Set `TRUTH_OUT` to choose
another destination. `FORCE=1` recreates all selected aerosol scene files and
moves the previous checkpoint aside.

The original 12-layer dataset made before correction of the altitude-profile
coordinate bug is retained under `sims_12layer_uncorrected/` for diagnostic
provenance only. The corrected 12-layer dataset is generated after the
corrected 16-layer run and stored under `sims_12layer/`; it must not be mixed
with the current 16-layer files in this directory.

`sim_wavelength.nc` stores the reported 0.1 cm-1 grids for O2 A, weak CO2, and
strong CO2. O2 A is solved with an additional 250 cm-1 shoulder on both sides;
shoulders are discarded from `hiressim_NNN.nc`. Each scene file stores the TOA
I,Q,U Rayleigh/noRS, Cabannes/RRS-elastic, and rotational-Raman radiances for
O2 A, and only the Rayleigh/noRS radiance for both CO2 bands. No redundant
Cabannes+RRS sum is stored.

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
