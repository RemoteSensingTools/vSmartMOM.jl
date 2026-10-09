# Land-glint truth-map extension

The former land-glint products were generated before the common ABSCO v5.2
truth/retrieval closure and are archived under
`RRS_XCO2/obsolete/pre_absco_closure_20260829/truth_map/land_glint/`. The
active geometry directories will be regenerated after the clean nadir map;
until then, this directory intentionally contains documentation only.

This tree contains off-nadir, principal-plane extensions of the 16-layer truth
map. Five native nine-stream Gauss-Legendre geometries are planned:

| Geometry directory | SZA | VZA | Relative azimuth |
|---|---:|---:|---:|
| `sza10p237292_vza10p237292_relaz00/` | 10.237292 deg | 10.237292 deg | 0 deg |
| `sza23p362328_vza23p362328_relaz00/` | 23.362328 deg | 23.362328 deg | 0 deg |
| `sza36p226627_vza36p226627_relaz00/` | 36.226627 deg | 36.226627 deg | 0 deg |
| `sza48p537727_vza48p537727_relaz00/` | 48.537727 deg | 48.537727 deg | 0 deg |
| `sza60_vza60_relaz00/` | 60 deg | 60 deg | 0 deg |

Here VZA is the outgoing-light direction with respect to zenith, SZA is the
incoming-light direction with respect to nadir, and land glint therefore uses
relative azimuth 0 degrees in the convention specified for this experiment.
The surface remains Lambertian. Only SIF-off states are included initially.

Each geometry has separate `aerosol_chunked/` and `no_aerosol/` products. The
original truth-map indices are preserved: aerosol/SIF-off scenes are 009--012,
025--028, 041--044, and 057--060; aerosol-free/SIF-off scenes are 001--004,
017--020, 033--036, and 049--052. The four indices within each surface are the
380, 400, 420, and 440 ppm CO2 cases.

Run one phase with, for example:

```bash
CUDA_VISIBLE_DEVICES=1 CUDA_DEVICE=0 \
  GLINT_STREAM_INDEX=9 GLINT_PHASE=aerosol \
  julia --project=. RRS_XCO2/scripts/generate_truth_map_land_glint.jl

CUDA_VISIBLE_DEVICES=1 CUDA_DEVICE=0 \
  GLINT_STREAM_INDEX=9 GLINT_PHASE=no_aerosol \
  julia --project=. RRS_XCO2/scripts/generate_truth_map_land_glint.jl
```

The `CUDA_VISIBLE_DEVICES=1` restriction above selects physical device 1 on a
two-GPU host; `CUDA_DEVICE=0` then selects the sole visible device. This avoids
even an initialization context on physical device 0.

For a resumable sequence of geometries on one GPU, use
`scripts/run_land_glint_native_sequence.jl`. Its
`GLINT_AEROSOL_STREAM_INDICES` and `GLINT_NO_AEROSOL_STREAM_INDICES`
variables accept comma-separated subsets of indices 5 through 9. The current
production dispatch serializes the four additional aerosol angles on wurst
after the active 60-degree aerosol job exits, while all five aerosol-free
geometries run independently on curry physical device 1.

### Native-angle O2 chunk size

All five viewing directions are selected from the nine Gauss-Legendre streams.
The wrapper derives a Float32 zenith angle whose `cosd` exactly round-trips to
the corresponding Float32 quadrature node. No zero-weight viewing node is
therefore appended: the IQU operator remains 27 x 27 for every geometry. The
aerosol runs are fixed to nine streams and 256 retained O2 points per chunk,
giving 11 O2 spectral chunks.

The 256-point Raman solve peaks near 45.5 GiB. It fits on wurst's approximately
46 GiB physical device 1 but not on curry's 40 GiB physical device 1. Native
11-chunk aerosol production is therefore serialized on wurst rather than
silently changing the requested chunk width on curry. Aerosol-free production
uses the smaller five-stream solve and fits comfortably on curry (about
18.3 GiB in the first 60-degree case).

The aerosol calculation uses the same 0.1 cm^-1 basis grids, 234 cm^-1 O2
Raman shoulders, eight-point strong-CO2 convolution shoulder, Float32
arithmetic, 16 layers, aerosol optical properties, delta-BGE truncation,
stream configuration, and resumable chunking as the accepted nadir production
truth map. Every checkpoint, wavelength file, and scene file records the actual
SZA, VZA, and relative azimuth.
