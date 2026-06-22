# Local-only tests

These are **not** part of the CI / GitHub-Actions suite (`test/runtests.jl`).
They need resources the CI runners don't have — local data, external LUTs, or a
GPU — so they can only be run on a machine that has the prerequisites.

## Run

```bash
julia --project=test test/local/runtests.jl        # everything below
julia --project=test test/local/gpu/runtests.jl    # just the GPU/Metal tests
```

## Contents

### `gpu/` — GPU / Metal tests
| File | Needs |
|------|-------|
| `test_mie_gpu.jl` | nothing (KernelAbstractions **CPU** backend) |
| `test_forward_raman_gpu.jl` | a functional **CUDA** device (guarded) |
| `test_jacobians_GPU.jl` | a functional **CUDA** device (guarded) |

### `test_parameters/` — data/feature-dependent configs
The forward/lin parity runner builds whatever your environment supports and
reports what the rest need.

| Group | Needs |
|-------|-------|
| H2O-override scenes (`2Band`, `3Band`, `S2Band`, `ParamsEMIT_MODTRANcomp*`, …) | a **config-level H2O override** — NOT yet implemented (H2O is driven by `atmospheric_profile.q`; listing H2O is rejected). Parked until the override lands. |
| `*_newLUT*`, `*SIF_grid*`, `O2Parameters*`, `CO2W*` | external **ABSCO/HITRAN LUTs** (`VSMARTMOM_HITRAN_LUT_DIR` or `/…/data/HITRAN_LUTs`) not shipped in the repo. |
| legacy scenes (`ThreeBandsParameters`, `O2ParametersVS`, …) | kept for reference. |

### scratch scripts (`scratch_*.jl`)
Dev exploration scripts — not run by any suite.

## Why this split exists

The CI suite must build every config in `test/test_parameters/` (the parity
test fails loudly if one doesn't). A shared-FS dev box can "see"
`/home/.../data` LUTs that CI's clean runners cannot, so configs needing local
data must NOT live in the CI folder — they live here.

> AtmosphericAbsorption.jl carries its own LBL/absorption unit tests, so the
> LUT-heavy absorption scenes are intentionally not duplicated in the CI suite.
