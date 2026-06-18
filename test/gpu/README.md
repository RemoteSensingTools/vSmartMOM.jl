# GPU / Metal tests

These tests are **not** part of the CI unit suite (`test/runtests.jl`). They are
GPU-related and either need GPU/Metal hardware or are GPU-kernel tests kept here
to keep the main suite hardware-free.

## Run

```bash
julia --project=test test/gpu/runtests.jl
```

## Contents

| File | Needs | Notes |
|------|-------|-------|
| `test_mie_gpu.jl` | nothing (KA **CPU** backend) | DoubleSingle/Neumaier precision layer + NAI2 Mie pipeline. Runs anywhere. |
| `test_forward_raman_gpu.jl` | a functional **CUDA** device | RRS forward model on GPU; skips cleanly without CUDA. |
| `test_jacobians_GPU.jl` | a functional **CUDA** device | Linearized RT on GPU; self-skips for missing CUDA / CUBLAS / scalar-indexing paths. |

The runner `cd`s to the `test/` root so the files' relative
`test_parameters/...` paths resolve (those configs stay in the main suite).
