# Run on GPU

**For:** users who want to run the solver on GPU backends.

**Next:** [Quick Start](quickstart.md), [Configure a Scene](IO/Overview.md), [Tutorial: GPU](tutorials/Tutorial_GPU.md).

GPU support is provided through optional backend extensions. CUDA is the mature
NVIDIA path; Metal is a first-pass Apple Silicon path for Float32 runs. The full
task page will cover:

- when `vSmartMOM.Architectures.GPU()` is available;
- when `vSmartMOM.Architectures.MetalGPU()` is available;
- how `array_type(model)` and architecture dispatch select CPU or GPU arrays;
- which workflows are GPU-safe today;
- memory and precision caveats for Float32 and Float64 runs;
- how to fall back cleanly to CPU.

For now, see the long-form [GPU tutorial](tutorials/Tutorial_GPU.md).

## Performance Notes

The main architecture switch is the `radiative_transfer.architecture` field in
the scene configuration, or `params.architecture = GPU()` after loading a
configuration. GPU runs are most useful when enough spectral points, layers, or
viewing geometries are batched to amortize kernel-launch and transfer costs.
Use `CPU()` for small debugging scenes and for workflows that are not yet
GPU-safe.

For Apple Silicon experiments, load Metal.jl and use `MetalGPU()` with
`Float32` scene parameters. The current Metal path focuses on the core batched
matrix multiply and inverse operations; fused kernels can be added after Mac
validation.
