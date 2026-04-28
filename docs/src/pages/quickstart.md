# Quick Start

**For:** first-time users who want one successful forward radiative-transfer run.

**Next:** [Configure a Scene](IO/Overview.md), [Compute Jacobians](jacobians.md), [Core RT Theory](vSmartMOM/CoreRTTheory.md).

This page is the 5-minute CPU path. It will be expanded into the canonical runnable example in the next documentation pass. For now, the long-form tutorial remains available at [Tutorial: Quick Start](tutorials/Tutorial_QuickStart.md).

The intended workflow is:

```julia
using vSmartMOM

params = read_parameters("path/to/scene.yaml")
params.architecture = vSmartMOM.Architectures.CPU()

model = model_from_parameters(params)
R, T = rt_run(model)
```

The next pass will replace the placeholder path with a small shipped CPU scene and add a docs smoke test for the exact code block.
