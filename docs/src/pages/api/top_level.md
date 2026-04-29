# Top-Level API

These are the functions and types most users call directly after `using vSmartMOM`.
They load a scene, construct an RT model, run the solver, and select the compute
architecture.

## Scene Loading

```@docs
vSmartMOM.read_parameters
vSmartMOM.parameters_from_file
vSmartMOM.parameters_from_source
vSmartMOM.parameters_from_dict
vSmartMOM.parameters_from_yaml
vSmartMOM.default_parameters
```

## Model Construction and Solver Entry Points

```@docs
vSmartMOM.model_from_parameters
vSmartMOM.rt_run
```

## Architectures

```@docs
vSmartMOM.CPU
vSmartMOM.GPU
vSmartMOM.MetalGPU
vSmartMOM.default_architecture
vSmartMOM.array_type
```
