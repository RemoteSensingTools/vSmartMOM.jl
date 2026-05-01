# Absorption API

`Absorption` provides HITRAN line-list reading, line-shape models, interpolation
models, and absorption cross-section evaluation.

## HITRAN Data Management

```@docs
vSmartMOM.artifact
vSmartMOM.fetch_hitran
vSmartMOM.fetch_hitran_by_ids
vSmartMOM.set_hitran_edition!
vSmartMOM.get_hitran_edition
vSmartMOM.available_hitran_editions
vSmartMOM.hitran_info
vSmartMOM.hitran_is_cached
```

## Cross-Section Workflow

```@docs
vSmartMOM.Absorption.compute_absorption_cross_section
vSmartMOM.Absorption.absorption_cross_section
vSmartMOM.Absorption.read_hitran
vSmartMOM.Absorption.make_hitran_model
vSmartMOM.Absorption.make_interpolation_model
```

## Line-Shape Types

```@docs
vSmartMOM.Absorption.AbstractBroadeningFunction
vSmartMOM.Absorption.Voigt
vSmartMOM.Absorption.Lorentz
vSmartMOM.Absorption.Doppler
```
