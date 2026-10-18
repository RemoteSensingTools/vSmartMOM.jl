# Absorption API

!!! warning "Internal absorption API is deprecated"
    AtmosphericAbsorption.jl is the supported absorption dependency for vSmartMOM going
    forward. The standalone line-by-line model here (`make_hitran_model`, `HitranModel`) is
    deprecated in 2.2 and will be removed in 3.0; use
    [AtmosphericAbsorption.jl](https://github.com/RemoteSensingTools/AtmosphericAbsorption.jl)
    for line-by-line calculations, ABSCO HDF tables, and GPU execution. YAML `LUTfiles`
    read OCO-2 ABSCO HDF files, AtmosphericAbsorption `.absco` tables, and vSmartMOM JLD2
    tables from `make_interpolation_model`, which remains the table builder in 2.x.

```@docs
vSmartMOM.Absorption
```

`Absorption` provides HITRAN line-list reading, line-shape models, interpolation
models, and absorption cross-section evaluation for legacy callers.

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

## HITRAN Data Structure

```@docs
vSmartMOM.Absorption.HitranTable
```

## Cross-Section Model Types

```@docs
vSmartMOM.Absorption.HitranModel
vSmartMOM.Absorption.InterpolationModel
```

## Complex Error Function Types

The complex error function evaluated inside the Voigt line shape is
selectable via the `CEF` field on `make_hitran_model` /
`absorption.CEF`.

```@docs
vSmartMOM.Absorption.HumlicekErrorFunction
vSmartMOM.Absorption.HumlicekWeidemann32VoigtErrorFunction
vSmartMOM.Absorption.HumlicekWeidemann32SDErrorFunction
vSmartMOM.Absorption.CPF12ErrorFunction
vSmartMOM.Absorption.ErfcHumliErrorFunctionVoigt
vSmartMOM.Absorption.ErfcHumliErrorFunctionSD
vSmartMOM.Absorption.ErfcErrorFunction
```
