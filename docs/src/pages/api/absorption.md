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
`radiative_transfer.absorption.CEF`.

```@docs
vSmartMOM.Absorption.HumlicekErrorFunction
vSmartMOM.Absorption.HumlicekWeidemann32VoigtErrorFunction
vSmartMOM.Absorption.HumlicekWeidemann32SDErrorFunction
vSmartMOM.Absorption.CPF12ErrorFunction
vSmartMOM.Absorption.ErfcHumliErrorFunctionVoigt
vSmartMOM.Absorption.ErfcHumliErrorFunctionSD
vSmartMOM.Absorption.ErfcErrorFunction
```
