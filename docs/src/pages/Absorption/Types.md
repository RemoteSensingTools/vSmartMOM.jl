# Absorption Module Methods & Types

## Downloading HITRAN Data

Use `artifact("CO2")`, `artifact("H2O")`, and related artifact names to resolve
packaged HITRAN `.par` files. See [HITRAN Data Management](HITRAN_Data.md) and
the [Library](../api_reference.md) for the canonical artifact helpers.

## Reading HITRAN Files

Use `read_hitran` to select molecule/isotopologue records and spectral windows
from HITRAN `.par` files.

## Defining Models with Cross-Section Parameters

Use `make_hitran_model` for line-by-line cross sections and
`make_interpolation_model` for lookup-table interpolation workflows.

## Computing Absorption Cross-Sections

Use `absorption_cross_section(model, ν, pressure, temperature)` to evaluate the
selected line-shape or interpolation model on a wavenumber grid.

## Types

### Hitran Data Structure Type

```@docs
HitranTable
```

### Broadening Function Types

The primary broadening models are `Doppler`, `Lorentz`, and `Voigt`. The
canonical docstrings are grouped in the [Library](../api_reference.md).

### Complex Error Function Types

```@docs
HumlicekErrorFunction
HumlicekWeidemann32VoigtErrorFunction
HumlicekWeidemann32SDErrorFunction
CPF12ErrorFunction
ErfcHumliErrorFunctionVoigt
ErfcHumliErrorFunctionSD
ErfcErrorFunction
```

### Cross-Section Model Types

```@docs
HitranModel
InterpolationModel
```
