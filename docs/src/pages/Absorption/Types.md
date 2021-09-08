# Absorption Module Methods & Types

## Downloading HITRAN Data

```@docs
artifact
```

## Reading HITRAN Files

```@docs
read_hitran
```

## Defining Models with Cross-Section Parameters

```@docs
make_hitran_model
make_interpolation_model
```

## Computing Absorption Cross-Sections

```@docs
absorption_cross_section
```

## Types

### Hitran Data Structure Type

```@docs
HitranTable
```

### Broadening Function Types
```@docs
Doppler
Lorentz
Voigt
```

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

