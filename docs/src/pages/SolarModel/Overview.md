# SolarModel

**For:** users supplying solar spectra or checking the package default solar transmission table.

**Next:** [Configure a Scene](../IO/Overview.md), [Core RT Theory](../vSmartMOM/CoreRTTheory.md), [API Reference](../api_reference.md).

The `SolarModel` module provides Planck-spectrum helpers and solar transmission readers used by vSmartMOM examples and RT configurations. It is intentionally small: it loads or generates spectral inputs, while the CoreRT solver applies those inputs through the solar source-function path.

## Spectrum Helpers

| Function | Grid | Output |
| --- | --- | --- |
| `planck_spectrum_wn` | wavenumber in `cm^-1` | black-body radiance in `mW m^-2 sr^-1 cm` |
| `planck_spectrum_wl` | wavelength in micrometers | black-body radiance in `W m^-2 sr^-1 um^-1` |
| `default_solar_spectrum_at_earth` | wavenumber in `cm^-1` | approximate photon spectrum multiplied by the default solar transmission |

The default solar spectrum helper is available as `vSmartMOM.SolarModel.default_solar_spectrum_at_earth`.

## Solar Transmission Tables

Use `solar_transmission_from_file` when you have an explicit two-column table:

```julia
using vSmartMOM

nu_grid = 12900.0:0.02:13100.0
solar = vSmartMOM.SolarModel.solar_transmission_from_file("solar.out", nu_grid)
```

Use `default_solar_transmission` when you want the package default table:

```julia
solar = vSmartMOM.SolarModel.default_solar_transmission(12900.0:0.02:13100.0)
```

## Default Data Resolution

`default_solar_transmission_path` resolves the default table in this order:

1. `ENV["VSMARTMOM_SOLAR_FILE"]`, when set;
2. the registered Julia `solar` artifact, when available;
3. a checksum-verified scratch-space cache downloaded from the legacy vSmartMOM host.

The fallback cache is relocatable and does not write into the package source tree.

## Useful APIs

- `planck_spectrum_wn`, `planck_spectrum_wl`
- `solar_transmission_from_file`
- `default_solar_transmission_path`, `default_solar_transmission`
- `vSmartMOM.SolarModel.default_solar_spectrum_at_earth`
