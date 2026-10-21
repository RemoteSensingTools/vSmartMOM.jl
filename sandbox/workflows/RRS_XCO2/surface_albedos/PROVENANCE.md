# ECOSTRESS-derived truth-map surface albedos

## Purpose

These files define the four Lambertian surface cases used by the 64-scene
OCO-2 truth map. They are derived from spectra downloaded from the NASA/JPL
ECOSTRESS Spectral Library into `~/data/ECOSTRESSlib` and copied into
`RRS_XCO2/data/ecostress/` so the calculation is reproducible on hosts where
the shared-data mount is unavailable.

Run the complete calculation with:

```bash
julia --project=. RRS_XCO2/scripts/fit_ecostress_surfaces.jl
python3 RRS_XCO2/scripts/plot_ecostress_surfaces.py
```

## Input interpretation

Every selected ECOSTRESS text file declares wavelength in micrometres and
reflectance in percent. The script converts these to nanometres and unitless
reflectance. Invalid/header lines are ignored, samples are sorted by increasing
wavelength, and linear interpolation is performed in wavelength space.

Generic classes are represented by unweighted category means:

- **Trees:** mean of all 15 downloaded live-tree VSWIR spectra. This reduces
  dependence on one species and is used both for the tree fraction in mixed
  scenes and for the forest scene.
- **Roofing:** mean of the downloaded metal, asphalt-shingle, and tile spectra.
- **Grassland:** the downloaded *Avena fatua* grass spectrum.
- **Concrete paving:** the downloaded asphaltic paving-concrete spectrum.
- **Sahara dust:** no downloaded entry is explicitly labeled Sahara. The
  implemented proxy is the mean of the two downloaded dry regional aridisols:
  a Syrian light-yellowish-brown calciorthid and a Jordanian pale-brown
  gypsiorthid. This limitation is explicit rather than silently relabeling one
  sample as Saharan dust. Replace these inputs if a true Sahara spectrum is
  supplied later.

Exact filenames are recorded in `ecostress_manifest.txt` and preserved as
inputs under `RRS_XCO2/data/ecostress/`.

## Mixture recipes

Mixing is linear in unitless Lambertian reflectance at every wavelength:

```text
urban  = 0.20 trees + 0.50 roofing + 0.30 concrete_paving
rural  = 0.50 grass + 0.30 roofing + 0.20 trees
desert = 1.00 sahara_dust_proxy
forest = 1.00 trees
```

The category mean is calculated before applying the scene area fraction. All
weights sum to one. No atmospheric correction, continuum removal, smoothing,
or extrapolation is applied.

## NetCDF products

`urban_surface_albedo.nc`, `rural_surface_albedo.nc`,
`desert_surface_albedo.nc`, and `forest_surface_albedo.nc` contain a common
420–2500 nm grid at 1 nm spacing. Each file stores:

- `wavelength(wavelength)` in nm;
- `surface_albedo(wavelength)`, the final unitless mixture;
- each category-mean component used by that scene, with its mixture weight as
  a variable attribute;
- global source, recipe, processing, and generation-history attributes.

The 1 nm grid reflects the information content of the ECOSTRESS inputs. The
surface is evaluated on the RT grid through its smooth Legendre representation;
the NetCDF data should not be interpreted as 0.1 cm⁻¹ laboratory sampling.

## vSmartMOM P0–P2 representation

For each surface and each OCO band, the ECOSTRESS mixture is evaluated on the
same ascending-wavenumber 0.1 cm⁻¹ basis grid used by the truth simulations.
The script solves an ordinary least-squares fit

```text
rho(x) = P0 + P1*x + P2*(3*x^2 - 1)/2,
```

where `x=-1` is the longest-wavelength band edge and `x=+1` the shortest.
`lambertian_legendre_inputs.dat` is the compact human-readable table used to
construct `LambertianSurfaceLegendre([P0,P1,P2])`. It also records RMSE and
maximum absolute fit error. `surface_spectra_and_fits.csv` retains every fine
band-grid target and fitted value for independent verification, and
`surface_legendre_fits.png` visualizes all twelve fits.

