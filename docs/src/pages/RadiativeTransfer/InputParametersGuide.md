# RT Parameters Guide

## **radiative_transfer** group

**spec_bands**: List of spectral bands to be included in the computation, written in the form, `ν_start:ν_step:ν_end`. Units in cm⁻¹. 

**surface**: Type of reflectance surface to model at bottom of atmosphere. Currently only supporting Lambertian surface types. Should be written in the form, `LambertianSurfaceScalar(albedo)`. 

**quadrature_type**: Quadrature type for RT streams. Should be one among [`RadauQuad`, `GaussQuadHemisphere`, `GaussQuadFullSphere`].

**polarization_type**: Level of polarization to simulate. Should be one among [`I`, `IQ`, `IQU`, `IQUV`].

**max_m**: Hard cutoff for maximum number of Fourier moments to loop over. Should be a positive integer. 

**Δ_angle**: Exclusion angle for forward peak in phase truncation. Units in `degrees`. 

**l_trunc**: Truncation length for legendre terms. Should be a positive integer. 

**depol**: Depolarization factor to use in Rayleigh calculations. 

**float_type**: Float type to use in the RT simulations. Should be one of [`Float64`, `Float32`]. 

**architecture**: Hardware architecture to use for calculations. Should be one among [`CPU`, `GPU`].

## **geometry** group

**sza**: Solar zenith angle, units in `degrees`. 

**vza**: List of viewing zenith angles, units in `degrees`. 

**vaz**: List of viewing azimuthal angles, units in `degrees`. Should be same length as **vza**. 

**obs_alt**: Altitude of observer, in `Pascal`. 

## **atmospheric_profile** group

**T**: Temperature at each layer *midpoint*, as list from TOA to BOA, in `Kelvin`. 

**p**: Pressure at each layer *boundary*, as list from TOA to BOA, in `hPa`. Should be one element greater in length than T. 

**q** (optional): Specific humidity at each layer *midpoint*. Units in `g/Kg`. 

**profile_reduction**: Length of profile reduction to be performed internally. Must be positive integer less than length of T, *or* -1 for no reduction. 

## **absorption** group (optional)

**molecules**: Which molecules should be used for absorption calculations in each band. Should be a list of lists of molecule names (outer list over number of bands, inner list to hold molecules). Use `RadiativeTransfer.Absorption.show_molecules()` to see a list of valid molecules. 

**vmr**: A dictionary of volume-mixing ratios that match every individual molecule listed in **molecules** to a vmr. A vmr can be either a single number, or an array to interpolate from TOA to BOA, according to the number of layers specified by the p/T grids. 

**broadening_function**: Type of broadening function to use in absorption calculations. Should be one of [`Voigt()`, `Lorentz()`, `Doppler()`]. Please note the parentheses. 

**CEF**: Complex error function for Voigt calculations. Default should generally be used, `HumlicekWeidemann32SDErrorFunction()`. Please see Absorption/complex_error_functions.jl for details. 

**wing_cutoff**: Wing cutoff to use in the line-by-line absorption calculations. Should be a positive number in units, `cm⁻¹`. 

## **scattering** group (optional)

**aerosols**: A list of scattering aerosols and their properties. Each aerosol should be a dictionary of key-value pairs for the following keys: 
- **τ_ref**: Reference τ
- **μ**: Log mean radius (µm)
- **σ**: Log stddev of radius (µm)
- **nᵣ**: Real part of refractive index
- **nᵢ**: Imag part of refractive index
- **p₀**: Pressure peak (Pa)
- **σp**: Pressure peak width (Pa)

**r_max**: Maximum aerosol particle radius for quadrature points/weights. Units in `µm`. 

**nquad_radius**: Number of quadrature points for aerosol radius. Should be a positive integer.

**λ_ref**: Reference wavelength. Units in `µm`. 

**decomp_type**: Fourier decomposition method to use. Should be one of [`NAI2()`, `PCW()`]. NAI2 should be used for most cases. Please note the parentheses. 