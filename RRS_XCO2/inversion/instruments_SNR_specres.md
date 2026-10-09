# O2 A-band instrument SNR and spectral-resolution comparison

Last reviewed: 2026-08-27

## Purpose

This note compares the signal-to-noise ratio (SNR) and spectral resolution of
current and planned satellite instruments that observe the molecular oxygen
A-band near 0.76 um. The comparison was assembled for the RRS/XCO2 retrieval
study, where the relevant question is whether rotational Raman scattering
(RRS) or other spectral-filling effects would be detectable above instrumental
noise.

FLEX/FLORIS is included even though FLEX is a vegetation-fluorescence mission,
not a greenhouse-gas column mission. Its O2 A-band measurements are directly
relevant because both SIF and RRS fill solar and atmospheric absorption
features.

## Comparability warning

SNR is not a single intrinsic number for a spectrometer. It depends on:

- scene radiance, surface reflectance, and illumination/viewing geometry;
- wavelength and absorption-line depth;
- spectral resolution and sampling;
- spatial footprint and integration time;
- detector/read/background noise and any spatial or spectral averaging.

Consequently, a larger native-sample SNR does not necessarily imply more
information. A wide spectral response collects more photons but smooths the
line structure that constrains atmospheric retrievals. Values below are the
published native-sample specifications and should not be compared without
their reference radiances and spectral widths.

## Instrument comparison

| Instrument | Status at review | O2 A-band resolution | Published A-band SNR | Spatial sampling | Interpretation |
|---|---|---:|---:|---:|---|
| OCO-2 | Operating | about 0.042 nm; resolving power greater than 17,000 | continuum SNR greater than 400; strongly radiance dependent | less than 3 km2 per sounding | High-resolution reference instrument. |
| OCO-3 | Operating | essentially the OCO-2 spare spectrometer | comparable to OCO-2 rather than a higher-SNR successor | ISS-dependent footprint | Not an independent sensitivity improvement over OCO-2. |
| TanSat/ACGS | Operating | about 0.044 nm | 360 at 5.8e19 photons s-1 m-2 sr-1 um-1 | about 2 x 3 km | Very close to OCO-2 resolution, with comparable rather than greater SNR. |
| MicroCarb | Operating; calibration/validation products still maturing | resolving power about 26,040, corresponding to about 0.029 nm at 763.5 nm | 61--106 at minimum, 194--320 at mean, and 494--753 at maximum reference radiance, depending on channel | about 4.5 x 9 km nominal | Bright-scene native SNR can exceed OCO-2 while retaining finer spectral resolution, but the larger footprint collects substantially more light. A released on-orbit A-band SNR validation is needed for a realized comparison. |
| CO2M/CO2I | Planned; first launch listed for late 2027 | 0.11 nm | 330 at 6.4e19 photons s-1 m-2 sr-1 um-1 | about 2.07 x 1.97 km | Comparable native SNR at the stated radiance, but about 2.6 times coarser than OCO-2. |
| GOSAT-GW/TANSO-3 | Operating | requirement less than 0.05 nm near 0.76 um | no sufficiently documented, validated A-band SNR comparison located | 10 km wide mode; selectable 1--3 km focus mode | Official archive states that its ATBD and validation results were not yet available as of July 2026, so a claim of superiority is premature. |
| TanSat-2 O2-A payload | Preliminary design | 0.12 nm FWHM; 0.04 nm sampling | 500 at 6.4e19 photons s-1 m-2 sr-1 um-1 | preliminary 2 km specification | Raw SNR exceeds OCO-2 at the common reference radiance, but the spectral response is almost three times wider. |
| FLEX/FLORIS | Prelaunch; ESA lists launch in 2026 | 0.28 nm FWHM and about 0.093--0.10 nm sampling in the central 759--769 nm region; about 0.474 nm in adjacent A-band intervals | requirements span roughly 115--455 through the central A-band and reach 1015 in brighter A-band continuum/shoulder intervals | 300 x 300 m | Headline continuum SNR is much higher than OCO-2, but central resolution is about 6.7 times coarser. Particularly important for RRS/SIF filling-in studies. |

## Common-radiance OCO-2 check

The CO2M and preliminary TanSat-2 specifications both use a reference radiance
of 6.4e19 photons s-1 m-2 sr-1 um-1. Evaluating OCO-2 Level-1B ATBD Eq. (3-8)
at the same radiance with the wavelength-dependent representative coefficients
stored in
`RRS_XCO2/inversion/instrument/representative_snr_coefficients.nc` gives an
OCO-2 A-band SNR range of 298--366, with a median of 338 across the synthetic
sampling grid.

This makes the native-sample interpretation more precise:

- CO2M's SNR of 330 is essentially comparable to representative OCO-2 at the
  same radiance, rather than clearly higher.
- TanSat-2's preliminary SNR of 500 is higher per native sample, but its
  0.12 nm response collects photons over a much wider spectral interval.

The local calculation uses the pointwise median of 32 OCO-2 coefficient
spectra (four nadir L1B files times eight footprints), so it is representative
rather than universal.

## Approximate equal-resolution comparison

For a rough photon-noise-limited comparison only, SNR scales as the square
root of spectral-bin width:

```text
SNR_at_OCO_width ~= SNR_native * sqrt(0.042 nm / FWHM_native).
```

Applying that heuristic gives:

| Instrument case | Native SNR and FWHM | Approximate SNR at 0.042 nm |
|---|---:|---:|
| CO2M | 330 at 0.11 nm | 204 |
| TanSat-2 preliminary | 500 at 0.12 nm | 296 |
| FLEX central A-band upper value | 455 at 0.28 nm | 176 |
| FLEX bright outer A-band | 1015 at 0.474 nm | 302 |

These scaled values are not instrument predictions: read noise, background
noise, throughput, sampling, and instrument-line-shape differences violate the
simple square-root model. They illustrate why the native SNR alone should not
be treated as an information-content ranking.

## Consequences for the RRS study

1. **MicroCarb is the strongest plausible higher-SNR, high-resolution GHG
   case.** A sensitivity experiment near SNR 600 is worthwhile, but should
   ultimately use MicroCarb's wavelength-dependent noise and instrument line
   shape rather than constant Gaussian noise.
2. **FLEX is the most directly relevant non-GHG case.** Its continuum SNR can
   reach about 1015, while its central A-band SNR is much lower and varies with
   line depth. Since both SIF and RRS appear as spectral filling, omitted RRS
   can map into retrieved SIF rather than merely producing random-looking
   residuals.
3. **FLEX must be simulated at its native resolution.** Convolve the existing
   high-resolution RRS truth spectra to the FLORIS response (central FWHM about
   0.28 nm), sample at about 0.0933 nm, and apply wavelength-dependent FLORIS
   noise. Simply replacing the OCO-2 SNR with 1015 would be physically
   misleading.
4. **The O2-B band may also matter for FLEX.** The present truth-map work covers
   O2-A only; a complete FLEX bias study would eventually need a corresponding
   O2-B RRS calculation because FLORIS retrieves fluorescence using both oxygen
   bands.

## Primary and mission sources

- OCO-2 on-orbit performance: [Crisp et al. (2017), NASA mission publication page](https://airbornescience.nasa.gov/acepwg/content/The_on-orbit_performance_of_the_Orbiting_Carbon_Observatory-2_OCO-2_instrument_and_its)
- OCO-2 A-band width and radiometric model: [Frankenberg et al. (2014), NASA NTRS](https://ntrs.nasa.gov/citations/20140012653)
- TanSat/ACGS specifications: [National Satellite Meteorological Center instrument page](https://www.nsmc.org.cn/nsmc/en/instrument/ACGS.html)
- MicroCarb detailed performance: [CNES MicroCarb team, IWGGMS-17 instrument presentation](https://cce-datasharing.gsfc.nasa.gov/files/conference_presentations/Talk_Jouglet_35_25.pdf)
- MicroCarb current mission status: [CNES MicroCarb mission page](https://cnes.fr/en/projects/microcarb)
- CO2M NIR performance requirements: [EUMETSAT CO2M performance-requirements document](https://www-cdn.eumetsat.int/files/2024-08/CO2M%20Ground-Based%20Network%20Reference%20Product%20Performance%20Requirements.pdf)
- CO2M schedule/status: [ESA CO2M instrument update, June 2026](https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Carbon_dioxide_monitoring_satellite_s_instrument_passes_vacuum_test)
- GOSAT-GW spectral configuration: [NIES GOSAT-GW/TANSO-3 overview](https://gosat-gw.nies.go.jp/en/gosat-gw02.html)
- GOSAT-GW documentation status: [NIES TANSO-3 technical-document archive](https://product.gosat-gw.nies.go.jp/document/technicaldoc/)
- TanSat-2 preliminary O2 specifications: [Zhao et al. (2025), Atmospheric Measurement Techniques](https://amt.copernicus.org/articles/18/3647/2025/)
- FLEX mission sampling and spatial specifications: [ESA FLEX facts and figures](https://www.esa.int/Applications/Observing_the_Earth/FutureEO/FLEX/Facts_and_figures)
- FLORIS spectral and SNR requirements: [Coppo et al. (2017), Fluorescence Imaging Spectrometer for ESA FLEX Mission](https://www.mdpi.com/2072-4292/9/7/649)
