# References

For BibTeX entries, use [`CITATION.bib`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/CITATION.bib) at the repository root.

## How To Cite

If you use the package itself, cite the JOSS software paper:

- Jeyaram, R., Sanghavi, S., and Frankenberg, C. (2022). vSmartMOM.jl: An open-source Julia package for atmospheric radiative transfer and remote sensing tools. *Journal of Open Source Software*, 7(80), 4575. https://doi.org/10.21105/joss.04575

Use the method references below when your work depends on a specific solver component.

## Matrix-Operator RT and Jacobians

- Sanghavi, S. V., Martonchik, J. V., Davis, A. B., and Diner, D. J. (2013). Linearization of a scalar matrix operator method radiative transfer model with respect to aerosol and surface properties. *Journal of Quantitative Spectroscopy and Radiative Transfer*, 116, 1-16. https://doi.org/10.1016/j.jqsrt.2012.10.021
- Sanghavi, S., Davis, A. B., and Eldering, A. (2014). vSmartMOM: A vector matrix operator method-based radiative transfer model linearized with respect to aerosol properties. *Journal of Quantitative Spectroscopy and Radiative Transfer*, 133, 412-433. https://doi.org/10.1016/j.jqsrt.2013.09.004
- Sanghavi, S. and Stephens, G. (2015). Adaptation of the delta-m and delta-fit truncation methods to vector radiative transfer: Effect of truncation on radiative transfer accuracy. *Journal of Quantitative Spectroscopy and Radiative Transfer*, 159, 53-68. https://doi.org/10.1016/j.jqsrt.2015.03.007

## Mie Phase-Matrix Fourier Decomposition

- Sanghavi, S. (2014). Revisiting the Fourier expansion of Mie scattering matrices in generalized spherical functions. *Journal of Quantitative Spectroscopy and Radiative Transfer*, 136, 16-27. https://doi.org/10.1016/j.jqsrt.2013.12.015 (NAI-2 vs PCW comparison; corrected Domke formalism.)

## Raman Scattering

- Sanghavi, S. (2022). Raman scattering in the Earth's atmosphere, Part I: Optical properties. *Journal of Quantitative Spectroscopy and Radiative Transfer*, 291, 108328. https://doi.org/10.1016/j.jqsrt.2022.108328
- Sanghavi, S. and Frankenberg, C. (2023). Raman scattering in the Earth's atmosphere, Part II: Radiative transfer modeling for remote sensing applications. *Journal of Quantitative Spectroscopy and Radiative Transfer*, 311, 108791. https://doi.org/10.1016/j.jqsrt.2023.108791

## Spectroscopic Data (HITRAN and related)

vSmartMOM computes gas absorption from [HITRAN](https://hitran.org) line
parameters, through [AtmosphericAbsorption.jl](https://github.com/RemoteSensingTools/AtmosphericAbsorption.jl).
The spectroscopic parameters behind every line-by-line cross-section come from
HITRAN, so **if you publish results computed with HITRAN data, please cite the
HITRAN edition you used**, together with any of the data sources below that
your calculation used. BibTeX entries are in `CITATION.bib`.

### Which HITRAN edition did I use?

- **Default (`"artifact"` edition):** HITRAN2016, distributed as Julia
  artifacts. This is what `artifact("CO2")` and YAML configurations use unless
  you change the edition.
- **Downloaded data (`set_hitran_edition!`, `fetch_hitran`):** hitran.org serves
  only the *current* edition — HITRAN2024 at the time of writing. The edition
  string is a cache and provenance label; it does not select an older edition
  from the server. `hitran_info(molecule)` reports the download date and
  checksum recorded with the data.
- **Your own `.par` files:** cite the edition the file came from.

### HITRAN

- Gordon, I. E., Rothman, L. S., Hargreaves, R. J., et al. (2026). The HITRAN2024
  molecular spectroscopic database. *Journal of Quantitative Spectroscopy and
  Radiative Transfer*, 353, 109807. https://doi.org/10.1016/j.jqsrt.2026.109807
- Gordon, I. E., Rothman, L. S., Hargreaves, R. J., et al. (2022). The HITRAN2020
  molecular spectroscopic database. *Journal of Quantitative Spectroscopy and
  Radiative Transfer*, 277, 107949. https://doi.org/10.1016/j.jqsrt.2021.107949
- Gordon, I. E., Rothman, L. S., Hill, C., et al. (2017). The HITRAN2016 molecular
  spectroscopic database. *Journal of Quantitative Spectroscopy and Radiative
  Transfer*, 203, 3-69. https://doi.org/10.1016/j.jqsrt.2017.06.038

### Partition sums

- Gamache, R. R., Vispoel, B., Rey, M., et al. (2021). Total internal partition sums
  for the HITRAN2020 database. *Journal of Quantitative Spectroscopy and Radiative
  Transfer*, 271, 107713. https://doi.org/10.1016/j.jqsrt.2021.107713 (TIPS-2021,
  used by AtmosphericAbsorption.jl and hence by the RT pipeline.)
- Gamache, R. R., Roller, C., Lopes, E., et al. (2017). Total internal partition sums
  for 166 isotopologues of 51 molecules important in planetary atmospheres:
  Application to HITRAN2016 and beyond. *Journal of Quantitative Spectroscopy and
  Radiative Transfer*, 203, 70-87. https://doi.org/10.1016/j.jqsrt.2017.03.045
  (TIPS-2017, used by the deprecated internal `vSmartMOM.Absorption` line engine.)

### HAPI, collision-induced absorption and the water-vapour continuum

- Kochanov, R. V., Gordon, I. E., Rothman, L. S., et al. (2016). HITRAN Application
  Programming Interface (HAPI): A comprehensive approach to working with
  spectroscopic data. *Journal of Quantitative Spectroscopy and Radiative
  Transfer*, 177, 15-30. https://doi.org/10.1016/j.jqsrt.2016.03.005 (Cite when
  using non-Voigt line profiles or parameters through AtmosphericAbsorption.jl.)
- Karman, T., Gordon, I. E., van der Avoird, A., et al. (2019). Update of the HITRAN
  collision-induced absorption section. *Icarus*, 328, 160-175.
  https://doi.org/10.1016/j.icarus.2019.02.034 (Cite when using HITRAN `.cia` files.)
- Mlawer, E. J., Payne, V. H., Moncet, J.-L., et al. (2012). Development and recent
  evaluation of the MT_CKD model of continuum absorption. *Philosophical
  Transactions of the Royal Society A*, 370, 2520-2556.
  https://doi.org/10.1098/rsta.2011.0295
- Mlawer, E. J., Cady-Pereira, K. E., Mascio, J., and Gordon, I. E. (2023). The
  inclusion of the MT_CKD water vapor continuum model in the HITRAN molecular
  spectroscopic database. *Journal of Quantitative Spectroscopy and Radiative
  Transfer*, 306, 108645. https://doi.org/10.1016/j.jqsrt.2023.108645
  (Cite both when using the MT_CKD water-vapour continuum.)

ExoMol line lists used through AtmosphericAbsorption.jl have their own citation
requirements; see that package's
[citing page](https://github.com/RemoteSensingTools/AtmosphericAbsorption.jl/blob/main/docs/src/citing.md).

## Background References

- Fell, F. (1997). *Validierung eines Modells zur Simulation des Strahlungstransportes in Atmosphäre und Ozean*.
- Grant, I. P. and Hunt, G. E. (1969). Discrete space theory of radiative transfer I. Fundamentals. *Proceedings of the Royal Society of London A*, 313(1513), 183-197.
- Plass, G. N., Kattawar, G. W., and Catchings, F. E. (1973). Matrix operator theory of radiative transfer. 1: Rayleigh scattering. *Applied Optics*, 12(2), 314-329.
