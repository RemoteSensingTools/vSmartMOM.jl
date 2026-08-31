# Three-aerosol OCO XCO2 experiment

This workflow defines OCO-like O2 A (757–773 nm), weak CO2 (1589–1622 nm),
and strong CO2 (2042–2084 nm) bands with three Mie aerosol modes:

- sulfate: median radius 0.10 µm, centered at 1.0 km;
- organic carbon: median radius 0.20 µm, centered at 1.2 km;
- background stratospheric sulfate: median radius 0.15 µm, centered at 12 km.

The assumed AODs, widths, and refractive indices are explicit in
`config/oco_grass_3aerosol.yaml`; they are scenario assumptions, not retrieved
site-specific values. Each band uses a Lambertian grassy-vegetation spectrum
represented by P0, P1, and P2. `scripts/fit_grass_surface.jl` reproduces the
coefficients from three representative reflectance anchors per band.

Set spectroscopy paths individually or set their common directory:

```bash
export VSMARTMOM_HITRAN_LUT_DIR=$HOME/data/HITRAN_LUTs
julia --project=. RRS_XCO2/scripts/fit_grass_surface.jl
julia --project=. RRS_XCO2/scripts/run_forward.jl
python RRS_XCO2/scripts/plot_results.py
```

The basis grid is uniformly spaced at 0.1 cm^-1 and the atmospheric profile is
reduced to 12 layers. Forward calculations use external-solar SFI: with five
weighted streams, one VZA, and IQU polarization the diffuse operator is 18×18;
the solar direction is retained only on the phase-source grid.

On wurst, benchmark physical CUDA device 1 and generate the native-grid plots:

```bash
bash RRS_XCO2/scripts/run_wurst_gpu1.sh
```

The benchmark records cold and repeated warm wall-clock timings and native
basis spectra. Set `RRS_XCO2_OUTPUT_DIR` to keep results from different hosts
separate. `plot_basis_results.py` writes one forward figure and
one numbered figure for every retrieval-state row of the available Stokes-I
Jacobian. Jacobian figures are collected under
`output/basis_jacobian_three_bands/`; filenames follow
`jacobian_row_NNN__parameter_name.png`, and titles state the parameter in
`∂I/∂x`. Every figure has vertically stacked O2 A, weak-CO2, and strong-CO2
panels showing the same physical state parameter in every band, including a
zero curve where a band has no sensitivity. Gas profiles are shared retrieval
parameters: one `co2_layer_01` row contains the O2-A, weak-CO2, and strong-CO2
responses rather than three band-specific CO2 rows. Instrument convolution and
OCO-2 resampling will be added after
an L1B calibration file is supplied.

## Linearization status

The shared state contains surface pressure, 21 aerosol columns (seven per
species), 12 H2O columns, 12 CO2 columns, and three Legendre surface columns.
Each gas column contributes to all three band-output blocks. The
forward and analytic-linearized calculations both use external-solar SFI, so
the diffuse IQU operator remains 18×18 and the direct solar coupling is carried
by rectangular source-column operators.
