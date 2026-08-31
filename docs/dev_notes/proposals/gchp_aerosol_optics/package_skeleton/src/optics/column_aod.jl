# Stage 1 (EXISTS on gchp-io): bins → column AOD. This is the validation target before the
# RT bridge. PORT TARGET: `compute_column_aod` in gchp-io:src/Aerosols/aod_diagnostics.jl
# and the bulk grid writer in gchp-io:src/IO/Benchmark/aod_diagnostic.jl.

"""
    column_aod(scene, wavelengths_um; rep=Sectional(), mixing=:volume, ri_db=...) -> Vector

Aerosol optical depth per wavelength for one column (extinction only; no phase matrix).
Loop layers → bins → effective RI → Mie extinction cross-section → integrate over the
size grid → sum × layer thickness. Use this to validate the speciation→RI→Mie chain
against GCHP's own AOD diagnostic before building per-layer optics.
"""
function column_aod end   # TODO (port from gchp-io)

"""
    write_grid_aod(out_path, gchp_path; wavelengths_um=[0.76, 1.6, 2.2], aerosol_scheme=:auto)

Open a GCHP file once, loop all (x, y, face) cells, write AOD(x, y, face, λ) to NetCDF.
PORT: gchp-io:src/IO/Benchmark/aod_diagnostic.jl (handles dimension ordering / chunking).
"""
function write_grid_aod end   # TODO (port from gchp-io)
