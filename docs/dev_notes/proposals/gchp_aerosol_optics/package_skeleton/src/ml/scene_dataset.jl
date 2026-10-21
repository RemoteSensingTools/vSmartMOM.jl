# Stage 6 — ML training dataset assembly (brief §8). Generate many scenes → (input state,
# output radiances). Design the schema BEFORE generating ~10^6 scenes.

"""
    SceneDataset

Container/spec for the ML training set. Decide and freeze:
  - OUTPUT: TOA reflectances/radiances over (wavelength grid × geometry sweep [× Stokes]).
            Store the wavelength and geometry grids ONCE; per-scene radiances are the payload.
  - INPUT:  the atmospheric + aerosol state (see `eof_reduce` in reduction.jl for the
            dimensionality problem). Store the EOF basis WITH the dataset so inputs are
            reversible/interpretable.
  - PROVENANCE: GCHP file, date, (x,y,face) indices, scheme, mixing rule, representation.
"""
struct SceneDataset
    # TODO: input spec, output spec (λ grid, geometry grid), basis handle, provenance
end

"""
    build_dataset(gchp_paths, base_params, sweep; rep, mixing, ri_db, reducer) -> SceneDataset

Loop GCHP columns (use `scenes(GCHPFile(...))`), build per-column RT setup
(`scene_to_rt_parameters`), run `vSmartMOM.rt_run` over the geometry/wavelength sweep,
reduce the input profiles (`eof_reduce`), and append (reduced-input, radiance-output) rows.
Stream to disk; do not hold 10^6 scenes in memory.
"""
function build_dataset end   # TODO (new — Stage 6)
