"""
Example: Using GEOSChem data with vSmartMOM

This example demonstrates the new clean API for reading GEOSChem NetCDF files
and running radiative transfer calculations.

The new IOSource-based approach provides:
- Type-safe data source specification
- Clean separation of I/O and RT logic
- Easy extensibility for other formats (WRF, GCHP, etc.)
- Full integration with existing vSmartMOM parameter system
"""

using vSmartMOM

# =============================================================================
# Approach 1: Using GeosChemSource directly (Recommended)
# =============================================================================

"""
    run_rt_with_geoschem_v2(file::String, idx::Int, idy::Int, idf::Int; 
                            spec_bands=nothing, kwargs...)

Clean API for running vSmartMOM with GEOSChem data using the new IOSource system.

# Arguments
- `file::String`: Path to GEOSChem NetCDF4 file
- `idx::Int`: X-dimension grid index
- `idy::Int`: Y-dimension grid index
- `idf::Int`: Cubed-sphere face index (1-6)

# Keyword Arguments
- `spec_bands`: Spectral bands as strings (e.g., ["(1e7/777):0.015:(1e7/757)"])
- `surface`: Surface BRDF strings, one per band (e.g., ["LambertianSurfaceScalar(0.15)"])
- `sza`, `vza`, `vaz`: Geometry angles
- `obs_alt`: Observer altitude
- `architecture`: "CPU()", "GPU()" or "default_architecture"
- Additional parameters to override defaults

# Returns
- Radiative transfer results

# Notes
The merged configuration starts from `DefaultParameters.yaml`, so it carries
the default aerosol of that file (one Mie aerosol with `τ_ref = 1.0` at
`λ_ref = 0.770` µm, `p₀ = 900` hPa, `σp = 50` hPa, `nquad_radius = 2500`).
That is a heavy, optically thick choice. To replace it, pass a whole
`scattering` block as a keyword (`scattering = Dict("aerosols" => [...], ...)`,
see the `scattering` schema page); to run without aerosols, pass
`scattering = nothing`, or `delete!(config, "scattering")` on the dict from
`load_geoschem_as_dict` below before parsing it.

# Example
```julia
# Minimal usage - uses mostly defaults
R = run_rt_with_geoschem_v2(
    "GEOSChem.Custom.20190101_0000z.nc4", 
    10, 20, 1
)

# With custom parameters
R = run_rt_with_geoschem_v2(
    "GEOSChem.Custom.20190101_0000z.nc4",
    10, 20, 1;
    spec_bands = ["(1e7/777):0.015:(1e7/757)"],
    sza = 45.0,
    vza = [0.0, 30.0, 60.0],
    vaz = [0.0, 0.0, 0.0],
    architecture = "GPU()"
)
```
"""
function run_rt_with_geoschem_v2(
    file::String, 
    idx::Int, 
    idy::Int, 
    idf::Int;
    # Overrideable parameters
    spec_bands = nothing,
    surface = nothing,
    sza = nothing,
    vza = nothing,
    vaz = nothing,
    obs_alt = nothing,
    architecture = nothing,
    quadrature_type = nothing,
    polarization_type = nothing,
    kwargs...
)
    # 1. Create the data source
    src = GeosChemSource(file, idx, idy, idf)
    
    # 2. Read GEOSChem data as a complete parameter Dict. `geoschem_to_dict`
    #    starts from the package defaults (DefaultParameters.yaml), so the
    #    radiative_transfer, geometry and scattering blocks are already set;
    #    the atmospheric profile and absorption come from the GEOSChem file.
    merged_dict = geoschem_to_dict(src)
    rt = merged_dict["radiative_transfer"]
    geom = merged_dict["geometry"]
    
    # 3. Overlay user choices. The parser reads the radiative_transfer entries
    #    as strings, e.g. spec_bands = ["(1e7/777):0.015:(1e7/757)"] and
    #    surface = ["LambertianSurfaceScalar(0.15)"].
    isnothing(spec_bands)        || (rt["spec_bands"] = spec_bands)
    isnothing(surface)           || (rt["surface"] = surface)
    isnothing(quadrature_type)   || (rt["quadrature_type"] = quadrature_type)
    isnothing(polarization_type) || (rt["polarization_type"] = polarization_type)
    isnothing(architecture)      || (rt["architecture"] = architecture)
    # One surface per band: repeat a single surface if several bands are given.
    n_bands = length(rt["spec_bands"])
    if length(rt["surface"]) == 1 && n_bands > 1
        rt["surface"] = fill(only(rt["surface"]), n_bands)
    end
    # `geoschem_to_dict` writes one-band molecule lists; repeat them (and any
    # other per-band absorption list) so every band has its own entry.
    if haskey(merged_dict, "absorption") && n_bands > 1
        ab = merged_dict["absorption"]
        for key in ("fixed_molecules", "variable_molecules", "molecules",
                    "h2o_absorption", "LUTfiles")
            haskey(ab, key) && length(ab[key]) == 1 || continue
            ab[key] = [deepcopy(only(ab[key])) for _ in 1:n_bands]
        end
    end
    
    # Geometry from user or defaults
    isnothing(sza)     || (geom["sza"] = sza)
    isnothing(vza)     || (geom["vza"] = vza)
    isnothing(vaz)     || (geom["vaz"] = vaz)
    isnothing(obs_alt) || (geom["obs_alt"] = obs_alt)
    
    # 4. Apply any additional kwargs overrides (whole top-level blocks);
    #    `nothing` removes a block, e.g. `scattering = nothing`.
    for (key, value) in kwargs
        if haskey(merged_dict, String(key))
            if isnothing(value)
                delete!(merged_dict, String(key))
            else
                merged_dict[String(key)] = value
            end
        end
    end
    
    # 5. Convert to vSmartMOM_Parameters
    parameters = parameters_from_dict(merged_dict)
    
    # 6. Create model and run RT
    model = model_from_parameters(parameters)
    R = rt_run(model)
    
    return R
end

# =============================================================================
# Approach 2: Direct Dict-based workflow (for advanced users)
# =============================================================================

"""
    load_geoschem_as_dict(file::String, idx::Int, idy::Int, idf::Int) -> Dict

Load GEOSChem data and return as a configuration Dict.
Useful for inspecting or modifying the configuration before creating parameters.

# Example
```julia
# Load and inspect
config = load_geoschem_as_dict("GEOSChem.Custom.20190101_0000z.nc4", 10, 20, 1)
println("Available molecules: ", config["absorption"]["fixed_molecules"])

# Modify before using
config["radiative_transfer"]["spec_bands"] = ["(1e7/777):0.015:(1e7/757)"]
config["geometry"]["sza"] = 45.0

# Create parameters
params = parameters_from_dict(config)
```
"""
function load_geoschem_as_dict(file::String, idx::Int, idy::Int, idf::Int)
    src = GeosChemSource(file, idx, idy, idf)
    return geoschem_to_dict(src)
end

# =============================================================================
# Approach 3: Backwards-compatible with your original function
# =============================================================================

"""
    run_rt_with_geoschem(file::String, idx::Int, idy::Int, idf::Int)

Original API maintained for backwards compatibility.
Now implemented cleanly using the new IOSource system.
"""
function run_rt_with_geoschem(file::String, idx::Int, idy::Int, idf::Int)
    # Simply delegates to the new implementation
    return run_rt_with_geoschem_v2(file, idx, idy, idf)
end

# =============================================================================
# Example Usage Patterns
# =============================================================================

# Pattern 1: Quick and simple (uses defaults)
# R = run_rt_with_geoschem_v2("GEOSChem.Custom.20190101_0000z.nc4", 10, 20, 1)

# Pattern 2: With custom RT parameters
# R = run_rt_with_geoschem_v2(
#     "GEOSChem.Custom.20190101_0000z.nc4", 10, 20, 1;
#     spec_bands = ["(1e7/777):0.015:(1e7/757)"],
#     sza = 60.0,
#     vza = [0.0, 30.0, 60.0],
#     vaz = [0.0, 0.0, 0.0],
#     architecture = "GPU()"
# )

# Pattern 3: Maximum control - modify Dict before conversion
# src = GeosChemSource("GEOSChem.Custom.20190101_0000z.nc4", 10, 20, 1)
# config_dict = geoschem_to_dict(src)
# # ... modify config_dict as needed ...
# params = parameters_from_dict(config_dict)
# model = model_from_parameters(params)
# R = rt_run(model)

# Pattern 4: Integrate with existing YAML workflow
# # Read GEOSChem for atmospheric profile
# geoschem_config = load_geoschem_as_dict("GEOSChem.Custom.20190101_0000z.nc4", 10, 20, 1)
# # Load RT parameters from YAML
# yaml_config = parameters_from_yaml("my_rt_params.yaml")
# # Merge: use YAML for RT settings, GEOSChem for atmosphere
# # ... custom merging logic ...
