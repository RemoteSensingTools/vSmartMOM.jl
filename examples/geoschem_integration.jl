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
- `spec_bands`: Spectral bands (e.g., [(1e7/777):0.015:(1e7/757)])
- `surface`: Surface BRDF (e.g., [LambertianSurfaceScalar(0.15)])
- `sza`, `vza`, `vaz`: Geometry angles
- `obs_alt`: Observer altitude
- `architecture`: CPU() or GPU()
- Additional parameters to override defaults

# Returns
- Radiative transfer results

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
    spec_bands = [(1e7/777):0.015:(1e7/757)],
    sza = 45.0,
    vza = [0.0, 30.0, 60.0],
    architecture = GPU()
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
    
    # 2. Load default parameters as a baseline
    default_params = default_parameters()
    
    # 3. Read GEOSChem data and convert to parameter Dict
    geoschem_dict = geoschem_to_dict(src)
    
    # 4. Merge with defaults (user can override radiative_transfer and geometry)
    # Start with defaults, then overlay GEOSChem atmospheric profile and absorption
    merged_dict = Dict{String, Any}()
    
    # Keep existing radiative_transfer section from defaults
    merged_dict["radiative_transfer"] = Dict{String, Any}(
        "spec_bands" => isnothing(spec_bands) ? 
            default_params.spec_bands : 
            [collect(sb) for sb in spec_bands],
        "surface" => isnothing(surface) ? 
            default_params.BRDF : 
            surface,
        "quadrature_type" => isnothing(quadrature_type) ? 
            "$(typeof(default_params.quadrature_type))()" : 
            quadrature_type,
        "polarization_type" => isnothing(polarization_type) ? 
            "$(typeof(default_params.polarization_type))()" : 
            polarization_type,
        "max_m" => default_params.max_m,
        "Δ_angle" => default_params.Δ_angle,
        "l_trunc" => default_params.l_trunc,
        "depol" => default_params.depol,
        "float_type" => "$(default_params.FT)",
        "architecture" => isnothing(architecture) ? 
            "default_architecture" : 
            "$(architecture)"
    )
    
    # Geometry from user or defaults
    merged_dict["geometry"] = Dict{String, Any}(
        "sza" => isnothing(sza) ? default_params.sza : sza,
        "vza" => isnothing(vza) ? default_params.vza : vza,
        "vaz" => isnothing(vaz) ? default_params.vaz : vaz,
        "obs_alt" => isnothing(obs_alt) ? default_params.obs_alt : obs_alt
    )
    
    # Use GEOSChem atmospheric profile and absorption
    merged_dict["atmospheric_profile"] = geoschem_dict["atmospheric_profile"]
    if haskey(geoschem_dict, "absorption")
        merged_dict["absorption"] = geoschem_dict["absorption"]
    end
    
    # Apply any additional kwargs overrides
    for (key, value) in kwargs
        if haskey(merged_dict, String(key))
            merged_dict[String(key)] = value
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
println("Available molecules: ", config["absorption"]["molecules"])

# Modify before using
config["radiative_transfer"]["spec_bands"] = [(1e7/777):0.015:(1e7/757)]
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
#     spec_bands = [(1e7/777):0.015:(1e7/757)],
#     sza = 60.0,
#     vza = [0.0, 30.0, 60.0],
#     vaz = [0.0, 0.0, 0.0],
#     architecture = GPU()
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
