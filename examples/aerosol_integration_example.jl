"""
Example: Complete aerosol integration workflow

Demonstrates how to read aerosol data, compute optical properties,
and prepare for vSmartMOM RT calculations.
"""

using vSmartMOM
using vSmartMOM.Aerosols

# ==============================================================================
# Example 1: TOMAS-15 Size-Resolved Aerosols
# ==============================================================================

function example_tomas15()
    println("\n" * "="^70)
    println("Example 1: TOMAS-15 Size-Resolved Aerosols")
    println("="^70 * "\n")
    
    # File paths
    config_file = "examples/aerosol_config_tomas15.yaml"
    netcdf_file = "GEOSChem.Custom.20190702_0000z.nc4"
    ri_database = "data/refractive_indices_database.yaml"
    
    # Check if files exist
    if !isfile(netcdf_file)
        @warn "NetCDF file not found: $netcdf_file"
        @warn "Skipping TOMAS-15 example"
        return
    end
    
    # Step 1: Read aerosol data
    println("Step 1: Reading TOMAS-15 aerosol data...")
    data = read_aerosol_data(config_file, netcdf_file)
    
    println("  Scheme: $(typeof(data.scheme))")
    println("  Species: $(data.scheme.species)")
    println("  Size bins: $(data.scheme.n_bins)")
    println("  Diameter range: $(data.scheme.diam_min) - $(data.scheme.diam_max) nm")
    println("  Vertical levels: $(length(data.coordinates["lev"]))")
    
    # Step 2: Load refractive index database
    println("\nStep 2: Loading refractive index database...")
    ri_db = load_refractive_index_database(ri_database)
    
    println("  Available species: $(list_species(ri_db))")
    
    # Step 3: Examine a specific species
    println("\nStep 3: Examining DUST aerosol...")
    dust_data = data.species_data["DUST"]
    dust_conc = dust_data.data["concentration"]  # [15 bins × 72 levels]
    
    println("  Concentration array size: $(size(dust_conc))")
    println("  Max concentration: $(maximum(dust_conc)) mol/mol")
    
    # Find level with maximum dust
    _, max_idx = findmax(vec(sum(dust_conc, dims=1)))
    println("  Peak dust level: $max_idx")
    
    # Step 4: Get refractive index
    println("\nStep 4: Retrieving refractive indices...")
    ri_key = data.scheme.refractive_indices["DUST"]
    println("  DUST uses RI key: $ri_key")
    
    wavelengths = [0.4, 0.55, 0.86, 1.6]  # UV, visible, NIR, SWIR
    n_dust = get_refractive_index(ri_db, ri_key, wavelengths)
    
    for (λ, n) in zip(wavelengths, n_dust)
        println("  λ = $λ μm: n = $(real(n)) + $(imag(n))i")
    end
    
    # Step 5: Compute optical properties
    println("\nStep 5: Computing optical properties...")
    opt_props = compute_optical_properties(data, wavelengths, ri_db)
    
    println("  Extinction: $(size(opt_props["extinction"]))")
    println("  Max extinction at 550 nm: $(maximum(opt_props["extinction"][:, 2])) km⁻¹")
    println("  Mean SSA at 550 nm: $(mean(opt_props["ssa"][:, 2]))")
    println("  Mean asymmetry parameter: $(mean(opt_props["asymmetry_parameter"][:, 2]))")
    
    # Step 6: Analyze spectral dependence
    println("\nStep 6: Spectral analysis at peak dust level...")
    level = max_idx
    
    println("  λ (μm) | Extinction (km⁻¹) | SSA | g")
    println("  " * "-"^50)
    for (i, λ) in enumerate(wavelengths)
        ext = opt_props["extinction"][level, i]
        ssa = opt_props["ssa"][level, i]
        g = opt_props["asymmetry_parameter"][level, i]
        println("  $(lpad(round(λ, digits=2), 5)) | $(lpad(round(ext, digits=4), 15)) | " *
                "$(lpad(round(ssa, digits=3), 3)) | $(lpad(round(g, digits=3), 3))")
    end
    
    println("\n✓ TOMAS-15 example completed successfully!")
    
    return data, opt_props, ri_db
end

# ==============================================================================
# Example 2: Two-Moment Bulk Aerosols
# ==============================================================================

function example_two_moment()
    println("\n" * "="^70)
    println("Example 2: Two-Moment Bulk Aerosols")
    println("="^70 * "\n")
    
    # File paths
    config_file = "examples/aerosol_config_two_moment.yaml"
    netcdf_file = "GEOSChem.Aerosols.20190702_0000z.nc4"
    ri_database = "data/refractive_indices_database.yaml"
    
    # Check if files exist
    if !isfile(netcdf_file)
        @warn "NetCDF file not found: $netcdf_file"
        @warn "This is expected - two-moment file not provided"
        @warn "Skipping two-moment example"
        return
    end
    
    # Step 1: Read aerosol data
    println("Step 1: Reading two-moment aerosol data...")
    data = read_aerosol_data(config_file, netcdf_file)
    
    println("  Scheme: $(typeof(data.scheme))")
    println("  Species: $(data.scheme.species)")
    println("  Vertical levels: $(length(data.coordinates["lev"]))")
    
    # Step 2: Load refractive index database
    println("\nStep 2: Loading refractive index database...")
    ri_db = load_refractive_index_database(ri_database)
    
    # Step 3: Examine sulfate aerosol
    println("\nStep 3: Examining sulfate aerosol...")
    so4_data = data.species_data["so4"]
    so4_aod = so4_data.data["aod"]
    so4_radius = so4_data.data["radius"]
    
    println("  AOD array size: $(size(so4_aod))")
    println("  Max AOD: $(maximum(so4_aod))")
    println("  Mean effective radius: $(mean(so4_radius)) μm")
    println("  Geometric std dev: $(data.scheme.sigma_g["so4"])")
    
    # Step 4: Compute optical properties
    println("\nStep 4: Computing optical properties...")
    wavelengths = [0.4, 0.55, 0.86, 1.6]
    opt_props = compute_optical_properties(data, wavelengths, ri_db)
    
    println("  Total AOD at 550 nm: $(sum(opt_props["extinction"][:, 2]))")
    println("  Mean SSA at 550 nm: $(mean(opt_props["ssa"][:, 2]))")
    
    println("\n✓ Two-moment example completed successfully!")
    
    return data, opt_props, ri_db
end

# ==============================================================================
# Example 3: Advanced Usage - Size Distribution Analysis
# ==============================================================================

function example_size_distribution_analysis()
    println("\n" * "="^70)
    println("Example 3: Size Distribution Analysis")
    println("="^70 * "\n")
    
    config_file = "examples/aerosol_config_tomas15.yaml"
    netcdf_file = "GEOSChem.Custom.20190702_0000z.nc4"
    
    if !isfile(netcdf_file)
        @warn "Skipping size distribution example"
        return
    end
    
    # Read data
    data = read_aerosol_data(config_file, netcdf_file)
    
    # Select a level (e.g., 800 hPa ~ level 20)
    level = 20
    
    println("Analysis at vertical level $level\n")
    
    # Get bin centers and widths
    bin_centers = data.scheme.bin_centers  # nm
    bin_edges = data.scheme.bin_edges
    bin_widths = diff(bin_edges)
    
    # Compute dN/dlogD for each species
    println("Species | Total Concentration | Mode Diameter")
    println("-"^55)
    
    for species_name in data.scheme.species
        conc = data.species_data[species_name].data["concentration"][:, level]
        
        # Convert to dN/dlogD
        dlog_D = log10.(bin_edges[2:end] ./ bin_edges[1:end-1])
        dN_dlogD = conc ./ dlog_D
        
        # Find mode (peak)
        max_val, max_idx = findmax(dN_dlogD)
        mode_diameter = bin_centers[max_idx]
        
        total_conc = sum(conc)
        
        println("$(rpad(species_name, 7)) | " *
                "$(lpad(round(total_conc, sigdigits=3), 18)) | " *
                "$(lpad(round(mode_diameter, digits=1), 14)) nm")
    end
    
    println("\n✓ Size distribution analysis completed!")
end

# ==============================================================================
# Example 4: Refractive Index Database Exploration
# ==============================================================================

function example_ri_database()
    println("\n" * "="^70)
    println("Example 4: Refractive Index Database Exploration")
    println("="^70 * "\n")
    
    ri_database = "data/refractive_indices_database.yaml"
    
    # Load database
    ri_db = load_refractive_index_database(ri_database)
    
    # Show full database info
    show_database_info(ri_db)
    
    # Compare species at 550 nm
    println("\nComparison at λ = 550 nm:")
    println("="^60)
    println("Species           | n_real | n_imag   | Type")
    println("-"^60)
    
    species_list = list_species(ri_db)
    for species in species_list
        n = get_refractive_index(ri_db, species, 0.55)
        
        # Classify by absorption
        if imag(n) < 1e-6
            type_str = "Non-absorbing"
        elseif imag(n) < 0.01
            type_str = "Weakly absorbing"
        elseif imag(n) < 0.1
            type_str = "Moderately absorbing"
        else
            type_str = "Strongly absorbing"
        end
        
        println("$(rpad(species, 17)) | " *
                "$(lpad(round(real(n), digits=4), 6)) | " *
                "$(lpad(round(imag(n), digits=6), 8)) | $type_str")
    end
    
    # Spectral variation for black carbon
    println("\n\nSpectral variation for black carbon:")
    println("="^60)
    wavelengths = [0.3, 0.4, 0.55, 0.86, 1.6, 2.5]
    
    println("λ (μm) | n_real | n_imag")
    println("-"^35)
    for λ in wavelengths
        try
            n = get_refractive_index(ri_db, "black_carbon", λ)
            println("$(lpad(round(λ, digits=2), 6)) | " *
                    "$(lpad(round(real(n), digits=4), 6)) | " *
                    "$(lpad(round(imag(n), digits=4), 6))")
        catch e
            println("$(lpad(round(λ, digits=2), 6)) | (out of range)")
        end
    end
    
    println("\n✓ Database exploration completed!")
end

# ==============================================================================
# Main execution
# ==============================================================================

function run_all_examples()
    println("\n")
    println("╔" * "="^68 * "╗")
    println("║" * " "^10 * "vSmartMOM Aerosol Framework Examples" * " "^20 * "║")
    println("╚" * "="^68 * "╝")
    
    try
        # Example 1: TOMAS-15
        example_tomas15()
    catch e
        @error "Error in TOMAS-15 example" exception=e
    end
    
    try
        # Example 2: Two-moment
        example_two_moment()
    catch e
        @error "Error in two-moment example" exception=e
    end
    
    try
        # Example 3: Size distribution
        example_size_distribution_analysis()
    catch e
        @error "Error in size distribution example" exception=e
    end
    
    try
        # Example 4: RI database
        example_ri_database()
    catch e
        @error "Error in RI database example" exception=e
    end
    
    println("\n")
    println("╔" * "="^68 * "╗")
    println("║" * " "^20 * "All Examples Completed!" * " "^26 * "║")
    println("╚" * "="^68 * "╝")
    println()
end

# Run examples if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    run_all_examples()
end
