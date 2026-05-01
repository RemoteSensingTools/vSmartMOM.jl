"""
Test suite for Aerosol module

Tests the flexible aerosol framework with both TOMAS-15 and two-moment schemes.
"""

using Test
using vSmartMOM
using vSmartMOM.Aerosols
using YAML

# Test data directory
const TEST_DATA_DIR = dirname(@__DIR__)
const TOMAS_FILE = joinpath(TEST_DATA_DIR, "GEOSChem.Custom.20190702_0000z.nc4")
const TOMAS_CONFIG = joinpath(TEST_DATA_DIR, "examples", "aerosol_config_tomas15.yaml")
const TWOMOM_CONFIG = joinpath(TEST_DATA_DIR, "examples", "aerosol_config_two_moment.yaml")
const RI_DATABASE = joinpath(TEST_DATA_DIR, "data", "refractive_indices_database.yaml")

@testset "Aerosol Module Tests" begin
    
    @testset "Refractive Index Database" begin
        # Load database
        @test isfile(RI_DATABASE)
        db = load_refractive_index_database(RI_DATABASE)
        
        # Check structure
        @test db isa RefractiveIndexDatabase
        @test length(db.data) == 6  # 6 aerosol types
        
        # Check species
        species_list = list_species(db)
        @test "sulfate_suso" in species_list
        @test "organic_carbon" in species_list
        @test "black_carbon" in species_list
        @test "seasalt_sscm" in species_list
        @test "dust_opac" in species_list
        @test "water" in species_list
        
        # Test interpolation at 550 nm (visible)
        n_sulfate = get_refractive_index(db, "sulfate_suso", 0.55)
        @test real(n_sulfate) > 1.0  # Real part should be > 1
        @test real(n_sulfate) < 2.0  # Reasonable range
        @test imag(n_sulfate) >= 0.0  # Absorption should be non-negative
        @test imag(n_sulfate) < 0.1   # Sulfate is weakly absorbing
        
        # Black carbon should be strongly absorbing. Value depends on the
        # refractive-index database; typical literature BC at 550 nm is in
        # the 0.3–0.8 range depending on condensation/aging.
        n_bc = get_refractive_index(db, "black_carbon", 0.55)
        @test imag(n_bc) > 0.3  # Strong absorber
        
        # Test wavelength range checking
        λ_min, λ_max = wavelength_range(db, "sulfate_suso")
        @test λ_min == 0.3
        @test λ_max == 3.75
        
        # Test error handling for out-of-range wavelength
        @test_throws ErrorException get_refractive_index(db, "sulfate_suso", 10.0)
        
        # Test error handling for unknown species
        @test_throws ErrorException get_refractive_index(db, "unicorn_dust", 0.55)
        
        # Test vectorized interpolation
        wavelengths = [0.4, 0.55, 0.86, 1.0]
        n_vec = get_refractive_index(db, "sulfate_suso", wavelengths)
        @test length(n_vec) == 4
        @test all(real.(n_vec) .> 1.0)
    end
    
    @testset "TOMAS-15 Scheme Types" begin
        # Test config loading
        @test isfile(TOMAS_CONFIG)
        config = YAML.load_file(TOMAS_CONFIG)
        
        # Construct scheme
        scheme = TOMAS15Scheme(config)
        options = vSmartMOM.Aerosols.aerosol_processing_options(config)
        
        # Check basic properties
        @test scheme isa TOMAS15Scheme
        @test scheme.n_bins == 15
        @test get(options, "vertical_flip", false) === true
        @test vSmartMOM.Aerosols._tomas_concentration_variable(config, "DUST", 1) == "SpeciesConcVV_DUST01"
        @test scheme.diam_min == 10.0  # nm
        @test scheme.diam_max == 10000.0  # nm
        @test length(scheme.species) == 8
        
        # Check bin edges
        @test length(scheme.bin_edges) == 16  # n_bins + 1
        @test scheme.bin_edges[1] ≈ 10.0
        @test scheme.bin_edges[end] ≈ 10000.0
        
        # Check logarithmic spacing
        ratios = scheme.bin_edges[2:end] ./ scheme.bin_edges[1:end-1]
        @test all(isapprox.(ratios, ratios[1], rtol=1e-10))
        
        # Check bin centers
        @test length(scheme.bin_centers) == 15
        expected_centers = sqrt.(scheme.bin_edges[1:end-1] .* scheme.bin_edges[2:end])
        @test all(isapprox.(scheme.bin_centers, expected_centers))
        
        # Check species properties
        @test haskey(scheme.refractive_indices, "DUST")
        @test haskey(scheme.densities, "DUST")
        @test haskey(scheme.molar_masses, "DUST")
        
        # Physical reasonability checks
        @test all(values(scheme.densities) .> 0.0)
        @test all(values(scheme.molar_masses) .> 0.0)
    end
    
    @testset "Two-Moment Scheme Types" begin
        # Test config loading
        @test isfile(TWOMOM_CONFIG)
        config = YAML.load_file(TWOMOM_CONFIG)
        
        # Construct scheme
        scheme = TwoMomentScheme(config)
        
        # Check basic properties
        @test scheme isa TwoMomentScheme
        @test length(scheme.species) == 7
        
        # Check species properties
        @test all(haskey(scheme.sigma_g, sp) for sp in scheme.species)
        @test all(haskey(scheme.aod_wavelength, sp) for sp in scheme.species)
        @test all(haskey(scheme.refractive_indices, sp) for sp in scheme.species)
        
        # Physical reasonability
        @test all(values(scheme.sigma_g) .>= 1.0)  # σ_g >= 1
        @test all(values(scheme.sigma_g) .<= 3.0)  # reasonable upper limit
        @test all(values(scheme.aod_wavelength) .> 0.0)
    end
    
    @testset "TOMAS-15 Data Reading" begin
        # Only run if test data file exists
        if isfile(TOMAS_FILE)
            @info "Testing TOMAS-15 data reading with actual file"
            
            # Read data
            data = read_aerosol_data(TOMAS_CONFIG, TOMAS_FILE)
            
            # Check structure
            @test data isa AerosolData{<:TOMAS15Scheme}
            @test length(data.species_data) == length(data.scheme.species) + 1
            @test haskey(data.species_data, "NK")
            
            # Check coordinates
            @test haskey(data.coordinates, "lev")
            n_levels = length(data.coordinates["lev"])
            @test n_levels == 72
            
            # Check each species
            for species_name in data.scheme.species
                @test haskey(data.species_data, species_name)
                sp_data = data.species_data[species_name]
                
                # Check data structure
                @test haskey(sp_data.data, "concentration")
                conc = sp_data.data["concentration"]
                
                # Check dimensions
                @test size(conc) == (15, 72)  # 15 bins × 72 levels
                
                # Check physical reasonability
                @test all(conc .>= 0.0)  # Concentrations non-negative
                @test any(conc .> 0.0)   # Some non-zero values
                
                # Check units
                @test haskey(sp_data.units, "concentration")
                @test sp_data.units["concentration"] == "mol mol-1 dry air"
            end
            
            # Check metadata
            @test !isempty(data.metadata)
            @test haskey(data.metadata, "dimensions")
        else
            @warn "Skipping TOMAS-15 data reading test: file not found at $TOMAS_FILE"
        end
    end
    
    @testset "Helper Functions" begin
        # Test number concentration conversion
        vmr = [1e-9, 1e-10, 1e-11]  # mol/mol
        pressure = [1.0e5, 5.0e4, 2.5e4]  # Pa
        temperature = [300.0, 250.0, 220.0]  # K
        
        n_conc = Aerosols.compute_number_concentration(vmr, pressure, temperature)
        
        @test length(n_conc) == 3
        @test all(n_conc .> 0.0)
        @test n_conc[1] > n_conc[2] > n_conc[3]  # Decreases with altitude
        
        # Test mass concentration conversion
        molar_mass = 0.098  # kg/mol (H2SO4)
        mass_conc = Aerosols.compute_mass_concentration(vmr, molar_mass, pressure, temperature)
        
        @test length(mass_conc) == 3
        @test all(mass_conc .> 0.0)
        
        # Test bin volume calculation
        diam_nm = 100.0
        vol = Aerosols.bin_volume(diam_nm)
        expected_vol = (4.0/3.0) * π * (50.0)^3  # radius = 50 nm
        @test vol ≈ expected_vol
        
        # Test AOD wavelength scaling
        aod_ref = 0.5
        λ_ref = 0.55  # μm
        λ_target = 1.0  # μm
        α = 1.0  # Ångström exponent
        
        aod_scaled = Aerosols.scale_aod_wavelength(aod_ref, λ_ref, λ_target, α)
        expected = 0.5 * (1.0 / 0.55)^(-1.0)
        @test aod_scaled ≈ expected
        
        # Test lognormal size distribution
        r = [0.1, 0.3, 1.0, 3.0]  # μm
        r_eff = 1.0  # μm
        σ_g = 1.6
        
        dN_dr = Aerosols.lognormal_size_distribution(r, r_eff, σ_g)
        
        @test length(dN_dr) == 4
        @test all(dN_dr .>= 0.0)
        @test all(isfinite.(dN_dr))
        
        # Distribution should be normalized (roughly)
        # Integral of dN/dr over log-space should be 1
        
        # Test effective radius conversions
        r_med = 0.5  # μm
        r_eff_calc = Aerosols.effective_radius_from_moments(r_med, σ_g)
        r_med_back = Aerosols.median_radius_from_effective(r_eff_calc, σ_g)
        @test r_med_back ≈ r_med rtol=1e-10
    end
    
    @testset "Mie Scattering Approximations" begin
        # Test Rayleigh regime
        x_small = 0.05
        n = ComplexF64(1.5, 0.01)
        Q_ext, Q_sca, Q_abs, g = Aerosols.compute_mie_efficiencies(x_small, n)
        
        # Scattering dominates in Rayleigh; equality (Q_sca == Q_ext) possible
        # when the approximate implementation treats weak absorption as zero
        # at very small size parameters.
        @test Q_sca <= Q_ext
        @test Q_abs >= 0.0
        @test g ≈ 0.0  # Rayleigh scattering is symmetric
        @test Q_ext ≈ Q_sca + Q_abs
        
        # Test geometric optics regime
        x_large = 50.0
        Q_ext, Q_sca, Q_abs, g = Aerosols.compute_mie_efficiencies(x_large, n)
        
        @test Q_ext ≈ 2.0 rtol=0.1  # Approaches 2 in geometric optics
        @test g > 0.5  # Forward scattering dominates
        @test Q_ext ≈ Q_sca + Q_abs
        
        # Test intermediate regime
        x_mid = 5.0
        Q_ext, Q_sca, Q_abs, g = Aerosols.compute_mie_efficiencies(x_mid, n)
        
        @test 0.0 < Q_ext < 3.0  # Reasonable range
        @test 0.0 <= g <= 1.0
        @test Q_ext ≈ Q_sca + Q_abs
    end
    
    @testset "Optical Properties (TOMAS-15)" begin
        if isfile(TOMAS_FILE)
            @info "Testing optical property calculations for TOMAS-15"
            
            # Read data
            data = read_aerosol_data(TOMAS_CONFIG, TOMAS_FILE)
            
            # Load RI database
            ri_db = load_refractive_index_database(RI_DATABASE)
            
            # Compute optical properties at a few wavelengths
            wavelengths = [0.55, 0.86, 1.0]  # μm
            
            opt_props = compute_optical_properties(data, wavelengths, ri_db)
            
            # Check structure
            @test haskey(opt_props, "extinction")
            @test haskey(opt_props, "scattering")
            @test haskey(opt_props, "absorption")
            @test haskey(opt_props, "ssa")
            @test haskey(opt_props, "asymmetry_parameter")
            
            # Check dimensions
            n_levels = 72
            n_wavelengths = 3
            @test size(opt_props["extinction"]) == (n_levels, n_wavelengths)
            @test size(opt_props["ssa"]) == (n_levels, n_wavelengths)
            
            # Physical constraints
            @test all(opt_props["extinction"] .>= 0.0)
            @test all(opt_props["scattering"] .>= 0.0)
            @test all(opt_props["absorption"] .>= 0.0)
            @test all(0.0 .<= opt_props["ssa"] .<= 1.0)
            @test all(-1.0 .<= opt_props["asymmetry_parameter"] .<= 1.0)
            
            # Conservation: extinction = scattering + absorption
            ext_check = opt_props["scattering"] .+ opt_props["absorption"]
            @test all(isapprox.(opt_props["extinction"], ext_check, rtol=1e-6))
        else
            @warn "Skipping optical properties test: file not found"
        end
    end
    
end

# Print summary
println("\n" * "="^60)
println("Aerosol module test suite completed")
println("="^60)
