#!/usr/bin/env python3
"""
Quick TOMAS Aerosol Data Explorer (Python version)

Explores the GEOSChem-TOMAS NetCDF file and creates visualizations.

Requirements:
    pip install netCDF4 numpy matplotlib scipy

Usage:
    python test/explore_tomas_aerosols.py
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit

def lognormal(r, N_total, r_med, sigma_g):
    """
    Log-normal size distribution (number concentration form)
    
    Parameters:
    -----------
    r : array
        Particle radius [μm]
    N_total : float
        Total concentration (integral of distribution)
    r_med : float
        Median radius [μm]
    sigma_g : float
        Geometric standard deviation (dimensionless)
    
    Returns:
    --------
    dN/dlogr : array
        Number (or mass) per log-radius bin
    """
    return (N_total / (np.sqrt(2*np.pi) * np.log(sigma_g))) * \
           np.exp(-0.5 * (np.log(r/r_med) / np.log(sigma_g))**2)

def fit_lognormal(r_centers, concentrations):
    """
    Fit a log-normal distribution to binned concentration data
    
    Returns (N_total, r_med, sigma_g) or None if fit fails
    """
    # Filter out zero/negative values
    valid = concentrations > 0
    if valid.sum() < 3:  # Need at least 3 points
        return None
    
    r_valid = r_centers[valid]
    c_valid = concentrations[valid]
    
    # Initial guess
    r_med_guess = r_valid[np.argmax(c_valid)]
    N_total_guess = np.trapz(c_valid, np.log(r_valid))
    sigma_g_guess = 2.0
    
    try:
        params, _ = curve_fit(
            lognormal, r_valid, c_valid,
            p0=[N_total_guess, r_med_guess, sigma_g_guess],
            bounds=([0, r_valid.min(), 1.01], [np.inf, r_valid.max(), 5.0]),
            maxfev=5000
        )
        return params
    except:
        return None

# Configuration
ncfile = "GEOSChem.Custom.20190702_0000z.nc4"
output_dir = Path("test/aerosol_exploration_output")
output_dir.mkdir(exist_ok=True)

# Define multiple land locations to analyze
# (idx, idy, face, name)
locations = [
    (2, 12, 4, "Central USA"),        # 36.8°N, 97.9°W (Oklahoma/Kansas)
    (15, 15, 1, "Amazon Basin"),      # Tropical South America
    (8, 18, 2, "Sahara Desert"),      # North Africa
    (12, 8, 3, "India/South Asia"),   # Indian subcontinent
    (5, 5, 5, "Northern China"),      # East Asia
]

print("=" * 70)
print("TOMAS Aerosol Data Exploration (Python)")
print("=" * 70)
print()

# Open dataset
ds = nc.Dataset(ncfile, 'r')

# Get location info for all sites
print("📍 Analyzing {} locations:".format(len(locations)))
for i, (idx_loc, idy_loc, idf_loc, name) in enumerate(locations, 1):
    lat_loc = ds['lats'][idf_loc, idy_loc, idx_loc]
    lon_loc = ds['lons'][idf_loc, idy_loc, idx_loc]
    print(f"   {i}. {name:20s}: {lat_loc:6.2f}°N, {lon_loc:6.2f}°E")
print()

# Use first location for detailed single-column analysis
idx, idy, idf, loc_name = locations[0]
lat = ds['lats'][idf, idy, idx]
lon = ds['lons'][idf, idy, idx]
print(f"📍 Detailed analysis for: {loc_name}")
print(f"   Location: idx={idx+1}, idy={idy+1}, face={idf+1}")
print(f"   Latitude:  {lat:.2f}°")
print(f"   Longitude: {lon:.2f}°")
print()

# TOMAS-15 bin definitions
# IMPORTANT: TOMAS bins are defined by DRY PARTICLE DIAMETER
# Size range: 10 nm to 10 μm diameter (15 logarithmically spaced bins)
# Reference: TOMAS documentation - 15 bins spanning 10 nm to 10 micron dry diameter

# Calculate bin edges: 10 nm to 10 μm diameter (15 bins → 16 edges)
diam_min_nm = 10.0      # 10 nm minimum dry diameter
diam_max_nm = 10000.0   # 10 μm maximum dry diameter
# 15 bins means 16 edges, logarithmically spaced
diam_edges_nm = diam_min_nm * (diam_max_nm / diam_min_nm) ** (np.arange(16) / 15)

# Convert to μm for plotting
diam_edges_um = diam_edges_nm / 1000.0  # nm → μm
diam_centers_um = np.sqrt(diam_edges_um[:-1] * diam_edges_um[1:])  # Geometric mean

# Also calculate radius for reference (radius = diameter/2)
radius_edges_um = diam_edges_um / 2.0
radius_centers_um = diam_centers_um / 2.0

print("🔬 TOMAS-15 Size Bins (Dry Particle Diameter):")
print("   Bin  | Diameter Range (μm)   | Center (μm)")
print("   " + "-" * 50)
for i in range(15):
    print(f"   {i+1:2d}   | {diam_edges_um[i]:8.4f} - {diam_edges_um[i+1]:8.4f} | {diam_centers_um[i]:8.4f}")
print()

# Aerosol species (mass concentration species only)
# NOTE: NK = Number concentration (#/cm³), NOT nitrate mass!
#       Actual nitrate mass is in SpeciesConcVV_NIT (no bins)
aerosol_types = [
    ("DUST", "Mineral Dust"),
    ("SS", "Sea Salt"),
    ("SF", "Sulfate"),
    ("ECIL", "BC Hydrophilic"),
    ("ECOB", "BC Hydrophobic"),
    ("OCIL", "OC Hydrophilic"),
    ("OCOB", "OC Hydrophobic"),
    # ("NK", "Number Conc."),  # Excluded - different units (#/cm³ not mol/mol)
    ("AW", "Aerosol Water")
]

# Get atmospheric structure
n_lev = ds.dimensions['lev'].size
dp = ds['Met_DELP'][0, :, idf, idy, idx]
sp = ds['Met_PS2WET'][0, idf, idy, idx]

# Pressure at half-levels (flip from BOA→TOA to TOA→BOA)
pressure_half = np.concatenate([[sp], sp - np.cumsum(dp)])
pressure_half_toa2boa = pressure_half[::-1]
pressure_centers = (pressure_half_toa2boa[:-1] + pressure_half_toa2boa[1:]) / 2

# Temperature (flip to TOA→BOA)
temperature = ds['Met_T'][0, :, idf, idy, idx][::-1]

print("🌍 Atmospheric Profile:")
print(f"   Levels: {n_lev}")
print(f"   Surface Pressure: {sp:.1f} hPa")
print(f"   Top Pressure: {pressure_half.min():.3f} hPa")
print(f"   Temperature range: {temperature.min():.1f} - {temperature.max():.1f} K")
print()

# Read all species data
print("📊 Reading Aerosol Species Profiles...")
print()

species_data = {}

for prefix, name in aerosol_types:
    # Preallocate: 15 bins × n_lev layers
    concentration = np.zeros((15, n_lev))
    
    for bin_idx in range(1, 16):
        var_name = f"SpeciesConcVV_{prefix}{bin_idx:02d}"
        if var_name in ds.variables:
            # VMR [mol/mol dry], flip BOA→TOA to TOA→BOA
            concentration[bin_idx-1, :] = ds[var_name][0, :, idf, idy, idx][::-1]
    
    species_data[prefix] = concentration
    
    # Statistics
    total_conc = concentration.sum()
    max_conc = concentration.max()
    
    if max_conc > 1e-30:
        print(f"   ✓ {name} ({prefix}):")
        print(f"      Total concentration: {total_conc:.2e} mol/mol")
        print(f"      Max concentration:   {max_conc:.2e} mol/mol")
        
        # Find peak
        bin_max, lev_max = np.unravel_index(concentration.argmax(), concentration.shape)
        print(f"      Peak at: Bin {bin_max+1} (d={diam_centers_um[bin_max]:.3f} μm), "
              f"Level {lev_max+1} (p={pressure_centers[lev_max]:.1f} hPa)")
    else:
        print(f"   ○ {name} ({prefix}): negligible")

# Don't close dataset yet - we need it for multi-location comparison
# ds.close()

print()
print("=" * 70)
print("Creating Visualizations...")
print("=" * 70)
print()

# Plot 1: Vertical profiles
print("📈 Plot 1: Vertical Profiles (Total Concentration)")
fig, ax = plt.subplots(figsize=(10, 7))

for prefix, name in aerosol_types:
    conc = species_data[prefix]
    total_per_layer = conc.sum(axis=0)  # Sum over all bins
    
    if total_per_layer.max() > 1e-30:
        ax.plot(total_per_layer, pressure_centers, label=name, linewidth=2)

ax.set_xlabel('Total Concentration [mol/mol]', fontsize=12)
ax.set_ylabel('Pressure [hPa]', fontsize=12)
ax.set_xscale('log')
ax.set_yscale('log')
ax.invert_yaxis()
ax.set_title('Aerosol Vertical Profiles', fontsize=14, fontweight='bold')
ax.legend(loc='lower right', fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_dir / "01_vertical_profiles.png", dpi=150)
print("   Saved: 01_vertical_profiles.png")
plt.close()

# Plot 2: Size distributions at 800 and 500 hPa with log-normal fits
print("📈 Plot 2: Size Distributions at 800 and 500 hPa (with log-normal fits)")

# Select pressure levels
lev_800 = np.argmin(np.abs(pressure_centers - 800))
lev_500 = np.argmin(np.abs(pressure_centers - 500))
print(f"   800 hPa: level {lev_800}, pressure = {pressure_centers[lev_800]:.1f} hPa")
print(f"   500 hPa: level {lev_500}, pressure = {pressure_centers[lev_500]:.1f} hPa")

# Mass species only (NK excluded - it's number concentration)
major_species = ["DUST", "SS", "SF", "ECIL", "OCIL", "OCOB", "ECOB", "AW"]

# Create subplot with two panels
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))

# Fine diameter grid for smooth log-normal curves
d_fine = np.logspace(np.log10(diam_centers_um.min()), 
                      np.log10(diam_centers_um.max()), 200)

# Calculate dlogD for TOMAS bins (bin width in log space)
dlogD = np.log10(diam_edges_um[1:] / diam_edges_um[:-1])

# Plot at 800 hPa
print(f"\n   Fitting log-normals at {pressure_centers[lev_800]:.0f} hPa:")
for species_prefix in major_species:
    conc = species_data[species_prefix]
    
    if conc.max() < 1e-30:
        continue
    
    conc_at_level = conc[:, lev_800]
    species_name = dict(aerosol_types)[species_prefix]
    
    # Convert to dN/dlogD (concentration per log-diameter bin)
    dN_dlogD = conc_at_level / dlogD
    
    # Plot data
    ax1.plot(diam_centers_um, dN_dlogD, 
            marker='o', linewidth=2, markersize=6,
            label=species_name, alpha=0.7)
    
    # Fit log-normal to dN/dlogD
    params = fit_lognormal(diam_centers_um, dN_dlogD)
    if params is not None:
        N_total, d_med, sigma_g = params
        conc_fit = lognormal(d_fine, N_total, d_med, sigma_g)
        ax1.plot(d_fine, conc_fit, '--', linewidth=1.5, alpha=0.6)
        print(f"     {species_prefix:6s}: d_med={d_med:.3f} μm, σ_g={sigma_g:.2f}")

ax1.set_xlabel('Dry Diameter [μm]', fontsize=12)
ax1.set_ylabel('dN/dlogD [mol/mol]', fontsize=12)
ax1.set_xscale('log')
# ax1.set_yscale('log')  # Linear y-axis
ax1.set_title(f'Size Distributions at {pressure_centers[lev_800]:.0f} hPa (Lower Troposphere)', 
              fontsize=13, fontweight='bold')
ax1.legend(fontsize=8, ncol=2, loc='best')
ax1.grid(True, alpha=0.3)

# Plot at 500 hPa
print(f"\n   Fitting log-normals at {pressure_centers[lev_500]:.0f} hPa:")
for species_prefix in major_species:
    conc = species_data[species_prefix]
    
    if conc.max() < 1e-30:
        continue
    
    conc_at_level = conc[:, lev_500]
    species_name = dict(aerosol_types)[species_prefix]
    
    # Convert to dN/dlogD
    dN_dlogD = conc_at_level / dlogD
    
    # Plot data
    ax2.plot(diam_centers_um, dN_dlogD, 
            marker='o', linewidth=2, markersize=6,
            label=species_name, alpha=0.7)
    
    # Fit log-normal to dN/dlogD
    params = fit_lognormal(diam_centers_um, dN_dlogD)
    if params is not None:
        N_total, d_med, sigma_g = params
        conc_fit = lognormal(d_fine, N_total, d_med, sigma_g)
        ax2.plot(d_fine, conc_fit, '--', linewidth=1.5, alpha=0.6)
        print(f"     {species_prefix:6s}: d_med={d_med:.3f} μm, σ_g={sigma_g:.2f}")

ax2.set_xlabel('Dry Diameter [μm]', fontsize=12)
ax2.set_ylabel('dN/dlogD [mol/mol]', fontsize=12)
ax2.set_xscale('log')
# ax2.set_yscale('log')  # Linear y-axis
ax2.set_title(f'Size Distributions at {pressure_centers[lev_500]:.0f} hPa (Mid-Troposphere)', 
              fontsize=13, fontweight='bold')
ax2.legend(fontsize=8, ncol=2, loc='best')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / "02_size_distributions_800_500hPa.png", dpi=150)
print(f"\n   Saved: 02_size_distributions_800_500hPa.png")
plt.close()

# Plot 3: 2D heatmaps
print("📈 Plot 3: Altitude-Size Heatmaps")

for species_prefix in major_species:
    conc = species_data[species_prefix]
    
    if conc.max() < 1e-30:
        continue
    
    # Log scale for visualization
    conc_log = np.log10(conc + 1e-35)
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    im = ax.pcolormesh(np.arange(1, 16), pressure_centers, conc_log.T,
                       shading='auto', cmap='viridis')
    
    ax.set_xlabel('Bin Number', fontsize=12)
    ax.set_ylabel('Pressure [hPa]', fontsize=12)
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_title(f'{species_prefix}: log₁₀(Concentration)', fontsize=14, fontweight='bold')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('log₁₀[mol/mol]', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_dir / f"03_heatmap_{species_prefix}.png", dpi=150)
    print(f"   Saved: 03_heatmap_{species_prefix}.png")
    plt.close()

# Plot 4: Species comparison at boundary layer
print("📈 Plot 4: Species Comparison at Boundary Layer")

bl_level = np.argmax(pressure_centers < 900)

fig, ax = plt.subplots(figsize=(10, 7))

for prefix, name in aerosol_types:
    conc = species_data[prefix]
    conc_at_bl = conc[:, bl_level]
    
    if conc_at_bl.max() > 1e-30:
        ax.plot(diam_centers_um, conc_at_bl, 
                marker='o', linewidth=2, markersize=6,
                label=name)

ax.set_xlabel('Dry Diameter [μm]', fontsize=12)
ax.set_ylabel('Concentration [mol/mol]', fontsize=12)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title(f'Size Distributions at {pressure_centers[bl_level]:.0f} hPa', 
             fontsize=14, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_dir / "04_all_species_comparison.png", dpi=150)
print("   Saved: 04_all_species_comparison.png")
plt.close()

# Plot 5: Multi-location comparison
print("📈 Plot 5: Size Distributions Across Multiple Locations")

# Select dominant species for comparison
species_to_compare = ["OCIL", "SF", "DUST", "SS"]

# Use 800 hPa level
lev_800 = np.argmin(np.abs(pressure_centers - 800))

# Create subplots: one per species
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
axes = axes.flatten()

for i, species_prefix in enumerate(species_to_compare):
    ax = axes[i]
    species_name = dict(aerosol_types)[species_prefix]
    
    # Loop over all locations
    for loc_idx, (idx, idy, idf, loc_name) in enumerate(locations):
        # Read concentration data for this location
        conc_loc = np.zeros((15, n_lev))
        for bin_idx in range(1, 16):
            var_name = f"SpeciesConcVV_{species_prefix}{bin_idx:02d}"
            if var_name in ds.variables:
                data = ds[var_name][0, :, idf, idy, idx][::-1]
                conc_loc[bin_idx-1, :] = data
        
        # Get concentration at 800 hPa
        conc_at_level = conc_loc[:, lev_800]
        
        if conc_at_level.max() > 1e-30:
            # Convert to dN/dlogD
            dN_dlogD = conc_at_level / dlogD
            
            # Plot
            ax.plot(diam_centers_um, dN_dlogD, 
                   marker='o', linewidth=2, markersize=5,
                   label=loc_name, alpha=0.8)
    
    ax.set_xlabel('Dry Diameter [μm]', fontsize=11)
    ax.set_ylabel('dN/dlogD [mol/mol]', fontsize=11)
    ax.set_xscale('log')
    ax.set_title(f'{species_name}', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, loc='best')
    ax.grid(True, alpha=0.3)

plt.suptitle(f'Aerosol Size Distributions at {pressure_centers[lev_800]:.0f} hPa - Multi-Location Comparison', 
             fontsize=15, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig(output_dir / "05_multi_location_comparison.png", dpi=150)
print("   Saved: 05_multi_location_comparison.png")
plt.close()

# Close dataset
ds.close()

print()
print("=" * 70)
print("✅ Exploration Complete!")
print("=" * 70)
print()
print(f"Outputs saved to: {output_dir}/")
print()
print("Key Findings:")
print(f"  • Dry diameter range: {diam_edges_um[0]:.3f} - {diam_edges_um[-1]:.1f} μm")
print(f"  • Number of bins: 15")
print(f"  • Vertical levels: {n_lev}")
active_species = sum(1 for prefix, _ in aerosol_types if species_data[prefix].max() > 1e-30)
print(f"  • Active species: {active_species}")
print()
print("Next Steps:")
print("  1. Review the generated plots")
print("  2. Identify which species/sizes dominate")
print("  3. Decide on fitting strategy (if needed)")
print("  4. Implement optical property calculations")
print()
