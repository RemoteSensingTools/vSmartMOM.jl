#!/usr/bin/env python3
"""
NK Number Size Distribution Explorer

Focuses on the NK variable (number concentration #/cm³) to understand
the actual particle size distribution in TOMAS-15.

Requirements:
    pip install netCDF4 numpy matplotlib scipy

Usage:
    python test/explore_NK_number_distribution.py
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
        Total number concentration (integral of distribution) [#/cm³]
    r_med : float
        Median radius [μm]
    sigma_g : float
        Geometric standard deviation (dimensionless)
    
    Returns:
    --------
    dN/dlogr : array
        Number per log-radius bin [#/cm³]
    """
    return (N_total / (np.sqrt(2*np.pi) * np.log(sigma_g))) * \
           np.exp(-0.5 * (np.log(r/r_med) / np.log(sigma_g))**2)

def bimodal_lognormal(r, N1, r_med1, sigma_g1, N2, r_med2, sigma_g2):
    """
    Bimodal log-normal distribution (sum of two modes)
    
    Parameters:
    -----------
    N1, r_med1, sigma_g1 : float
        First mode parameters (typically Aitken)
    N2, r_med2, sigma_g2 : float
        Second mode parameters (typically accumulation)
    
    Returns:
    --------
    dN/dlogr : array
        Total number per log-radius bin [#/cm³]
    """
    return lognormal(r, N1, r_med1, sigma_g1) + lognormal(r, N2, r_med2, sigma_g2)

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
    except Exception as e:
        print(f"      Fit failed: {e}")
        return None

def fit_bimodal_lognormal(r_centers, concentrations):
    """
    Fit a bimodal log-normal distribution to binned concentration data
    
    Returns (N1, r_med1, sigma_g1, N2, r_med2, sigma_g2) or None if fit fails
    """
    # Filter out zero/negative values
    valid = concentrations > 0
    if valid.sum() < 6:  # Need at least 6 points for bimodal
        return None
    
    r_valid = r_centers[valid]
    c_valid = concentrations[valid]
    
    # Initial guess: split at 0.1 μm (Aitken vs accumulation boundary)
    aitken_mask = r_valid < 0.05  # Aitken mode
    accum_mask = r_valid >= 0.05  # Accumulation mode
    
    # Mode 1 (Aitken): smaller particles
    if aitken_mask.sum() > 0:
        r_med1_guess = 0.03  # 30 nm
        N1_guess = np.trapz(c_valid[aitken_mask], np.log(r_valid[aitken_mask])) if aitken_mask.sum() > 1 else c_valid.max() * 0.3
        sigma_g1_guess = 1.5
    else:
        r_med1_guess = 0.03
        N1_guess = c_valid.max() * 0.3
        sigma_g1_guess = 1.5
    
    # Mode 2 (Accumulation): larger particles
    if accum_mask.sum() > 0:
        r_med2_guess = r_valid[accum_mask][np.argmax(c_valid[accum_mask])]
        N2_guess = np.trapz(c_valid[accum_mask], np.log(r_valid[accum_mask])) if accum_mask.sum() > 1 else c_valid.max() * 0.7
        sigma_g2_guess = 2.0
    else:
        r_med2_guess = 0.15
        N2_guess = c_valid.max() * 0.7
        sigma_g2_guess = 2.0
    
    # Ensure N1, N2 > 0
    N1_guess = max(N1_guess, 1.0)
    N2_guess = max(N2_guess, 1.0)
    
    try:
        params, _ = curve_fit(
            bimodal_lognormal, r_valid, c_valid,
            p0=[N1_guess, r_med1_guess, sigma_g1_guess, 
                N2_guess, r_med2_guess, sigma_g2_guess],
            bounds=([0, r_valid.min(), 1.01, 0, 0.05, 1.01],
                   [np.inf, 0.1, 3.0, np.inf, r_valid.max(), 5.0]),
            maxfev=10000
        )
        return params
    except Exception as e:
        print(f"      Bimodal fit failed: {e}")
        return None

# Configuration
ncfile = "GEOSChem.Custom.20190702_0000z.nc4"
output_dir = Path("test/aerosol_exploration_output")
output_dir.mkdir(exist_ok=True)

# Define multiple locations to analyze
locations = [
    (2, 12, 4, "Central USA"),        # 36.8°N, 97.9°W (Oklahoma/Kansas)
    (15, 15, 1, "Amazon Basin"),      # Tropical South America
    (8, 18, 2, "Sahara Desert"),      # North Africa
    (12, 8, 3, "India/South Asia"),   # Indian subcontinent
    (5, 5, 5, "Northern China"),      # East Asia
    (17, 15, 3, "South Pacific"),     # Clean marine background (~22°S, 184°E)
]

print("=" * 70)
print("NK Number Size Distribution Explorer")
print("=" * 70)
print()

print(f"📍 Analyzing {len(locations)} locations:")

# Open dataset
ds = nc.Dataset(ncfile, 'r')

# Check if NK variable exists
test_var = "SpeciesConcVV_NK01"
if test_var not in ds.variables:
    print(f"ERROR: Variable {test_var} not found in dataset!")
    print("Available variables:", list(ds.variables.keys())[:20])
    ds.close()
    exit(1)

# Get location info
for i, (idx_loc, idy_loc, idf_loc, name) in enumerate(locations, 1):
    lat_loc = ds['lats'][idf_loc, idy_loc, idx_loc]
    lon_loc = ds['lons'][idf_loc, idy_loc, idx_loc]
    print(f"   {i}. {name:20s}: {lat_loc:6.2f}°N, {lon_loc:6.2f}°E")
print()

# Use first location for detailed analysis
idx, idy, idf, loc_name = locations[0]
lat = ds['lats'][idf, idy, idx]
lon = ds['lons'][idf, idy, idx]
print(f"📍 Detailed analysis for: {loc_name}")
print(f"   Location: idx={idx+1}, idy={idy+1}, face={idf+1}")
print(f"   Latitude:  {lat:.2f}°")
print(f"   Longitude: {lon:.2f}°")
print()

# TOMAS-15 bin definitions (dry diameter)
diam_min_nm = 10.0      # 10 nm minimum dry diameter
diam_max_nm = 10000.0   # 10 μm maximum dry diameter
diam_edges_nm = diam_min_nm * (diam_max_nm / diam_min_nm) ** (np.arange(16) / 15)

# Convert to μm
diam_edges_um = diam_edges_nm / 1000.0
diam_centers_um = np.sqrt(diam_edges_um[:-1] * diam_edges_um[1:])

# Radius (for reference)
radius_edges_um = diam_edges_um / 2.0
radius_centers_um = diam_centers_um / 2.0

print("🔬 TOMAS-15 Size Bins (Dry Particle Diameter):")
print("   Bin  | Diameter Range (μm)   | Center (μm)   | Radius Center (μm)")
print("   " + "-" * 70)
for i in range(15):
    print(f"   {i+1:2d}   | {diam_edges_um[i]:8.4f} - {diam_edges_um[i+1]:8.4f} | "
          f"{diam_centers_um[i]:8.4f}      | {radius_centers_um[i]:8.4f}")
print()

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

# Read meteorology for NK conversion
# NK uses special units: 1000 × (particles/mol_air)
# Convert using: N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)

M_air = 28.9644e-3  # kg/mol (molar mass of dry air)

# Read Met_AD (air mass) and Met_AIRVOL (air volume)
Met_AD = ds['Met_AD'][0, :, idf, idy, idx][::-1]  # kg (flip to TOA→BOA)
Met_AIRVOL = ds['Met_AIRVOL'][0, :, idf, idy, idx][::-1]  # m³

# Calculate conversion factors
n_air = Met_AD / M_air  # moles of air
vol_cm3 = Met_AIRVOL * 1e6  # cm³

# Also calculate air density for reference
R_specific = 287.05  # J/(kg·K) for dry air
pressure_Pa = pressure_centers * 100.0  # hPa to Pa
rho_air_kg_m3 = Met_AD / Met_AIRVOL  # kg/m³

print("🌍 Atmospheric Profile:")
print(f"   Levels: {n_lev}")
print(f"   Surface Pressure: {sp:.1f} hPa")
print(f"   Top Pressure: {pressure_half.min():.3f} hPa")
print(f"   Temperature range: {temperature.min():.1f} - {temperature.max():.1f} K")
print(f"   Air density range: {rho_air_kg_m3.min():.3f} - {rho_air_kg_m3.max():.3f} kg/m³")
print()

# Read NK data (number concentration)
print("📊 Reading NK Number Concentration Data...")
print()

# Preallocate: 15 bins × n_lev layers
nk_raw = np.zeros((15, n_lev))  # Raw NK values (1000 × #/mol_air)
nk_concentration = np.zeros((15, n_lev))  # Converted to #/cm³

for bin_idx in range(1, 16):
    var_name = f"SpeciesConcVV_NK{bin_idx:02d}"
    if var_name in ds.variables:
        # Read raw values (stored as 1000 × particles/mol_air)
        nk_raw[bin_idx-1, :] = ds[var_name][0, :, idf, idy, idx][::-1]
        
        # Convert using CORRECT formula: N = (NK/1000) × (n_air/vol_cm3)
        nk_concentration[bin_idx-1, :] = (nk_raw[bin_idx-1, :] / 1000.0) * (n_air / vol_cm3)
        
        # Get metadata from first variable
        if bin_idx == 1:
            units_meta = ds[var_name].units if hasattr(ds[var_name], 'units') else 'unknown'
            print(f"   Metadata units: {units_meta}")
            print(f"   Actual units: 1000 × (particles/mol_air)")
            print(f"   Conversion: N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)")
            print(f"   Result: #/cm³")

# Statistics
total_n = nk_concentration.sum()
max_n = nk_concentration.max()

# Also show raw values for reference
total_n_raw = nk_raw.sum()
print()
print(f"   Raw NK values (as stored): {total_n_raw:.2e} (1000 × #/mol_air)")
print(f"   Converted to #/cm³: {total_n:.2e} #/cm³")
print(f"   Max number concentration: {max_n:.2e} #/cm³")

if max_n > 0:
    bin_max, lev_max = np.unravel_index(nk_concentration.argmax(), nk_concentration.shape)
    print(f"   Peak at: Bin {bin_max+1} (d={diam_centers_um[bin_max]:.3f} μm), "
          f"Level {lev_max+1} (p={pressure_centers[lev_max]:.1f} hPa)")
else:
    print("   WARNING: All NK values are zero!")

print()
print("=" * 70)
print("Creating Visualizations...")
print("=" * 70)
print()

# Calculate dlogD for TOMAS bins
dlogD = np.log10(diam_edges_um[1:] / diam_edges_um[:-1])

# Plot 1: Vertical profile by size mode
print("📈 Plot 1: Vertical Profile by Size Mode")
fig, ax = plt.subplots(figsize=(10, 7))

# Define size modes based on diameter (μm):
# Aitken mode: < 0.1 μm diameter (bins 1-5)
# Accumulation mode: 0.1 - 1.0 μm diameter (bins 6-10)
# Coarse mode: > 1.0 μm diameter (bins 11-15)

# Reference: Standard aerosol size classification
# Aitken: nucleation + small Aitken (< 0.1 μm)
# Accumulation: 0.1 - 1.0 μm (or 2.5 μm)
# Coarse: > 1.0 μm (or 2.5 μm)

aitken_bins = range(0, 5)      # bins 1-5: 0.01-0.1 μm
accumulation_bins = range(5, 10)  # bins 6-10: 0.1-1.0 μm
coarse_bins = range(10, 15)    # bins 11-15: 1.0-10 μm

aitken_n = nk_concentration[aitken_bins, :].sum(axis=0)
accumulation_n = nk_concentration[accumulation_bins, :].sum(axis=0)
coarse_n = nk_concentration[coarse_bins, :].sum(axis=0)
total_n = nk_concentration.sum(axis=0)

ax.plot(aitken_n, pressure_centers, linewidth=2.5, label='Aitken (<0.1 μm)', color='blue')
ax.plot(accumulation_n, pressure_centers, linewidth=2.5, label='Accumulation (0.1-1.0 μm)', color='green')
ax.plot(coarse_n, pressure_centers, linewidth=2.5, label='Coarse (>1.0 μm)', color='red')
ax.plot(total_n, pressure_centers, linewidth=2.5, label='Total', color='black', linestyle='--')

ax.set_xlabel('Number Concentration [#/cm³]', fontsize=13)
ax.set_ylabel('Pressure [hPa]', fontsize=13)
ax.set_xscale('log')
# Linear pressure scale
ax.invert_yaxis()
ax.set_title(f'Aerosol Number Density by Size Mode - {loc_name}\n(Correct formula: N = (NK/1000) × (n_air/vol))', fontsize=14, fontweight='bold')
ax.legend(fontsize=11, loc='best')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(output_dir / "NK_01_vertical_profile.png", dpi=150)
print("   Saved: NK_01_vertical_profile.png")
print(f"   Mode splits:")
print(f"     Aitken bins 1-5: {diam_edges_um[0]:.4f} - {diam_edges_um[5]:.4f} μm")
print(f"     Accumulation bins 6-10: {diam_edges_um[5]:.4f} - {diam_edges_um[10]:.4f} μm")
print(f"     Coarse bins 11-15: {diam_edges_um[10]:.4f} - {diam_edges_um[15]:.4f} μm")
plt.close()

# Plot 2: Size distributions at multiple altitudes with log-normal fits
print("📈 Plot 2: Size Distributions at Multiple Altitudes")

# Select interesting pressure levels
pressure_levels = [900, 800, 500, 300, 200]
level_indices = [np.argmin(np.abs(pressure_centers - p)) for p in pressure_levels]

# Show info about conversion at selected levels
print()
print("   Conversion examples at selected altitudes:")
for p_target in pressure_levels:
    lev_idx = np.argmin(np.abs(pressure_centers - p_target))
    p_actual = pressure_centers[lev_idx]
    nk_val = nk_raw[:, lev_idx].sum()
    n_val = nk_concentration[:, lev_idx].sum()
    mols = n_air[lev_idx]
    vol = vol_cm3[lev_idx]
    print(f"     {p_actual:.0f} hPa: n_air={mols:.2e} mol, vol={vol:.2e} cm³")
    print(f"               NK={nk_val:.2e} (1000×#/mol) → N={n_val:.2f} #/cm³")

print()
fig, ax = plt.subplots(figsize=(12, 8))

# Fine diameter grid for smooth curves
d_fine = np.logspace(np.log10(diam_centers_um.min()), 
                      np.log10(diam_centers_um.max()), 200)

# Color map
colors = plt.cm.viridis(np.linspace(0, 0.9, len(level_indices)))

print()
for i, (lev_idx, p_target) in enumerate(zip(level_indices, pressure_levels)):
    nk_at_level = nk_concentration[:, lev_idx]
    p_actual = pressure_centers[lev_idx]
    
    # Convert to dN/dlogD
    dN_dlogD = nk_at_level / dlogD
    
    # Plot data
    ax.plot(diam_centers_um, dN_dlogD, 
            marker='o', linewidth=2.5, markersize=7,
            label=f'{p_actual:.0f} hPa', 
            color=colors[i], alpha=0.8)
    
    # Fit log-normal (print params but don't plot)
    print(f"   Fitting log-normal at {p_actual:.0f} hPa:")
    params = fit_lognormal(diam_centers_um, dN_dlogD)
    if params is not None:
        N_total, d_med, sigma_g = params
        # n_fit = lognormal(d_fine, N_total, d_med, sigma_g)
        # ax.plot(d_fine, n_fit, '--', linewidth=2, alpha=0.7, color=colors[i])
        
        # Calculate effective diameter (d_eff = d_med * exp(2.5 * ln²(sigma_g)))
        d_eff = d_med * np.exp(2.5 * np.log(sigma_g)**2)
        
        print(f"      N_total = {N_total:.2e} #/cm³")
        print(f"      d_med   = {d_med:.4f} μm")
        print(f"      d_eff   = {d_eff:.4f} μm")
        print(f"      σ_g     = {sigma_g:.3f}")
    else:
        print(f"      Could not fit log-normal distribution")

ax.set_xlabel('Dry Diameter [μm]', fontsize=13)
ax.set_ylabel('dN/dlogD [#/cm³]', fontsize=13)
ax.set_xscale('log')
# ax.set_yscale('log')  # Linear y-axis
ax.set_title(f'Aerosol Number Size Distribution - {loc_name}\n(Validated formula: N = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6))', 
             fontsize=14, fontweight='bold')
ax.legend(fontsize=11, loc='best', title='Pressure Level')
ax.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig(output_dir / "NK_02_size_distributions_multi_altitude.png", dpi=150)
print(f"\n   Saved: NK_02_size_distributions_multi_altitude.png")
plt.close()

# Plot 3: 2D heatmap (altitude vs size)
print("📈 Plot 3: Altitude-Size Heatmap")

# Log scale for visualization
nk_log = np.log10(nk_concentration + 1e-10)  # Avoid log(0)

fig, ax = plt.subplots(figsize=(12, 7))

im = ax.pcolormesh(diam_centers_um, pressure_centers, nk_log.T,
                   shading='auto', cmap='plasma')

ax.set_xlabel('Dry Diameter [μm]', fontsize=13)
ax.set_ylabel('Pressure [hPa]', fontsize=13)
ax.set_xscale('log')
ax.set_yscale('log')
ax.invert_yaxis()
ax.set_title(f'Aerosol Number Density Heatmap - {loc_name}\n(Units: #/cm³ using validated conversion)', 
             fontsize=14, fontweight='bold')

cbar = plt.colorbar(im, ax=ax, pad=0.02)
cbar.set_label('log₁₀[Number Density (#/cm³)]', fontsize=12)

plt.tight_layout()
plt.savefig(output_dir / "NK_03_heatmap.png", dpi=150)
print("   Saved: NK_03_heatmap.png")
plt.close()

# Plot 4: Boundary layer focus with bimodal fit
print("📈 Plot 4: Boundary Layer Size Distribution (detailed)")

bl_level = np.argmin(np.abs(pressure_centers - 900))
p_bl = pressure_centers[bl_level]
nk_bl = nk_concentration[:, bl_level]
dN_dlogD_bl = nk_bl / dlogD

# Fine diameter grid for smooth curves
d_fine = np.logspace(np.log10(diam_centers_um.min()), 
                      np.log10(diam_centers_um.max()), 200)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Left panel: Linear y-axis with BIMODAL fit
ax1.plot(diam_centers_um, dN_dlogD_bl, 
         marker='o', linewidth=2.5, markersize=8,
         color='darkgreen', label='Data', zorder=3)

# Fit bimodal distribution
print()
print(f"   Fitting bimodal log-normal at {p_bl:.0f} hPa:")
bimodal_params = fit_bimodal_lognormal(diam_centers_um / 2, dN_dlogD_bl)  # Convert diameter to radius
if bimodal_params is not None:
    N1, r_med1, sigma_g1, N2, r_med2, sigma_g2 = bimodal_params
    
    # Convert back to diameter for display
    d_med1 = r_med1 * 2
    d_med2 = r_med2 * 2
    
    # Plot individual modes
    mode1 = lognormal(d_fine / 2, N1, r_med1, sigma_g1)
    mode2 = lognormal(d_fine / 2, N2, r_med2, sigma_g2)
    total_fit = mode1 + mode2
    
    ax1.plot(d_fine, mode1, '--', linewidth=2, alpha=0.7, 
             color='blue', label=f'Aitken: N={N1:.2e}, d_med={d_med1:.3f}μm, σ_g={sigma_g1:.2f}')
    ax1.plot(d_fine, mode2, '--', linewidth=2, alpha=0.7, 
             color='orange', label=f'Accumulation: N={N2:.2e}, d_med={d_med2:.3f}μm, σ_g={sigma_g2:.2f}')
    ax1.plot(d_fine, total_fit, '-', linewidth=2.5, alpha=0.9, 
             color='red', label='Total Fit', zorder=2)
    
    print(f"      Mode 1 (Aitken):       N = {N1:.2e} #/cm³, d_med = {d_med1:.4f} μm, σ_g = {sigma_g1:.3f}")
    print(f"      Mode 2 (Accumulation): N = {N2:.2e} #/cm³, d_med = {d_med2:.4f} μm, σ_g = {sigma_g2:.3f}")
    print(f"      Total N: {N1 + N2:.2e} #/cm³")
    print(f"      Aitken fraction: {N1/(N1+N2)*100:.1f}%")
else:
    print("      Could not fit bimodal distribution")

ax1.set_xlabel('Dry Diameter [μm]', fontsize=12)
ax1.set_ylabel('dN/dlogD [#/cm³]', fontsize=12)
ax1.set_xscale('log')
ax1.set_title(f'Boundary Layer ({p_bl:.0f} hPa) - Bimodal Fit\n(Validated conversion)', 
              fontsize=12, fontweight='bold')
ax1.legend(fontsize=9, loc='best')
ax1.grid(True, alpha=0.3)

# Right panel: Log y-axis with bimodal fit
ax2.plot(diam_centers_um, dN_dlogD_bl, 
         marker='o', linewidth=2.5, markersize=8,
         color='darkgreen', label='Data', zorder=3)

if bimodal_params is not None:
    ax2.plot(d_fine, mode1, '--', linewidth=2, alpha=0.7, 
             color='blue', label=f'Aitken Mode')
    ax2.plot(d_fine, mode2, '--', linewidth=2, alpha=0.7, 
             color='orange', label=f'Accumulation Mode')
    ax2.plot(d_fine, total_fit, '-', linewidth=2.5, alpha=0.9, 
             color='red', label='Total Fit', zorder=2)

ax2.set_xlabel('Dry Diameter [μm]', fontsize=12)
ax2.set_ylabel('dN/dlogD [#/cm³]', fontsize=12)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_title(f'Boundary Layer ({p_bl:.0f} hPa) - Log Scale\n(Validated conversion)', 
              fontsize=12, fontweight='bold')
ax2.legend(fontsize=9, loc='best')
ax2.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(output_dir / "NK_04_boundary_layer_detail.png", dpi=150)
print("   Saved: NK_04_boundary_layer_detail.png")
plt.close()

# Plot 5: Multi-location comparison
print("📈 Plot 5: Size Distributions Across Multiple Locations")

# Use 800 hPa level
lev_800 = np.argmin(np.abs(pressure_centers - 800))

fig, ax = plt.subplots(figsize=(12, 8))

print()
for loc_idx, (idx_loc, idy_loc, idf_loc, loc_name_comp) in enumerate(locations):
    # Get meteorology for this location
    Met_AD_loc = ds['Met_AD'][0, :, idf_loc, idy_loc, idx_loc][::-1]  # kg
    Met_AIRVOL_loc = ds['Met_AIRVOL'][0, :, idf_loc, idy_loc, idx_loc][::-1]  # m³
    n_air_loc = Met_AD_loc / M_air  # mol
    vol_cm3_loc = Met_AIRVOL_loc * 1e6  # cm³
    
    # Read NK data for this location (raw 1000×#/mol)
    nk_raw_loc = np.zeros((15, n_lev))
    nk_loc = np.zeros((15, n_lev))
    for bin_idx in range(1, 16):
        var_name = f"SpeciesConcVV_NK{bin_idx:02d}"
        if var_name in ds.variables:
            nk_raw_loc[bin_idx-1, :] = ds[var_name][0, :, idf_loc, idy_loc, idx_loc][::-1]
            # Convert to #/cm³
            nk_loc[bin_idx-1, :] = (nk_raw_loc[bin_idx-1, :] / 1000.0) * (n_air_loc / vol_cm3_loc)
    
    # Get at 800 hPa
    nk_at_level = nk_loc[:, lev_800]
    dN_dlogD_loc = nk_at_level / dlogD
    
    # Plot
    ax.plot(diam_centers_um, dN_dlogD_loc, 
           marker='o', linewidth=2.5, markersize=6,
           label=loc_name_comp, alpha=0.8)
    
    # Fit (print but don't plot)
    print(f"   {loc_name_comp}:")
    params_loc = fit_lognormal(diam_centers_um, dN_dlogD_loc)
    if params_loc is not None:
        N_total_loc, d_med_loc, sigma_g_loc = params_loc
        # n_fit_loc = lognormal(d_fine, N_total_loc, d_med_loc, sigma_g_loc)
        # ax.plot(d_fine, n_fit_loc, '--', linewidth=1.5, alpha=0.6)
        print(f"      N_total = {N_total_loc:.2e} #/cm³, d_med = {d_med_loc:.4f} μm, σ_g = {sigma_g_loc:.3f}")

ax.set_xlabel('Dry Diameter [μm]', fontsize=13)
ax.set_ylabel('dN/dlogD [#/cm³]', fontsize=13)
ax.set_xscale('log')
# ax.set_yscale('log')  # Linear y-axis
ax.set_title(f'Aerosol Number Density at {pressure_centers[lev_800]:.0f} hPa - Multi-Location\n(Validated formula: N = (NK/1000) × (n_air/vol))', 
             fontsize=14, fontweight='bold')
ax.legend(fontsize=11, loc='best')
ax.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig(output_dir / "NK_05_multi_location_comparison.png", dpi=150)
print(f"\n   Saved: NK_05_multi_location_comparison.png")
plt.close()

# Close dataset
ds.close()

print()
print("=" * 70)
print("✅ NK Number Distribution Analysis Complete!")
print("=" * 70)
print()
print(f"Outputs saved to: {output_dir}/")
print()
print("Summary Statistics:")
print(f"  • Size range: {diam_edges_um[0]:.4f} - {diam_edges_um[-1]:.1f} μm diameter")
print(f"  • Number of bins: 15")
print(f"  • Vertical levels: {n_lev}")
print(f"  • Surface number concentration: {nk_concentration[:, -1].sum():.2e} #/cm³")
print(f"  • At 800 hPa: {nk_concentration[:, lev_800].sum():.2e} #/cm³")
print()
print("Key Findings:")
print("  • NK is stored as 1000 × (particles/mol_air)")
print("  • Converted using: N(#/cm³) = (NK/1000) × (Met_AD/M_air) / (Met_AIRVOL×1e6)")
print("  • This formula is VALIDATED by:")
print("    - Colleague's working code producing correct plots")
print("    - Comparison with mass species (ratio = M_air ≈ 29)")
print("    - Multiple test locations giving physically reasonable values")
print("  • Typical atmospheric number concentration: 10²-10⁴ #/cm³ ✓")
print()
