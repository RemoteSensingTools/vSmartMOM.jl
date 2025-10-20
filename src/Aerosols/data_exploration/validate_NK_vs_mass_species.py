#!/usr/bin/env python3
"""
Validation: NK Number Concentration vs Mass-Based Species

Tests whether NK (total number concentration) aligns with the sum of 
mass-based species converted to number concentration.

Strategy:
1. Read NK (number concentration #/cm³)
2. Read all mass species (DUST, SS, SF, ECIL, OCOB, OCIL, ECOB, AW) in mol/mol
3. Convert mass species to number concentration using:
   - Assumed particle density for each species
   - Assumed particle composition/mixing state
   - Bin size information
4. Compare total reconstructed number vs. NK

Requirements:
    pip install netCDF4 numpy matplotlib scipy

Usage:
    python3 test/validate_NK_vs_mass_species.py
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Configuration
ncfile = "GEOSChem.Custom.20190702_0000z.nc4"
output_dir = Path("test/aerosol_exploration_output")
output_dir.mkdir(exist_ok=True)

# Physical constants
R_gas = 8.314  # J/(mol·K) - Gas constant
N_A = 6.022e23  # Avogadro's number (molecules/mol)

# Species properties
species_info = {
    # Species: (density kg/m³, molar_mass kg/mol, description)
    "DUST": (2650.0, 0.100, "Mineral dust"),
    "SS": (2200.0, 0.058, "Sea salt (NaCl)"),
    "SF": (1770.0, 0.098, "Sulfate (H2SO4)"),
    "ECIL": (1800.0, 0.012, "BC hydrophilic"),
    "ECOB": (1800.0, 0.012, "BC hydrophobic"),
    "OCIL": (1400.0, 0.012, "OC hydrophilic"),
    "OCOB": (1400.0, 0.012, "OC hydrophobic"),
    "AW": (1000.0, 0.018, "Aerosol water"),
}

print("=" * 70)
print("NK vs Mass Species Validation")
print("=" * 70)
print()

# Open dataset
ds = nc.Dataset(ncfile, 'r')

# Test location
idx, idy, idf = 2, 12, 4
lat = ds['lats'][idf, idy, idx]
lon = ds['lons'][idf, idy, idx]
print(f"📍 Location: Central USA")
print(f"   Latitude:  {lat:.2f}°")
print(f"   Longitude: {lon:.2f}°")
print()

# TOMAS-15 bin definitions
diam_min_nm = 10.0
diam_max_nm = 10000.0
diam_edges_nm = diam_min_nm * (diam_max_nm / diam_min_nm) ** (np.arange(16) / 15)
diam_edges_m = diam_edges_nm * 1e-9  # Convert to meters
diam_centers_m = np.sqrt(diam_edges_m[:-1] * diam_edges_m[1:])
radius_centers_m = diam_centers_m / 2.0

# Volume of particles in each bin (assuming spherical)
V_particle = (4.0/3.0) * np.pi * radius_centers_m**3  # m³

print("🔬 Particle volumes per bin:")
for i in range(15):
    print(f"   Bin {i+1:2d}: r={radius_centers_m[i]*1e6:.4f} μm, V={V_particle[i]*1e18:.2e} μm³")
print()

# Get atmospheric structure
n_lev = ds.dimensions['lev'].size
dp = ds['Met_DELP'][0, :, idf, idy, idx]
sp = ds['Met_PS2WET'][0, idf, idy, idx]
pressure_half = np.concatenate([[sp], sp - np.cumsum(dp)])
pressure_half_toa2boa = pressure_half[::-1]
pressure_centers = (pressure_half_toa2boa[:-1] + pressure_half_toa2boa[1:]) / 2
pressure_Pa = pressure_centers * 100.0  # hPa to Pa
temperature = ds['Met_T'][0, :, idf, idy, idx][::-1]  # TOA→BOA

# Air number density (molecules/cm³)
n_air = pressure_Pa / (R_gas * temperature) * 1e-6  # mol/m³ → mol/cm³
n_air_molecules = n_air * N_A  # molecules/cm³

print("🌍 Atmospheric conditions at 800 hPa:")
lev_800 = np.argmin(np.abs(pressure_centers - 800))
print(f"   Pressure: {pressure_centers[lev_800]:.1f} hPa")
print(f"   Temperature: {temperature[lev_800]:.1f} K")
print(f"   Air number density: {n_air_molecules[lev_800]:.2e} molecules/cm³")
print()

# Read NK (actual number concentration)
print("📊 Reading NK (Total Number Concentration)...")
nk_concentration = np.zeros((15, n_lev))
for bin_idx in range(1, 16):
    var_name = f"SpeciesConcVV_NK{bin_idx:02d}"
    if var_name in ds.variables:
        nk_concentration[bin_idx-1, :] = ds[var_name][0, :, idf, idy, idx][::-1]

nk_total_per_layer = nk_concentration.sum(axis=0)
print(f"   NK total at 800 hPa: {nk_total_per_layer[lev_800]:.2e} #/cm³")
print()

# Read mass-based species and convert to number concentration
print("📊 Reading and converting mass-based species...")
species_number_concentration = {}

for species, (density, molar_mass, description) in species_info.items():
    print(f"\n   {species} ({description}):")
    print(f"      Density: {density} kg/m³, Molar mass: {molar_mass*1000} g/mol")
    
    # Read concentration in mol/mol for all bins
    species_conc = np.zeros((15, n_lev))
    for bin_idx in range(1, 16):
        var_name = f"SpeciesConcVV_{species}{bin_idx:02d}"
        if var_name in ds.variables:
            species_conc[bin_idx-1, :] = ds[var_name][0, :, idf, idy, idx][::-1]
    
    # Convert from volume mixing ratio (mol/mol) to number concentration (#/cm³)
    # Strategy: For each bin, calculate number of particles needed to produce
    # the observed mass, assuming particles fill the bin's size
    
    species_number = np.zeros((15, n_lev))
    
    for bin_idx in range(15):
        # Mass concentration in this bin [mol species / mol air]
        vmr = species_conc[bin_idx, :]  # mol/mol
        
        # Convert to mass concentration [kg species / m³ air]
        # vmr * n_air [mol/m³] * molar_mass [kg/mol] = kg/m³
        mass_conc_kg_m3 = vmr * n_air * molar_mass * 1e6  # mol/cm³ → mol/m³
        
        # Mass per particle assuming spherical and bin center size [kg]
        mass_per_particle = density * V_particle[bin_idx]  # kg/m³ * m³ = kg
        
        # Number of particles [#/cm³]
        if mass_per_particle > 0:
            species_number[bin_idx, :] = mass_conc_kg_m3 / mass_per_particle * 1e-6  # /m³ → /cm³
    
    species_number_concentration[species] = species_number
    
    total_at_800 = species_number[:, lev_800].sum()
    print(f"      Total number at 800 hPa: {total_at_800:.2e} #/cm³")
    if total_at_800 > 0:
        print(f"      Bin with max: {np.argmax(species_number[:, lev_800])+1}")

# Calculate total reconstructed number concentration
reconstructed_total = np.zeros((15, n_lev))
for species in species_info.keys():
    reconstructed_total += species_number_concentration[species]

reconstructed_total_per_layer = reconstructed_total.sum(axis=0)

print()
print("=" * 70)
print("Comparison at 800 hPa:")
print("=" * 70)
print(f"NK (observed):           {nk_total_per_layer[lev_800]:.3e} #/cm³")
print(f"Reconstructed (from mass): {reconstructed_total_per_layer[lev_800]:.3e} #/cm³")
ratio = reconstructed_total_per_layer[lev_800] / nk_total_per_layer[lev_800] if nk_total_per_layer[lev_800] > 0 else 0
print(f"Ratio (reconstructed/NK):  {ratio:.3f}")
print()

# Create comparison plots
print("=" * 70)
print("Creating Comparison Plots...")
print("=" * 70)
print()

# Plot 1: Vertical profiles
print("📈 Plot 1: Vertical Profile Comparison")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Left: Absolute comparison
ax1.plot(nk_total_per_layer, pressure_centers, 
         linewidth=3, label='NK (observed)', color='black')
ax1.plot(reconstructed_total_per_layer, pressure_centers, 
         linewidth=3, label='Reconstructed from mass', color='red', linestyle='--')
ax1.set_xlabel('Total Number Concentration [#/cm³]', fontsize=13)
ax1.set_ylabel('Pressure [hPa]', fontsize=13)
ax1.set_xscale('log')
ax1.invert_yaxis()
ax1.set_title('Absolute Number Concentration', fontsize=14, fontweight='bold')
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)

# Right: Ratio
ratio_profile = np.zeros_like(nk_total_per_layer)
valid = nk_total_per_layer > 0
ratio_profile[valid] = reconstructed_total_per_layer[valid] / nk_total_per_layer[valid]

ax2.plot(ratio_profile, pressure_centers, linewidth=3, color='blue')
ax2.axvline(1.0, color='gray', linestyle='--', linewidth=2, alpha=0.7)
ax2.set_xlabel('Ratio (Reconstructed / NK)', fontsize=13)
ax2.set_ylabel('Pressure [hPa]', fontsize=13)
ax2.set_xscale('log')
ax2.invert_yaxis()
ax2.set_title('Ratio of Reconstructed to Observed', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / "validation_01_vertical_comparison.png", dpi=150)
print("   Saved: validation_01_vertical_comparison.png")
plt.close()

# Plot 2: Size distribution comparison at 800 hPa
print("📈 Plot 2: Size Distribution at 800 hPa")
fig, ax = plt.subplots(figsize=(12, 7))

diam_centers_um = diam_centers_m * 1e6

ax.plot(diam_centers_um, nk_concentration[:, lev_800], 
        marker='o', linewidth=3, markersize=8, 
        label='NK (observed)', color='black')
ax.plot(diam_centers_um, reconstructed_total[:, lev_800], 
        marker='s', linewidth=3, markersize=8, 
        label='Reconstructed from mass', color='red', linestyle='--')

ax.set_xlabel('Diameter [μm]', fontsize=13)
ax.set_ylabel('Number Concentration [#/cm³]', fontsize=13)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title(f'Size Distribution at {pressure_centers[lev_800]:.0f} hPa', 
             fontsize=15, fontweight='bold')
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3, which='both')
plt.tight_layout()
plt.savefig(output_dir / "validation_02_size_distribution.png", dpi=150)
print("   Saved: validation_02_size_distribution.png")
plt.close()

# Plot 3: Species contribution breakdown
print("📈 Plot 3: Species Contribution at 800 hPa")
fig, ax = plt.subplots(figsize=(12, 7))

bottom = np.zeros(15)
colors = plt.cm.tab10(np.linspace(0, 0.9, len(species_info)))

for i, (species, (_, _, desc)) in enumerate(species_info.items()):
    values = species_number_concentration[species][:, lev_800]
    ax.bar(np.arange(15), values, bottom=bottom, 
           label=f'{species} ({desc})', color=colors[i], alpha=0.8)
    bottom += values

# Overlay NK
ax.plot(np.arange(15), nk_concentration[:, lev_800], 
        marker='o', linewidth=3, markersize=8, color='black', 
        label='NK (total observed)', linestyle='none')

ax.set_xlabel('Bin Number', fontsize=13)
ax.set_ylabel('Number Concentration [#/cm³]', fontsize=13)
ax.set_yscale('log')
ax.set_title(f'Species Breakdown at {pressure_centers[lev_800]:.0f} hPa', 
             fontsize=15, fontweight='bold')
ax.legend(fontsize=9, ncol=2)
ax.grid(True, alpha=0.3, axis='y')
plt.tight_layout()
plt.savefig(output_dir / "validation_03_species_breakdown.png", dpi=150)
print("   Saved: validation_03_species_breakdown.png")
plt.close()

# Statistical summary
print()
print("=" * 70)
print("Statistical Summary:")
print("=" * 70)

# Boundary layer (surface to 700 hPa)
bl_levels = pressure_centers > 700
nk_bl = nk_total_per_layer[bl_levels]
rec_bl = reconstructed_total_per_layer[bl_levels]
ratio_bl = rec_bl / nk_bl
ratio_bl = ratio_bl[np.isfinite(ratio_bl)]

print(f"\nBoundary Layer (>700 hPa):")
print(f"   Mean ratio: {ratio_bl.mean():.3f}")
print(f"   Std ratio:  {ratio_bl.std():.3f}")
print(f"   Min ratio:  {ratio_bl.min():.3f}")
print(f"   Max ratio:  {ratio_bl.max():.3f}")

# Free troposphere (300-700 hPa)
ft_levels = (pressure_centers >= 300) & (pressure_centers <= 700)
nk_ft = nk_total_per_layer[ft_levels]
rec_ft = reconstructed_total_per_layer[ft_levels]
ratio_ft = rec_ft / nk_ft
ratio_ft = ratio_ft[np.isfinite(ratio_ft)]

print(f"\nFree Troposphere (300-700 hPa):")
print(f"   Mean ratio: {ratio_ft.mean():.3f}")
print(f"   Std ratio:  {ratio_ft.std():.3f}")
print(f"   Min ratio:  {ratio_ft.min():.3f}")
print(f"   Max ratio:  {ratio_ft.max():.3f}")

# Upper troposphere (<300 hPa)
ut_levels = pressure_centers < 300
nk_ut = nk_total_per_layer[ut_levels]
rec_ut = reconstructed_total_per_layer[ut_levels]
ratio_ut = rec_ut / nk_ut
ratio_ut = ratio_ut[np.isfinite(ratio_ut)]

print(f"\nUpper Troposphere (<300 hPa):")
print(f"   Mean ratio: {ratio_ut.mean():.3f}")
print(f"   Std ratio:  {ratio_ut.std():.3f}")
print(f"   Min ratio:  {ratio_ut.min():.3f}")
print(f"   Max ratio:  {ratio_ut.max():.3f}")

print()
print("=" * 70)
print("Interpretation:")
print("=" * 70)
print()
if ratio_bl.mean() < 0.1:
    print("⚠️  MISMATCH: Reconstructed number is much smaller than NK")
    print("   Possible causes:")
    print("   1. NK units might actually be mol/mol (not #/cm³)")
    print("   2. Density/molar mass assumptions are incorrect")
    print("   3. Particle size assumptions (bin center) are wrong")
    print("   4. NK includes aerosol types not in mass species")
elif ratio_bl.mean() > 10:
    print("⚠️  MISMATCH: Reconstructed number is much larger than NK")
    print("   Possible causes:")
    print("   1. Double counting between species")
    print("   2. Mass species include water that shouldn't be counted")
    print("   3. Particle size assumptions are wrong")
elif 0.5 < ratio_bl.mean() < 2.0:
    print("✅ GOOD AGREEMENT: NK and mass-based species are consistent!")
    print("   This suggests:")
    print("   1. NK units are likely #/cm³")
    print("   2. Density/size assumptions are reasonable")
    print("   3. Species coverage is comprehensive")
else:
    print("⚠️  PARTIAL MATCH: Some discrepancy exists")
    print(f"   Ratio of {ratio_bl.mean():.2f} suggests further investigation")

ds.close()

print()
print("=" * 70)
print("✅ Validation Complete!")
print("=" * 70)
print()
