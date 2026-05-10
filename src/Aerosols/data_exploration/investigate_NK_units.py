#!/usr/bin/env python3
"""
Comprehensive Investigation of NK Units and Meaning

Tests:
1. Does NK = sum of all mass species (in mol/mol)?
2. What other variables exist in the NetCDF file?
3. Convert NK from mol/mol to #/cm³ properly
4. Check consistency across different interpretations

Usage:
    python3 test/investigate_NK_units.py
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
R_gas = 8.314  # J/(mol·K)
N_A = 6.022e23  # Avogadro's number

print("=" * 80)
print("COMPREHENSIVE NK UNIT INVESTIGATION")
print("=" * 80)
print()

# Open dataset
ds = nc.Dataset(ncfile, 'r')

# Test location
idx, idy, idf = 2, 12, 4
lat = ds['lats'][idf, idy, idx]
lon = ds['lons'][idf, idy, idx]

print(f"📍 Location: Central USA ({lat:.2f}°N, {lon:.2f}°E)")
print()

# Get atmospheric structure
n_lev = ds.dimensions['lev'].size
dp = ds['Met_DELP'][0, :, idf, idy, idx]
sp = ds['Met_PS2WET'][0, idf, idy, idx]
pressure_half = np.concatenate([[sp], sp - np.cumsum(dp)])
pressure_half_toa2boa = pressure_half[::-1]
pressure_centers = (pressure_half_toa2boa[:-1] + pressure_half_toa2boa[1:]) / 2
pressure_Pa = pressure_centers * 100.0
temperature = ds['Met_T'][0, :, idf, idy, idx][::-1]

# Air number density
n_air_molpm3 = pressure_Pa / (R_gas * temperature)  # mol/m³
n_air_molpcm3 = n_air_molpm3 * 1e-6  # mol/cm³
n_air_molecules_cm3 = n_air_molpcm3 * N_A  # molecules/cm³

lev_800 = np.argmin(np.abs(pressure_centers - 800))
print(f"🌍 At ~800 hPa (level {lev_800}):")
print(f"   Pressure: {pressure_centers[lev_800]:.1f} hPa")
print(f"   Temperature: {temperature[lev_800]:.1f} K")
print(f"   Air density: {n_air_molpm3[lev_800]:.2e} mol/m³")
print(f"   Air density: {n_air_molpcm3[lev_800]:.2e} mol/cm³")
print(f"   Air density: {n_air_molecules_cm3[lev_800]:.2e} molecules/cm³")
print()

# =============================================================================
# TEST 1: Read NK and check its metadata
# =============================================================================
print("=" * 80)
print("TEST 1: NK Variable Metadata")
print("=" * 80)
print()

nk_var = ds['SpeciesConcVV_NK01']
print(f"Variable: SpeciesConcVV_NK01")
print(f"  Dimensions: {nk_var.dimensions}")
print(f"  Shape: {nk_var.shape}")
print(f"  Data type: {nk_var.dtype}")
print()
print("  Attributes:")
for attr in nk_var.ncattrs():
    print(f"    {attr}: {nk_var.getncattr(attr)}")
print()

# Read all NK bins
nk_concentration_raw = np.zeros((15, n_lev))
for bin_idx in range(1, 16):
    var_name = f"SpeciesConcVV_NK{bin_idx:02d}"
    if var_name in ds.variables:
        nk_concentration_raw[bin_idx-1, :] = ds[var_name][0, :, idf, idy, idx][::-1]

print(f"NK at 800 hPa (raw values from file):")
print(f"  Bin 1: {nk_concentration_raw[0, lev_800]:.6e}")
print(f"  Bin 5: {nk_concentration_raw[4, lev_800]:.6e}")
print(f"  Bin 10: {nk_concentration_raw[9, lev_800]:.6e}")
print(f"  Total (sum over bins): {nk_concentration_raw[:, lev_800].sum():.6e}")
print()

# =============================================================================
# TEST 2: Read all mass species and sum them
# =============================================================================
print("=" * 80)
print("TEST 2: Sum of All Mass-Based Species")
print("=" * 80)
print()

species_list = ["DUST", "SS", "SF", "ECIL", "ECOB", "OCIL", "OCOB", "AW"]
species_sum = np.zeros((15, n_lev))

for species in species_list:
    species_conc = np.zeros((15, n_lev))
    for bin_idx in range(1, 16):
        var_name = f"SpeciesConcVV_{species}{bin_idx:02d}"
        if var_name in ds.variables:
            species_conc[bin_idx-1, :] = ds[var_name][0, :, idf, idy, idx][::-1]
    species_sum += species_conc
    print(f"  {species:6s} at 800 hPa, bin 1: {species_conc[0, lev_800]:.6e}")

print()
print(f"Sum of all species at 800 hPa:")
print(f"  Bin 1: {species_sum[0, lev_800]:.6e}")
print(f"  Bin 5: {species_sum[4, lev_800]:.6e}")
print(f"  Bin 10: {species_sum[9, lev_800]:.6e}")
print(f"  Total (sum over bins): {species_sum[:, lev_800].sum():.6e}")
print()

# Compare NK vs sum of species
print("Comparison: NK vs Sum of Species (at 800 hPa)")
print("  Bin |        NK        |    Sum Species   |    Ratio (NK/Sum)")
print("  " + "-" * 65)
for bin_idx in range(15):
    nk_val = nk_concentration_raw[bin_idx, lev_800]
    sum_val = species_sum[bin_idx, lev_800]
    ratio = nk_val / sum_val if sum_val > 1e-30 else 0
    print(f"  {bin_idx+1:2d}  | {nk_val:16.6e} | {sum_val:16.6e} | {ratio:10.4f}")

print()
ratio_total = nk_concentration_raw[:, lev_800].sum() / species_sum[:, lev_800].sum()
print(f"Total ratio (NK/Sum): {ratio_total:.6f}")
print()

if 0.9 < ratio_total < 1.1:
    print("✅ NK ≈ Sum of species! NK is likely total aerosol in mol/mol units")
elif ratio_total < 0.1:
    print("❌ NK << Sum of species. NK might be something else entirely")
elif ratio_total > 10:
    print("❌ NK >> Sum of species. NK includes additional components")
else:
    print("⚠️  NK and Sum differ by factor of ~{:.1f}".format(ratio_total))

# =============================================================================
# TEST 3: Check what other variables exist in the file
# =============================================================================
print()
print("=" * 80)
print("TEST 3: Search for Other Relevant Variables")
print("=" * 80)
print()

print("Looking for variables that might be number concentration...")
print()

relevant_patterns = ['numb', 'num', 'NK', 'particle', 'count', 'N_']
relevant_vars = []

for var_name in ds.variables.keys():
    var_lower = var_name.lower()
    if any(pattern.lower() in var_lower for pattern in relevant_patterns):
        relevant_vars.append(var_name)

if relevant_vars:
    print("Found potentially relevant variables:")
    for var in relevant_vars[:20]:  # Limit to first 20
        var_obj = ds[var]
        units = var_obj.getncattr('units') if 'units' in var_obj.ncattrs() else 'no units'
        long_name = var_obj.getncattr('long_name') if 'long_name' in var_obj.ncattrs() else 'no description'
        print(f"  {var:40s} | {units:20s} | {long_name}")
else:
    print("No obvious number concentration variables found")

print()
print("Checking SpeciesConcVV prefix variables:")
speciesconcvv_vars = [v for v in ds.variables.keys() if v.startswith('SpeciesConcVV_')]
unique_species = set()
for var in speciesconcvv_vars:
    # Extract species name (remove bin number)
    species = var.replace('SpeciesConcVV_', '').rstrip('0123456789')
    unique_species.add(species)

print(f"  Found {len(unique_species)} unique species:")
for species in sorted(unique_species):
    print(f"    {species}")

# =============================================================================
# TEST 4: If NK is in mol/mol, convert it properly to #/cm³
# =============================================================================
print()
print("=" * 80)
print("TEST 4: Convert NK from mol/mol to #/cm³")
print("=" * 80)
print()

# Interpretation 1: NK is in mol/mol, directly convert
nk_number_cm3_v1 = nk_concentration_raw * n_air_molpcm3[:, np.newaxis].T * N_A

print(f"Interpretation 1: NK (mol/mol) × n_air (mol/cm³) × N_A")
print(f"  At 800 hPa:")
print(f"    Raw NK (mol/mol): {nk_concentration_raw[:, lev_800].sum():.3e}")
print(f"    Converted (#/cm³): {nk_number_cm3_v1[:, lev_800].sum():.3e}")
print()

# =============================================================================
# TEST 5: Create comparison plots
# =============================================================================
print("=" * 80)
print("Creating Diagnostic Plots...")
print("=" * 80)
print()

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: NK vs Sum of species (raw values)
ax = axes[0, 0]
bins = np.arange(1, 16)
width = 0.35
ax.bar(bins - width/2, nk_concentration_raw[:, lev_800], width, label='NK', alpha=0.8)
ax.bar(bins + width/2, species_sum[:, lev_800], width, label='Sum Species', alpha=0.8)
ax.set_xlabel('Bin Number', fontsize=11)
ax.set_ylabel('Concentration [mol/mol]', fontsize=11)
ax.set_yscale('log')
ax.set_title(f'NK vs Species Sum at {pressure_centers[lev_800]:.0f} hPa', fontsize=12, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')

# Plot 2: Ratio NK/Sum
ax = axes[0, 1]
ratio_by_bin = np.zeros(15)
for i in range(15):
    if species_sum[i, lev_800] > 1e-30:
        ratio_by_bin[i] = nk_concentration_raw[i, lev_800] / species_sum[i, lev_800]
ax.plot(bins, ratio_by_bin, marker='o', linewidth=2, markersize=8)
ax.axhline(1.0, color='red', linestyle='--', linewidth=2, alpha=0.7)
ax.set_xlabel('Bin Number', fontsize=11)
ax.set_ylabel('Ratio (NK / Sum Species)', fontsize=11)
ax.set_title('Ratio by Size Bin', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

# Plot 3: Vertical profile of NK vs Sum
ax = axes[1, 0]
nk_total_profile = nk_concentration_raw.sum(axis=0)
species_total_profile = species_sum.sum(axis=0)
ax.plot(nk_total_profile, pressure_centers, linewidth=2.5, label='NK', color='blue')
ax.plot(species_total_profile, pressure_centers, linewidth=2.5, label='Sum Species', color='red', linestyle='--')
ax.set_xlabel('Total Concentration [mol/mol]', fontsize=11)
ax.set_ylabel('Pressure [hPa]', fontsize=11)
ax.set_xscale('log')
ax.invert_yaxis()
ax.set_title('Vertical Profile: NK vs Species Sum', fontsize=12, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 4: Vertical profile of ratio
ax = axes[1, 1]
ratio_profile = np.zeros_like(nk_total_profile)
valid = species_total_profile > 1e-30
ratio_profile[valid] = nk_total_profile[valid] / species_total_profile[valid]
ax.plot(ratio_profile, pressure_centers, linewidth=2.5, color='purple')
ax.axvline(1.0, color='gray', linestyle='--', linewidth=2, alpha=0.7)
ax.set_xlabel('Ratio (NK / Sum Species)', fontsize=11)
ax.set_ylabel('Pressure [hPa]', fontsize=11)
ax.set_xlim([0, 2])
ax.invert_yaxis()
ax.set_title('Vertical Profile of Ratio', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / "investigation_NK_vs_species.png", dpi=150)
print("   Saved: investigation_NK_vs_species.png")
plt.close()

# =============================================================================
# FINAL SUMMARY
# =============================================================================
print()
print("=" * 80)
print("FINAL SUMMARY")
print("=" * 80)
print()

print("Key Findings:")
print()
print(f"1. NK metadata says units are: mol mol-1 dry")
print(f"2. NK/Species ratio at 800 hPa: {ratio_total:.4f}")
print()

if 0.95 < ratio_total < 1.05:
    print("🎯 CONCLUSION: NK = Sum of all aerosol species")
    print()
    print("   Interpretation:")
    print("   - NK is NOT number concentration in #/cm³")
    print("   - NK is total aerosol MASS in mol/mol units")
    print("   - NK is simply the sum: DUST + SS + SF + ECIL + ECOB + OCIL + OCOB + AW")
    print()
    print("   To get actual particle number concentration:")
    print("   - You need particle size distribution assumptions")
    print("   - Convert mass to number using density and bin sizes")
    print("   - Or look for a different variable with actual number data")
    print()
    print(f"   If we convert NK (mol/mol) → (#/cm³):")
    print(f"   NK × n_air × N_A = {nk_number_cm3_v1[:, lev_800].sum():.2e} particles/cm³")
    print(f"   (But this treats mol of aerosol as mol of particles - likely wrong!)")
    
elif ratio_total < 0.5:
    print("❓ CONCLUSION: NK < Sum of species")
    print()
    print("   NK might be:")
    print("   - A subset of species (missing some components)")
    print("   - Scaled differently")
    print("   - A different quantity entirely")
    
elif ratio_total > 2.0:
    print("❓ CONCLUSION: NK > Sum of species")
    print()
    print("   NK might include:")
    print("   - Additional aerosol types not in our species list")
    print("   - Water content counted differently")
    print("   - Other atmospheric components")

else:
    print(f"⚠️  CONCLUSION: NK and species differ by ~{ratio_total:.2f}×")
    print()
    print("   Could be:")
    print("   - Slight differences in what's included")
    print("   - Numerical precision issues")
    print("   - Different aggregation methods")

print()
print("Recommendation:")
print("   For optical calculations, use the individual species (DUST, SS, SF, etc.)")
print("   with their known densities and refractive indices, rather than NK.")
print()

ds.close()

print("=" * 80)
print("✅ Investigation Complete!")
print("=" * 80)
print()
