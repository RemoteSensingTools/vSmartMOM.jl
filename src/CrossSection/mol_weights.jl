# File that stores all isotope information
const file_path = String(@__DIR__) * "/iso_info.nc"

# Read in all the variables beforehand
const global_ids = ncread(file_path, "global_id")
const isonames = ncread(file_path, "isoname")
const abundances = ncread(file_path, "abundance")
const weights = ncread(file_path, "mol_weight")
const names = ncread(file_path, "mol_name")

# Check that the value is not an "unfilled" value (using -1 and "" as defaults)
check_exists(value) = (value in (-1, "")) ? error("No matching (mol, iso) pair") : value

# Functions to get molecule/isotope information 
mol_globalID(mol, iso) = check_exists(global_ids[mol, iso])
mol_isoname(mol, iso) = check_exists(isonames[mol, iso])
mol_abundance(mol, iso) = check_exists(abundances[mol, iso])
mol_weight(mol, iso) = check_exists(weights[mol, iso])
mol_name(mol, iso) = check_exists(names[mol, iso])
