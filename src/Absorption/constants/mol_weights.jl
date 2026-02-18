#=
 
This file contains helper functions to load and store isotope information relevant 
to absorption cross section calculation
 
=#

# File that stores all isotope information
const file_path = String(@__DIR__) * "/iso_info.nc"

# Read in all the variables beforehand
const mols = ncread(file_path, "molecule")
const isos = ncread(file_path, "isotope")
const global_ids = ncread(file_path, "global_id")
const isonames = ncread(file_path, "isoname")
const abundances = ncread(file_path, "abundance")
const weights = ncread(file_path, "mol_weight")
const mol_names = ncread(file_path, "mol_name")

"""Throw if value is unfilled (-1 or ""). Used for invalid (mol, iso) lookups."""
check_exists(value) = (value in (-1, "")) ? error("No matching (mol, iso) pair") : value

"""Get HITRAN global ID for molecule mol and isotopologue iso."""
mol_globalID(mol, iso) = check_exists(global_ids[mol, iso])

"""Get isotopologue name (e.g. \"1H2-16O\") for molecule mol and isotopologue iso."""
mol_isoname(mol, iso) = check_exists(isonames[mol, iso])

"""Get natural abundance for molecule mol and isotopologue iso."""
mol_abundance(mol, iso) = check_exists(abundances[mol, iso])

"""Get molar mass (kg/mol) for molecule mol and isotopologue iso."""
mol_weight(mol, iso) = check_exists(weights[mol, iso])

"""Get molecule name for molecule mol and isotopologue iso."""
mol_name(mol, iso) = check_exists(names[mol, iso])

"""Get HITRAN molecule ID from molecule name string."""
mol_number(name) = mols[findfirst(x->x==name, mol_names[:,1])]

"""List all unique molecule names in the isotope database."""
show_molecules() = unique(mol_names)

"""Find molecule names containing search_str."""
search_molecules(search_str) = filter(x->occursin(search_str, x), unique(mol_names)) 