# NetCDF I/O utilities for vSmartMOM
# Provides shared functionality for reading NetCDF-based atmospheric data

module NetCDFIO

using NCDatasets

# Access parent modules
const Formats = Base.parentmodule(Base.parentmodule(@__MODULE__)).Formats
const Sources = Base.parentmodule(Base.parentmodule(@__MODULE__)).Sources
const CoreRT = Base.parentmodule(Base.parentmodule(Base.parentmodule(@__MODULE__))).CoreRT

using .Formats: IOSource
import .Formats: load_config
using .Sources: NetCDFSource, GeosChemSource, NetCDFGridSource
using .CoreRT: compute_atmos_profile_fields

# Include format-specific readers
include("GeosChem.jl")

# Re-export for convenience
export NetCDFSource, GeosChemSource, NetCDFGridSource
export geoschem_to_dict, read_geoschem_profile
export load_config  # Re-export the extended method

end # module NetCDFIO
