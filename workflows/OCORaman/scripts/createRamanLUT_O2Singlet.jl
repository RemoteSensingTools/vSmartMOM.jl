#!/usr/bin/env julia

# O2 1.27 micron singlet Raman/Cabannes/Rayleigh LUT driver.
#
# This wrapper sets band-specific defaults and then reuses the chunk-capable
# O2A LUT driver. All usual chunk controls remain available:
#   PSURFS, SZA_IDXS, ALBEDO_IDXS, FULL_VZAS, VZA_IDXS, RAZS, OUTDIR, OUT_NC.
#
# Dry-run example:
#   CUDA_DEVICE=1 DRY_RUN=1 julia --project=. \
#       workflows/OCORaman/scripts/createRamanLUT_O2Singlet.jl
#
# One pressure slice example:
#   CUDA_DEVICE=1 PSURFS=1000 \
#       julia --project=. workflows/OCORaman/scripts/createRamanLUT_O2Singlet.jl

function set_default_env!(key::AbstractString, value::AbstractString)
    haskey(ENV, key) || (ENV[key] = value)
    return ENV[key]
end

const _SCRIPT_DIR = @__DIR__
const _REPO_ROOT = normpath(joinpath(_SCRIPT_DIR, "..", "..", ".."))
const _HITRAN_LUT_DIR = set_default_env!(
    "VSMARTMOM_HITRAN_LUT_DIR", "/home/sanghavi/data/HITRAN_LUTs")

set_default_env!("SCRIPT_PATH_OVERRIDE", @__FILE__)
set_default_env!("BAND_TAG", "O2Singlet")
set_default_env!("BAND_TITLE", "O2 1.27 micron singlet band (1265-1275 nm)")
set_default_env!("BAND_OUTPUT_PREFIX", "RamanLUT_O2Singlet_float32_gpu")
set_default_env!("PARAM_YAML", joinpath(
    _REPO_ROOT, "workflows", "OCORaman", "config", "paramsRamanLUT_O2Singlet.yaml"))
set_default_env!("WL_EDGES_NM", "1275,1265")
set_default_env!("OUTDIR", joinpath(homedir(), "data", "RamanSIFgrid", "O2Singlet"))
set_default_env!("FULL_VZAS", "0.0")
set_default_env!("VZA_IDXS", "1")
set_default_env!("RAZS", "0.0")

# The local O2 ABSCO v5.2 file spans only the A-band, so use HITRAN O2 here.
set_default_env!("USE_O2_ABSCO", "0")
set_default_env!("O2_LUT", joinpath(_HITRAN_LUT_DIR, "O2.jld2"))
set_default_env!("CO2_LUT", joinpath(_HITRAN_LUT_DIR, "CO2.jld2"))
set_default_env!("H2O_LUT", joinpath(_HITRAN_LUT_DIR, "H2O.jld2"))
set_default_env!("FIXED_MOLECULES", "O2,CO2")
set_default_env!("VARIABLE_MOLECULES", "")

include(joinpath(_SCRIPT_DIR, "chunked_createRamanLUT_O2A.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    run!()
end
