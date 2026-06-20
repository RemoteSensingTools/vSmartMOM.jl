#!/usr/bin/env bash
set -euo pipefail

# Remainder wrapper for the Curry psurf=750 hPa slice.
#
# Progress snapshot from 2026-06-18 14:45 PDT:
#   complete: SZA indices 1:2, all 21 albedos
#             SZA index 3, albedo indices 1:14 (rho=0.00:0.65)
#   first remaining at snapshot: SZA index 3, albedo index 15 (rho=0.70)
#
# If you wait until the currently running rho=0.70 scene finishes cleanly,
# launch with START_ALBEDO_IDX=16 to avoid duplicating that scene.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export PSURF="${PSURF:-750}"
export CUDA_DEVICE="${CUDA_DEVICE:-1}"
export START_SZA_IDX="${START_SZA_IDX:-3}"
export START_ALBEDO_IDX="${START_ALBEDO_IDX:-15}"
export END_SZA_IDX="${END_SZA_IDX:-14}"
export CHUNK_ALBEDOS="${CHUNK_ALBEDOS:-4}"
export OUTDIR="${OUTDIR:-/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_chunked_psurf0750}"

exec "${SCRIPT_DIR}/run_remaining_RamanLUT_O2A_chunked.sh"
