#!/usr/bin/env bash
set -euo pipefail

# Remainder wrapper for the Wurst psurf=1000 hPa slice.
#
# Progress snapshot from 2026-06-18 14:45 PDT:
#   complete: SZA indices 1:3, all 21 albedos
#   first remaining: SZA index 4, albedo index 1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export PSURF="${PSURF:-1000}"
export CUDA_DEVICE="${CUDA_DEVICE:-1}"
export START_SZA_IDX="${START_SZA_IDX:-4}"
export START_ALBEDO_IDX="${START_ALBEDO_IDX:-1}"
export END_SZA_IDX="${END_SZA_IDX:-14}"
export CHUNK_ALBEDOS="${CHUNK_ALBEDOS:-4}"
export OUTDIR="${OUTDIR:-/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_chunked_psurf1000}"

exec "${SCRIPT_DIR}/run_remaining_RamanLUT_O2A_chunked.sh"
