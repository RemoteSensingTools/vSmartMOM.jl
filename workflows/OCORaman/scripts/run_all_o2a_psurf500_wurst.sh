#!/usr/bin/env bash
set -euo pipefail

# Full chunked wrapper for the Wurst psurf=500 hPa slice.
#
# Intended use:
#   Start this on Wurst on GPU 0, in parallel with the psurf=1000 hPa run on
#   GPU 1, or after GPU 1 is free.
#
# This is a full-slice run: all 14 SZA nodes and all 21 albedos.  The shared
# chunk launcher still checks completed scenes before each chunk, so restarting
# this script will skip any valid scenes already written.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export PSURF="${PSURF:-500}"
export CUDA_DEVICE="${CUDA_DEVICE:-0}"
export START_SZA_IDX="${START_SZA_IDX:-1}"
export START_ALBEDO_IDX="${START_ALBEDO_IDX:-1}"
export END_SZA_IDX="${END_SZA_IDX:-14}"
export CHUNK_ALBEDOS="${CHUNK_ALBEDOS:-4}"
export SKIP_COMPLETED="${SKIP_COMPLETED:-1}"
export TEE_CHUNK_LOGS="${TEE_CHUNK_LOGS:-0}"
export OUTDIR="${OUTDIR:-/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_chunked_psurf0500}"

exec "${SCRIPT_DIR}/run_remaining_RamanLUT_O2A_chunked.sh"
