#!/usr/bin/env bash
set -euo pipefail

# Sequentially run several psurf=1000 hPa off-nadir VZA blocks on one GPU.
#
# Required:
#   VZA_BLOCKS="2:3 4:5 6:7"
#
# Typical:
#   CUDA_DEVICE=0 VZA_BLOCKS="2:3 4:5 6:7" \
#     workflows/OCORaman/scripts/run_o2a_psurf1000_vza_blocks.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${REPO_ROOT}"

: "${VZA_BLOCKS:?Set VZA_BLOCKS, e.g. VZA_BLOCKS=\"2:3 4:5 6:7\"}"

CUDA_DEVICE="${CUDA_DEVICE:-0}"
OUTROOT="${OUTROOT:-/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_vza_psurf1000}"
CHUNK_ALBEDOS="${CHUNK_ALBEDOS:-21}"
JULIA_FLAGS="${JULIA_FLAGS:---pkgimages=no}"
TEE_CHUNK_LOGS="${TEE_CHUNK_LOGS:-0}"

mkdir -p "${OUTROOT}/logs"

printf 'O2A psurf=1000 multi-VZA-block run\n'
printf '  repo: %s\n' "${REPO_ROOT}"
printf '  CUDA_DEVICE: %s\n' "${CUDA_DEVICE}"
printf '  VZA_BLOCKS: %s\n' "${VZA_BLOCKS}"
printf '  OUTROOT: %s\n' "${OUTROOT}"
printf '  CHUNK_ALBEDOS: %s\n' "${CHUNK_ALBEDOS}"
printf '  JULIA_FLAGS: %s\n' "${JULIA_FLAGS}"

for block in ${VZA_BLOCKS}; do
  printf '\n[%s] Starting VZA block %s on CUDA device %s\n' \
    "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${block}" "${CUDA_DEVICE}"
  env \
    VZA_IDXS="${block}" \
    CUDA_DEVICE="${CUDA_DEVICE}" \
    OUTROOT="${OUTROOT}" \
    CHUNK_ALBEDOS="${CHUNK_ALBEDOS}" \
    JULIA_FLAGS="${JULIA_FLAGS}" \
    TEE_CHUNK_LOGS="${TEE_CHUNK_LOGS}" \
    workflows/OCORaman/scripts/run_o2a_psurf1000_vza_block.sh
  printf '[%s] Finished VZA block %s on CUDA device %s\n' \
    "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${block}" "${CUDA_DEVICE}"
done

printf '\n[%s] Finished all VZA blocks on CUDA device %s\n' \
  "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${CUDA_DEVICE}"
