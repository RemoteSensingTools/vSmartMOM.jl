#!/usr/bin/env bash
set -euo pipefail

# Sequentially run several off-nadir VZA blocks for a selected pressure.
#
# Required:
#   PSURF=750
#   VZA_BLOCKS="1 2 3 4"
#
# Typical:
#   CUDA_DEVICE=1 PSURF=750 VZA_BLOCKS="$(seq -s ' ' 1 14)" \
#     workflows/OCORaman/scripts/run_o2a_psurf_vza_blocks.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${REPO_ROOT}"

: "${PSURF:?Set PSURF, e.g. PSURF=750 or PSURF=1000}"
: "${VZA_BLOCKS:?Set VZA_BLOCKS, e.g. VZA_BLOCKS=\"1 2 3\"}"

CUDA_DEVICE="${CUDA_DEVICE:-0}"
CHUNK_ALBEDOS="${CHUNK_ALBEDOS:-21}"
JULIA_FLAGS="${JULIA_FLAGS:---pkgimages=no}"
TEE_CHUNK_LOGS="${TEE_CHUNK_LOGS:-0}"

printf -v PSURF_TAG '%04.0f' "${PSURF}"
OUTROOT="${OUTROOT:-/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_vza_psurf${PSURF_TAG}}"
mkdir -p "${OUTROOT}/logs"

printf 'O2A psurf=%s multi-VZA-block run\n' "${PSURF}"
printf '  repo: %s\n' "${REPO_ROOT}"
printf '  CUDA_DEVICE: %s\n' "${CUDA_DEVICE}"
printf '  VZA_BLOCKS: %s\n' "${VZA_BLOCKS}"
printf '  OUTROOT: %s\n' "${OUTROOT}"
printf '  CHUNK_ALBEDOS: %s\n' "${CHUNK_ALBEDOS}"
printf '  JULIA_FLAGS: %s\n' "${JULIA_FLAGS}"

for block in ${VZA_BLOCKS}; do
  printf '\n[%s] Starting psurf=%s VZA block %s on CUDA device %s\n' \
    "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${PSURF}" "${block}" "${CUDA_DEVICE}"
  env \
    PSURF="${PSURF}" \
    VZA_IDXS="${block}" \
    CUDA_DEVICE="${CUDA_DEVICE}" \
    OUTROOT="${OUTROOT}" \
    CHUNK_ALBEDOS="${CHUNK_ALBEDOS}" \
    JULIA_FLAGS="${JULIA_FLAGS}" \
    TEE_CHUNK_LOGS="${TEE_CHUNK_LOGS}" \
    workflows/OCORaman/scripts/run_o2a_psurf_vza_block.sh
  printf '[%s] Finished psurf=%s VZA block %s on CUDA device %s\n' \
    "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${PSURF}" "${block}" "${CUDA_DEVICE}"
done

printf '\n[%s] Finished all psurf=%s VZA blocks on CUDA device %s\n' \
  "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${PSURF}" "${CUDA_DEVICE}"
