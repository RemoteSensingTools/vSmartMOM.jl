#!/usr/bin/env bash
set -euo pipefail

# Run one psurf=1000 hPa off-nadir VZA block for the expanded O2A Raman LUT.
#
# This wrapper intentionally uses a dedicated output directory per VZA block.
# That keeps restart checks isolated from the completed nadir LUT and from
# other VZA blocks.
#
# Required:
#   VZA_IDXS=2:3          # one-based indices into FULL_VZAS
#
# Typical:
#   CUDA_DEVICE=0 nohup workflows/OCORaman/scripts/run_o2a_psurf1000_vza_block.sh \
#       > /home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_vza_psurf1000/vza002-003/logs/nohup.out 2>&1 &

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${REPO_ROOT}"

: "${VZA_IDXS:?Set VZA_IDXS, e.g. VZA_IDXS=2:3}"

CUDA_DEVICE="${CUDA_DEVICE:-0}"
FULL_VZAS="${FULL_VZAS:-sza}"
RAZS="${RAZS:-0:15:345}"
CHUNK_ALBEDOS="${CHUNK_ALBEDOS:-21}"
JULIA_FLAGS="${JULIA_FLAGS:---pkgimages=no}"
OUTROOT="${OUTROOT:-/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_vza_psurf1000}"
TEE_CHUNK_LOGS="${TEE_CHUNK_LOGS:-0}"
# The generic chunk completeness checker is only keyed by psurf/SZA/albedo and
# is not VZA-aware.  Keep VZA blocks non-skipping by default; block allocation
# markers and locks prevent duplicate VZA work.
SKIP_COMPLETED="${SKIP_COMPLETED:-0}"

vza_tag() {
  printf 'vza%s' "$1" | tr ',' '_' | tr ':' '-'
}

TAG="$(vza_tag "${VZA_IDXS}")"
OUTDIR="${OUTDIR:-${OUTROOT}/${TAG}}"
LOGDIR="${LOGDIR:-${OUTDIR}/logs}"
mkdir -p "${OUTDIR}" "${LOGDIR}"

ALLOC_FILE="${OUTDIR}/.allocated_cuda_device"
if [[ -f "${ALLOC_FILE}" ]]; then
  ALLOCATED_DEVICE="$(tr -d '[:space:]' < "${ALLOC_FILE}")"
  if [[ -n "${ALLOCATED_DEVICE}" && "${ALLOCATED_DEVICE}" != "${CUDA_DEVICE}" ]]; then
    printf '[%s] Skipping psurf=1000 VZA block %s on CUDA device %s; allocated to CUDA device %s\n' \
      "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${TAG}" "${CUDA_DEVICE}" "${ALLOCATED_DEVICE}"
    exit 0
  fi
fi

LOCKDIR="${OUTDIR}/.run_lock"
if ! mkdir "${LOCKDIR}" 2>/dev/null; then
  printf '[%s] Skipping psurf=1000 VZA block %s on CUDA device %s; lock exists at %s\n' \
    "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${TAG}" "${CUDA_DEVICE}" "${LOCKDIR}"
  if [[ -f "${LOCKDIR}/owner.txt" ]]; then
    cat "${LOCKDIR}/owner.txt"
  fi
  exit 0
fi
{
  printf 'host=%s\n' "$(hostname)"
  printf 'pid=%s\n' "$$"
  printf 'started=%s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
  printf 'psurf=1000\n'
  printf 'vza_idxs=%s\n' "${VZA_IDXS}"
  printf 'cuda_device=%s\n' "${CUDA_DEVICE}"
} > "${LOCKDIR}/owner.txt"

STAMP="$(date +%Y%m%d_%H%M%S)"
MONLOG="${LOGDIR}/${TAG}_gpu${CUDA_DEVICE}_${STAMP}_memory.csv"

printf 'O2A psurf=1000 VZA-block run\n'
printf '  repo: %s\n' "${REPO_ROOT}"
printf '  VZA_IDXS: %s\n' "${VZA_IDXS}"
printf '  FULL_VZAS: %s\n' "${FULL_VZAS}"
printf '  RAZS: %s\n' "${RAZS}"
printf '  CUDA_DEVICE: %s\n' "${CUDA_DEVICE}"
printf '  OUTDIR: %s\n' "${OUTDIR}"
printf '  LOGDIR: %s\n' "${LOGDIR}"
printf '  CHUNK_ALBEDOS: %s\n' "${CHUNK_ALBEDOS}"
printf '  JULIA_FLAGS: %s\n' "${JULIA_FLAGS}"
printf '  SKIP_COMPLETED: %s\n' "${SKIP_COMPLETED}"
printf '  memory log: %s\n' "${MONLOG}"

(
  while true; do
    printf '%s,' "$(date '+%F %T %Z')"
    nvidia-smi --query-gpu=index,utilization.gpu,memory.used,memory.total,power.draw \
      --format=csv,noheader,nounits | tr '\n' ';'
    printf '\n'
    sleep 5
  done
) > "${MONLOG}" &
MONPID=$!

cleanup() {
  kill "${MONPID}" 2>/dev/null || true
  wait "${MONPID}" 2>/dev/null || true
  rm -rf "${LOCKDIR}"
}
trap cleanup EXIT INT TERM

env \
  PSURF=1000 \
  START_SZA_IDX=1 \
  START_ALBEDO_IDX=1 \
  END_SZA_IDX=14 \
  CHUNK_ALBEDOS="${CHUNK_ALBEDOS}" \
  SKIP_COMPLETED="${SKIP_COMPLETED}" \
  TEE_CHUNK_LOGS="${TEE_CHUNK_LOGS}" \
  CUDA_HOME="${CUDA_HOME:-/usr/local/cuda}" \
  CUDA_DEVICE="${CUDA_DEVICE}" \
  FULL_VZAS="${FULL_VZAS}" \
  VZA_IDXS="${VZA_IDXS}" \
  RAZS="${RAZS}" \
  OUTDIR="${OUTDIR}" \
  LOGDIR="${LOGDIR}" \
  JULIA_FLAGS="${JULIA_FLAGS}" \
  workflows/OCORaman/scripts/run_remaining_RamanLUT_O2A_chunked.sh

printf '\n[%s] Finished VZA block %s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${TAG}"
printf '  memory log: %s\n' "${MONLOG}"
