#!/usr/bin/env bash
set -euo pipefail

# Run one off-nadir VZA block for the expanded O2A Raman LUT.
#
# Required:
#   PSURF=750
#   VZA_IDXS=1          # one-based indices into FULL_VZAS; use one index for
#                       # single-VZA memory-conservative runs
#
# Typical:
#   CUDA_DEVICE=1 PSURF=750 VZA_IDXS=1 \
#     workflows/OCORaman/scripts/run_o2a_psurf_vza_block.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${REPO_ROOT}"

: "${PSURF:?Set PSURF, e.g. PSURF=750 or PSURF=1000}"
: "${VZA_IDXS:?Set VZA_IDXS, e.g. VZA_IDXS=1 or VZA_IDXS=2:3}"

CUDA_DEVICE="${CUDA_DEVICE:-0}"
FULL_VZAS="${FULL_VZAS:-sza}"
RAZS="${RAZS:-0:15:345}"
CHUNK_ALBEDOS="${CHUNK_ALBEDOS:-21}"
JULIA_FLAGS="${JULIA_FLAGS:---pkgimages=no}"
TEE_CHUNK_LOGS="${TEE_CHUNK_LOGS:-0}"
ALLOW_NADIR_VZA="${ALLOW_NADIR_VZA:-0}"
# The generic chunk completeness checker is only keyed by psurf/SZA/albedo and
# is not VZA-aware.  Keep VZA blocks non-skipping by default; block allocation
# markers and locks prevent duplicate VZA work.
SKIP_COMPLETED="${SKIP_COMPLETED:-0}"

if [[ "${ALLOW_NADIR_VZA}" != "1" ]]; then
  IFS=',:' read -ra _vza_parts <<< "${VZA_IDXS}"
  for _vza_part in "${_vza_parts[@]}"; do
    if [[ "${_vza_part}" == "1" ]]; then
      cat >&2 <<'EOF'
Refusing to run VZA index 1 because it is the nadir-viewing case already
covered by the nadir LUT. Set ALLOW_NADIR_VZA=1 only for an intentional
diagnostic recomputation.
EOF
      exit 2
    fi
  done
fi

printf -v PSURF_TAG '%04.0f' "${PSURF}"
OUTROOT="${OUTROOT:-/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_vza_psurf${PSURF_TAG}}"

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
    printf '[%s] Skipping psurf=%s VZA block %s on CUDA device %s; allocated to CUDA device %s\n' \
      "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${PSURF}" "${TAG}" "${CUDA_DEVICE}" "${ALLOCATED_DEVICE}"
    exit 0
  fi
fi

ALLOC_HOST_FILE="${OUTDIR}/.allocated_host"
if [[ -f "${ALLOC_HOST_FILE}" ]]; then
  ALLOCATED_HOST="$(tr -d '[:space:]' < "${ALLOC_HOST_FILE}")"
  CURRENT_HOST="$(hostname)"
  if [[ -n "${ALLOCATED_HOST}" && "${ALLOCATED_HOST}" != "${CURRENT_HOST}" ]]; then
    printf '[%s] Skipping psurf=%s VZA block %s on host %s; allocated to host %s\n' \
      "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${PSURF}" "${TAG}" "${CURRENT_HOST}" "${ALLOCATED_HOST}"
    exit 0
  fi
fi

LOCKDIR="${OUTDIR}/.run_lock"
if ! mkdir "${LOCKDIR}" 2>/dev/null; then
  printf '[%s] Skipping psurf=%s VZA block %s on CUDA device %s; lock exists at %s\n' \
    "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${PSURF}" "${TAG}" "${CUDA_DEVICE}" "${LOCKDIR}"
  if [[ -f "${LOCKDIR}/owner.txt" ]]; then
    cat "${LOCKDIR}/owner.txt"
  fi
  exit 0
fi
{
  printf 'host=%s\n' "$(hostname)"
  printf 'pid=%s\n' "$$"
  printf 'started=%s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
  printf 'psurf=%s\n' "${PSURF}"
  printf 'vza_idxs=%s\n' "${VZA_IDXS}"
  printf 'cuda_device=%s\n' "${CUDA_DEVICE}"
} > "${LOCKDIR}/owner.txt"

STAMP="$(date +%Y%m%d_%H%M%S)"
MONLOG="${LOGDIR}/${TAG}_psurf${PSURF_TAG}_gpu${CUDA_DEVICE}_${STAMP}_memory.csv"

printf 'O2A psurf=%s VZA-block run\n' "${PSURF}"
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

export PSURF
export START_SZA_IDX=1
export START_ALBEDO_IDX=1
export END_SZA_IDX=14
export CHUNK_ALBEDOS
export SKIP_COMPLETED
export TEE_CHUNK_LOGS
export CUDA_HOME="${CUDA_HOME:-/usr/local/cuda}"
export CUDA_DEVICE
export FULL_VZAS
export VZA_IDXS
export RAZS
export OUTDIR
export LOGDIR
export JULIA_FLAGS

workflows/OCORaman/scripts/run_remaining_RamanLUT_O2A_chunked.sh

printf '\n[%s] Finished psurf=%s VZA block %s\n' \
  "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${PSURF}" "${TAG}"
printf '  memory log: %s\n' "${MONLOG}"
