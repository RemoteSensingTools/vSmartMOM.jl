#!/usr/bin/env bash
set -euo pipefail

# Run one nadir-only O2 singlet Raman LUT pressure slice.
#
# Required/typical environment:
#   PSURF=1000 CUDA_DEVICE=1 workflows/OCORaman/scripts/run_o2singlet_nadir_psurf.sh
#
# Optional:
#   OUTDIR=/home/sanghavi/data/RamanSIFgrid/O2Singlet
#   JULIA_FLAGS="--pkgimages=no"
#   FORCE=1      # overwrite existing NetCDF
#   DRY_RUN=1    # print plan only

PSURF="${PSURF:-1000}"
CUDA_DEVICE="${CUDA_DEVICE:-1}"
ARCH="${ARCH:-GPU}"
JULIA_FLAGS="${JULIA_FLAGS:---pkgimages=no}"

case "${PSURF}" in
  1000|1000.0) PSURF_TAG="1000" ;;
  750|750.0)   PSURF_TAG="0750" ;;
  500|500.0)   PSURF_TAG="0500" ;;
  *)           PSURF_TAG="$(printf "%04.0f" "${PSURF}")" ;;
esac

ROOT="${OUTDIR:-/home/sanghavi/data/RamanSIFgrid/O2Singlet}"
mkdir -p "${ROOT}/logs"

OUT_NC="${OUT_NC:-${ROOT}/o2singlet_raman_lut_psurf${PSURF_TAG}_nadir.nc}"
LOG="${LOG:-${ROOT}/logs/o2singlet_raman_lut_psurf${PSURF_TAG}_nadir_$(date +%Y%m%d_%H%M%S).log}"

# Keep exactly one writer attached to each NetCDF. Concurrent writers can leave
# an apparently valid filename pointing at stale or incomplete HDF5 contents.
LOCK_FILE="${OUT_NC}.lock"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
  echo "ERROR: another process already holds the output lock: ${LOCK_FILE}" >&2
  exit 2
fi

RESUME_FLAG="${RESUME:-0}"
if [[ -f "${OUT_NC}" && "${FORCE:-0}" != "1" ]]; then
  RESUME_FLAG="1"
fi

echo "O2 singlet nadir LUT run"
echo "  psurf       = ${PSURF}"
echo "  CUDA_DEVICE = ${CUDA_DEVICE}"
echo "  OUT_NC      = ${OUT_NC}"
echo "  LOG         = ${LOG}"
echo "  LOCK_FILE   = ${LOCK_FILE}"
echo "  RESUME      = ${RESUME_FLAG}"

env \
  PSURFS="${PSURF}" \
  OUTDIR="${ROOT}" \
  OUT_NC="${OUT_NC}" \
  CUDA_DEVICE="${CUDA_DEVICE}" \
  ARCH="${ARCH}" \
  RESUME="${RESUME_FLAG}" \
  FORCE="${FORCE:-0}" \
  DRY_RUN="${DRY_RUN:-0}" \
  SIF_STATES="0" \
  FULL_VZAS="0.0" \
  VZA_IDXS="1" \
  RAZS="0.0" \
  julia ${JULIA_FLAGS} --project=. workflows/OCORaman/scripts/createRamanLUT_O2Singlet.jl \
  2>&1 | tee "${LOG}"
