#!/usr/bin/env bash
set -euo pipefail

# Wait for an already-running pressure slice, verify that it completed, then
# run the remaining O2-singlet nadir pressure slices sequentially.
#
# Example for a currently running 1000 hPa process:
#   WAIT_PID=12345 \
#   WAIT_LOG=/path/to/o2singlet_raman_lut_psurf1000_nadir_TIMESTAMP.log \
#   CUDA_DEVICE=1 \
#   nohup workflows/OCORaman/scripts/run_o2singlet_nadir_sequence.sh \
#     > /path/to/logs/o2singlet_nadir_sequence.out 2>&1 &

WAIT_PID="${WAIT_PID:-}"
WAIT_LOG="${WAIT_LOG:-}"
POLL_SECONDS="${POLL_SECONDS:-60}"
CUDA_DEVICE="${CUDA_DEVICE:-1}"
PSURF_SEQUENCE="${PSURF_SEQUENCE:-750 500}"
OUTDIR="${OUTDIR:-/home/sanghavi/data/RamanSIFgrid/O2Singlet}"
LAUNCHER="workflows/OCORaman/scripts/run_o2singlet_nadir_psurf.sh"

[[ -n "${WAIT_PID}" ]] || { echo "ERROR: WAIT_PID is required" >&2; exit 2; }
[[ "${WAIT_PID}" =~ ^[0-9]+$ ]] || { echo "ERROR: WAIT_PID must be numeric" >&2; exit 2; }
[[ -n "${WAIT_LOG}" ]] || { echo "ERROR: WAIT_LOG is required" >&2; exit 2; }
[[ -f "${WAIT_LOG}" ]] || { echo "ERROR: WAIT_LOG does not exist: ${WAIT_LOG}" >&2; exit 2; }
[[ -x "${LAUNCHER}" ]] || { echo "ERROR: launcher is not executable: ${LAUNCHER}" >&2; exit 2; }

echo "O2 singlet nadir pressure sequence"
echo "  waiting for PID = ${WAIT_PID}"
echo "  completion log  = ${WAIT_LOG}"
echo "  next pressures  = ${PSURF_SEQUENCE}"
echo "  CUDA device     = ${CUDA_DEVICE}"

while kill -0 "${WAIT_PID}" 2>/dev/null; do
  sleep "${POLL_SECONDS}"
done

if ! grep -Fq "Wrote NetCDF" "${WAIT_LOG}"; then
  echo "ERROR: preceding run exited without a successful NetCDF completion marker" >&2
  echo "       sequence stopped before starting ${PSURF_SEQUENCE}" >&2
  exit 1
fi

echo "Preceding run completed successfully. Starting queued pressure slices."
for psurf in ${PSURF_SEQUENCE}; do
  echo "Starting psurf=${psurf} hPa at $(date --iso-8601=seconds)"
  PSURF="${psurf}" \
  CUDA_DEVICE="${CUDA_DEVICE}" \
  OUTDIR="${OUTDIR}" \
    "${LAUNCHER}"
  echo "Completed psurf=${psurf} hPa at $(date --iso-8601=seconds)"
done

echo "O2 singlet nadir pressure sequence complete at $(date --iso-8601=seconds)"
