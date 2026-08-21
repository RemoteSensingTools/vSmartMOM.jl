#!/usr/bin/env bash
set -euo pipefail

# Periodically promote validated chunk/subset Raman LUT scenes into the
# psurf-level NetCDF containers.  The Julia promoter skips unreadable/open
# chunks and only copies complete finite scenes.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${REPO_ROOT}"

DATA_ROOT="${DATA_ROOT:-/home/sanghavi/data/RamanSIFgrid/O2ABand}"
PSURFS="${PSURFS:-1000,750}"
INTERVAL_SECONDS="${INTERVAL_SECONDS:-3600}"
MIN_AGE_SECONDS="${MIN_AGE_SECONDS:-60}"
JULIA_BIN="${JULIA_BIN:-julia}"
JULIA_FLAGS="${JULIA_FLAGS:---pkgimages=no}"
DRY_RUN="${DRY_RUN:-0}"
OVERWRITE="${OVERWRITE:-0}"
COPY_OPTICAL="${COPY_OPTICAL:-1}"
LOGDIR="${LOGDIR:-${DATA_ROOT}/logs}"
PROMOTE_LOG="${PROMOTE_LOG:-${LOGDIR}/hourly_promote_RamanLUT_chunks.log}"
LOCKFILE="${LOCKFILE:-${LOGDIR}/hourly_promote_RamanLUT_chunks.lock}"
PIDFILE="${PIDFILE:-${LOGDIR}/hourly_promote_RamanLUT_chunks.pid}"

mkdir -p "${LOGDIR}"
echo "$$" > "${PIDFILE}"

# shellcheck disable=SC2206
JULIA_FLAGS_ARRAY=(${JULIA_FLAGS})

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')" "$*" | tee -a "${PROMOTE_LOG}"
}

run_once() {
  flock -n "${LOCKFILE}" \
    env \
      DATA_ROOT="${DATA_ROOT}" \
      PSURFS="${PSURFS}" \
      MIN_AGE_SECONDS="${MIN_AGE_SECONDS}" \
      DRY_RUN="${DRY_RUN}" \
      OVERWRITE="${OVERWRITE}" \
      COPY_OPTICAL="${COPY_OPTICAL}" \
      "${JULIA_BIN}" "${JULIA_FLAGS_ARRAY[@]}" --project=. \
        workflows/OCORaman/scripts/promote_RamanLUT_chunks.jl
}

log "Starting hourly Raman LUT chunk promoter"
log "repo=${REPO_ROOT}"
log "DATA_ROOT=${DATA_ROOT}"
log "PSURFS=${PSURFS}"
log "INTERVAL_SECONDS=${INTERVAL_SECONDS}"
log "MIN_AGE_SECONDS=${MIN_AGE_SECONDS}"
log "DRY_RUN=${DRY_RUN}"
log "OVERWRITE=${OVERWRITE}"
log "COPY_OPTICAL=${COPY_OPTICAL}"

while true; do
  log "Promotion pass starting"
  if run_once >> "${PROMOTE_LOG}" 2>&1; then
    log "Promotion pass finished"
  else
    rc=$?
    log "Promotion pass failed or skipped by lock; rc=${rc}"
  fi
  sleep "${INTERVAL_SECONDS}"
done
