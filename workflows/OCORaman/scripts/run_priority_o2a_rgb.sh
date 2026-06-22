#!/usr/bin/env bash
set -euo pipefail

# Generic priority launcher for manuscript-critical O2 A-band Raman LUT scenes:
# all requested SZA values, but only the RGB-figure Lambertian albedos
# rho = 0.0, 0.3, 1.0. These are full-grid albedo indices 1, 7, 21.
#
# Each SZA is written to its own NetCDF chunk. This keeps failures local and
# lets the priority figure be assembled before the dense albedo LUT completes.
#
# Typical launch via psurf-specific thin wrappers:
#   CUDA_HOME=/usr/local/cuda CUDA_DEVICE=1 \
#   nohup workflows/OCORaman/scripts/run_priority_o2a_psurf1000_rgb.sh \
#     > /home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_priority_psurf1000_rgb/logs/nohup_driver.log 2>&1 &
#
#   CUDA_HOME=/usr/local/cuda CUDA_DEVICE=1 \
#   nohup workflows/OCORaman/scripts/run_priority_o2a_psurf750_rgb.sh \
#     > /home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_priority_psurf0750_rgb/logs/nohup_driver.log 2>&1 &
#
# Optional environment:
#   PSURF=1000 or 750
#   START_SZA_IDX=1
#   END_SZA_IDX=14
#   ALBEDO_IDXS=1,7,21
#   SKIP_COMPLETED=1
#   OUTDIR=/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_priority_psurfXXXX_rgb
#   JULIA_BIN=julia
#   JULIA_FLAGS="--pkgimages=no"
#   DRY_RUN=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${REPO_ROOT}"

PSURF="${PSURF:-1000}"
CUDA_DEVICE="${CUDA_DEVICE:-1}"
START_SZA_IDX="${START_SZA_IDX:-1}"
END_SZA_IDX="${END_SZA_IDX:-14}"
ALBEDO_IDXS="${ALBEDO_IDXS:-1,7,21}"
JULIA_BIN="${JULIA_BIN:-julia}"
JULIA_FLAGS="${JULIA_FLAGS:---pkgimages=no}"
DRY_RUN="${DRY_RUN:-0}"
SKIP_COMPLETED="${SKIP_COMPLETED:-1}"
# shellcheck disable=SC2206
JULIA_FLAGS_ARRAY=(${JULIA_FLAGS})

PSURF_TAG="$(printf '%04.0f' "${PSURF}")"
OUTDIR="${OUTDIR:-/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_priority_psurf${PSURF_TAG}_rgb}"
LOGDIR="${LOGDIR:-${OUTDIR}/logs}"

mkdir -p "${OUTDIR}" "${LOGDIR}"

if (( START_SZA_IDX < 1 || START_SZA_IDX > END_SZA_IDX || END_SZA_IDX > 14 )); then
  echo "SZA range ${START_SZA_IDX}:${END_SZA_IDX} is outside 1:14" >&2
  exit 2
fi

printf 'Priority O2A Raman LUT run\n'
printf '  repo: %s\n' "${REPO_ROOT}"
printf '  psurf: %s hPa\n' "${PSURF}"
printf '  CUDA_DEVICE: %s\n' "${CUDA_DEVICE}"
printf '  SZA indices: %s:%s\n' "${START_SZA_IDX}" "${END_SZA_IDX}"
printf '  albedo indices: %s  (rho = 0.0, 0.3, 1.0 for default grid)\n' "${ALBEDO_IDXS}"
printf '  output dir: %s\n' "${OUTDIR}"
printf '  log dir: %s\n' "${LOGDIR}"
printf '  Julia flags: %s\n' "${JULIA_FLAGS:-<none>}"
printf '  DRY_RUN: %s\n' "${DRY_RUN}"
printf '  skip completed: %s\n' "${SKIP_COMPLETED}"

albedo_tag() {
  printf '%s' "$1" | tr ',' '_' | tr ':' '-'
}

for sza_idx in $(seq "${START_SZA_IDX}" "${END_SZA_IDX}"); do
  requested_albedos="${ALBEDO_IDXS}"
  if [[ "${SKIP_COMPLETED}" != "0" ]]; then
    requested_albedos="$(
      env \
        PSURF="${PSURF}" \
        SZA_IDX="${sza_idx}" \
        ALBEDO_IDXS="${ALBEDO_IDXS}" \
        OUTDIR="${OUTDIR}" \
        "${JULIA_BIN}" "${JULIA_FLAGS_ARRAY[@]}" --project=. workflows/OCORaman/scripts/priority_missing_RamanLUT_O2A.jl
    )"
    requested_albedos="$(printf '%s' "${requested_albedos}" | tail -n 1)"
  fi

  if [[ -z "${requested_albedos}" ]]; then
    printf '\n[%s] Skipping psurf=%s sza_idx=%03d; requested priority albedos already complete\n' \
      "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${PSURF}" "${sza_idx}" | tee -a "${LOGDIR}/driver.log"
    continue
  fi

  chunk_albedo_tag="$(albedo_tag "${requested_albedos}")"
  chunk_tag="$(printf 'priority_psurf%s_sza%03d_alb%s' "${PSURF_TAG}" "${sza_idx}" "${chunk_albedo_tag}")"
  out_nc="${OUTDIR}/${chunk_tag}.nc"
  log_file="${LOGDIR}/${chunk_tag}.log"

  resume_flag=()
  if [[ -f "${out_nc}" ]]; then
    resume_flag=("RESUME=1")
  fi

  printf '\n[%s] Running priority chunk %s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${chunk_tag}" | tee -a "${LOGDIR}/driver.log"
  printf '  OUT_NC=%s\n' "${out_nc}" | tee -a "${LOGDIR}/driver.log"

  env \
    CUDA_DEVICE="${CUDA_DEVICE}" \
    PSURFS="${PSURF}" \
    SZA_IDXS="${sza_idx}" \
    ALBEDO_IDXS="${requested_albedos}" \
    CHUNK_TAG="${chunk_tag}" \
    OUTDIR="${OUTDIR}" \
    OUT_NC="${out_nc}" \
    DRY_RUN="${DRY_RUN}" \
    "${resume_flag[@]}" \
    "${JULIA_BIN}" "${JULIA_FLAGS_ARRAY[@]}" --project=. workflows/OCORaman/scripts/chunked_createRamanLUT_O2A.jl \
    2>&1 | tee -a "${log_file}"
done

printf '\n[%s] Finished priority O2A Raman LUT run for psurf=%s\n' \
  "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${PSURF}" | tee -a "${LOGDIR}/driver.log"
