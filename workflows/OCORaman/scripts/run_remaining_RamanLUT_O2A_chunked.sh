#!/usr/bin/env bash
set -euo pipefail

# Sequential chunk launcher for the remaining O2 A-band Raman LUT scenes.
#
# This is deliberately a thin wrapper around chunked_createRamanLUT_O2A.jl:
# each chunk writes a separate NetCDF file, and existing chunk files are
# resumed instead of overwritten.
#
# Required/typical environment:
#   PSURF=1000 or 750
#   CUDA_DEVICE=1
#   START_SZA_IDX=<1-based full-grid SZA index>
#   START_ALBEDO_IDX=<1-based full-grid albedo index for START_SZA_IDX>
#
# Optional:
#   END_SZA_IDX=14
#   CHUNK_ALBEDOS=4
#   SKIP_COMPLETED=1
#   OUTDIR=/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_chunked_psurf${PSURF}
#   LOGDIR=$OUTDIR/logs
#   JULIA_BIN=julia
#   JULIA_FLAGS="--pkgimages=no"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
cd "${REPO_ROOT}"

: "${PSURF:?Set PSURF, e.g. PSURF=1000 or PSURF=750}"
: "${START_SZA_IDX:?Set START_SZA_IDX, e.g. 4}"
: "${START_ALBEDO_IDX:?Set START_ALBEDO_IDX, e.g. 1}"

CUDA_DEVICE="${CUDA_DEVICE:-1}"
END_SZA_IDX="${END_SZA_IDX:-14}"
CHUNK_ALBEDOS="${CHUNK_ALBEDOS:-4}"
OUTDIR="${OUTDIR:-/home/sanghavi/data/RamanSIFgrid/o2a_raman_lut_chunked_psurf${PSURF}}"
LOGDIR="${LOGDIR:-${OUTDIR}/logs}"
JULIA_BIN="${JULIA_BIN:-julia}"
JULIA_FLAGS="${JULIA_FLAGS:---pkgimages=no}"
SKIP_COMPLETED="${SKIP_COMPLETED:-1}"
# shellcheck disable=SC2206
JULIA_FLAGS_ARRAY=(${JULIA_FLAGS})

mkdir -p "${OUTDIR}" "${LOGDIR}"

if (( START_SZA_IDX < 1 || START_SZA_IDX > END_SZA_IDX )); then
  echo "START_SZA_IDX=${START_SZA_IDX} is outside 1:${END_SZA_IDX}" >&2
  exit 2
fi
if (( START_ALBEDO_IDX < 1 || START_ALBEDO_IDX > 21 )); then
  echo "START_ALBEDO_IDX=${START_ALBEDO_IDX} is outside 1:21" >&2
  exit 2
fi
if (( CHUNK_ALBEDOS < 1 || CHUNK_ALBEDOS > 21 )); then
  echo "CHUNK_ALBEDOS=${CHUNK_ALBEDOS} is outside 1:21" >&2
  exit 2
fi

printf 'Chunked Raman LUT remainder run\n'
printf '  repo: %s\n' "${REPO_ROOT}"
printf '  psurf: %s hPa\n' "${PSURF}"
printf '  CUDA_DEVICE: %s\n' "${CUDA_DEVICE}"
printf '  SZA indices: %s:%s\n' "${START_SZA_IDX}" "${END_SZA_IDX}"
printf '  first-SZA albedo start index: %s\n' "${START_ALBEDO_IDX}"
printf '  albedos per chunk: %s\n' "${CHUNK_ALBEDOS}"
printf '  output dir: %s\n' "${OUTDIR}"
printf '  log dir: %s\n' "${LOGDIR}"
printf '  Julia flags: %s\n' "${JULIA_FLAGS:-<none>}"
printf '  skip completed: %s\n' "${SKIP_COMPLETED}"

albedo_tag() {
  printf '%s' "$1" | tr ',' '_' | tr ':' '-'
}

for sza_idx in $(seq "${START_SZA_IDX}" "${END_SZA_IDX}"); do
  alb_start=1
  if (( sza_idx == START_SZA_IDX )); then
    alb_start="${START_ALBEDO_IDX}"
  fi

  alb_idx="${alb_start}"
  while (( alb_idx <= 21 )); do
    alb_end=$(( alb_idx + CHUNK_ALBEDOS - 1 ))
    if (( alb_end > 21 )); then
      alb_end=21
    fi

    requested_albedos="${alb_idx}:${alb_end}"
    if [[ "${SKIP_COMPLETED}" != "0" ]]; then
      requested_albedos="$(
        env \
          PSURF="${PSURF}" \
          SZA_IDX="${sza_idx}" \
          ALBEDO_IDXS="${requested_albedos}" \
          OUTDIR="${OUTDIR}" \
          "${JULIA_BIN}" "${JULIA_FLAGS_ARRAY[@]}" --project=. workflows/OCORaman/scripts/priority_missing_RamanLUT_O2A.jl
      )"
      requested_albedos="$(printf '%s' "${requested_albedos}" | tail -n 1)"
    fi

    if [[ -z "${requested_albedos}" ]]; then
      printf '\n[%s] Skipping psurf=%s sza_idx=%03d alb=%03d-%03d; scenes already complete\n' \
        "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${PSURF}" "${sza_idx}" "${alb_idx}" "${alb_end}" | tee -a "${LOGDIR}/driver.log"
      alb_idx=$(( alb_end + 1 ))
      continue
    fi

    chunk_albedo_tag="$(albedo_tag "${requested_albedos}")"
    chunk_tag="$(printf 'remain_psurf%04.0f_sza%03d_alb%s' "${PSURF}" "${sza_idx}" "${chunk_albedo_tag}")"
    out_nc="${OUTDIR}/${chunk_tag}.nc"
    log_file="${LOGDIR}/${chunk_tag}.log"

    resume_flag=()
    if [[ -f "${out_nc}" ]]; then
      resume_flag=("RESUME=1")
    fi

    printf '\n[%s] Running chunk %s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${chunk_tag}" | tee -a "${LOGDIR}/driver.log"
    printf '  OUT_NC=%s\n' "${out_nc}" | tee -a "${LOGDIR}/driver.log"

    env \
      CUDA_DEVICE="${CUDA_DEVICE}" \
      PSURFS="${PSURF}" \
      SZA_IDXS="${sza_idx}" \
      ALBEDO_IDXS="${requested_albedos}" \
      CHUNK_TAG="${chunk_tag}" \
      OUTDIR="${OUTDIR}" \
      OUT_NC="${out_nc}" \
      "${resume_flag[@]}" \
      "${JULIA_BIN}" "${JULIA_FLAGS_ARRAY[@]}" --project=. workflows/OCORaman/scripts/chunked_createRamanLUT_O2A.jl \
      2>&1 | tee -a "${log_file}"

    alb_idx=$(( alb_end + 1 ))
  done
done

printf '\n[%s] Finished chunked Raman LUT remainder run for psurf=%s\n' \
  "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${PSURF}" | tee -a "${LOGDIR}/driver.log"
