#!/usr/bin/env bash

# Run one disjoint direct-RT partition of the corrected full-column SIF-on
# O2 truth restart. Prepare the staging tree first with
# prepare_sif_truth_restart.jl. This launcher always resumes checkpoints and
# never discards completed chunks.

set -euo pipefail
umask 077

usage() {
    echo "usage: $0 {none|aerosol} {1:16|17:32|33:48|49:64} [physical_gpu]" >&2
}

aerosol_filter="${1:-}"
state_range="${2:-}"
physical_gpu="${3:-}"
[[ "${aerosol_filter}" == none || "${aerosol_filter}" == aerosol ]] || {
    usage
    exit 2
}
case "${state_range}" in
    1:16|17:32|33:48|49:64) ;;
    *) usage; exit 2 ;;
esac

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
default_repo_root="$(cd "${script_dir}/../.." && pwd)"
repo_root="${REPO_ROOT:-${default_repo_root}}"
truth_root="${FULL_COLUMN_TRUTH_ROOT:-${repo_root}/RRS_XCO2/truth_map}"
restart_root="${SIF_RESTART_ROOT:-${truth_root}/.sif_v2_restart}"
state_table="${restart_root}/true_states.dat"
output_root="${restart_root}/$([[ "${aerosol_filter}" == aerosol ]] && echo aerosol || echo clear)"
program="${repo_root}/RRS_XCO2/scripts/regenerate_o2_preserve_co2.jl"
julia_bin="${JULIA_BIN:-julia}"
first_state="${state_range%%:*}"
last_state="${state_range##*:}"
chunk_points=2735
[[ "${aerosol_filter}" == aerosol ]] && chunk_points=64

for required in "${state_table}" "${program}"; do
    [[ -f "${required}" ]] || {
        echo "missing required restart input: ${required}" >&2
        exit 1
    }
done

if [[ "${julia_bin}" == */* ]]; then
    [[ -x "${julia_bin}" ]] || {
        echo "Julia executable is missing or not executable: ${julia_bin}" >&2
        exit 1
    }
else
    command -v "${julia_bin}" >/dev/null 2>&1 || {
        echo "Julia executable is not on PATH: ${julia_bin}" >&2
        exit 1
    }
fi

# Slurm has already isolated the allocated GPU and exposes it as logical zero.
# On Curry/Wurst, require an explicit physical device and mask to it.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    [[ -n "${CUDA_VISIBLE_DEVICES:-}" ]] || {
        echo "Slurm did not expose an allocated CUDA device" >&2
        exit 1
    }
    logical_gpu=0
else
    [[ "${physical_gpu}" =~ ^[0-9]+$ ]] || {
        echo "physical_gpu is required outside Slurm" >&2
        usage
        exit 2
    }
    gpu_rows="$(nvidia-smi -i "${physical_gpu}" \
        --query-compute-apps=pid,process_name --format=csv,noheader,nounits \
        2>/dev/null || true)"
    [[ -z "${gpu_rows//[[:space:]]/}" ]] || {
        echo "physical GPU ${physical_gpu} is busy; refusing overlap" >&2
        printf '%s\n' "${gpu_rows}" >&2
        exit 1
    }
    export CUDA_VISIBLE_DEVICES="${physical_gpu}"
    logical_gpu=0
fi

mkdir -p "${restart_root}/logs" "${restart_root}/.claims"
claim_dir="${restart_root}/.claims/${aerosol_filter}_${first_state}_${last_state}.claim"
if ! mkdir "${claim_dir}"; then
    echo "restart partition is already claimed: ${claim_dir}" >&2
    [[ -f "${claim_dir}/owner.txt" ]] && sed -n '1,30p' "${claim_dir}/owner.txt" >&2
    exit 1
fi

cleanup_claim() {
    rm -f "${claim_dir}/owner.txt"
    rmdir "${claim_dir}" 2>/dev/null || true
}
trap cleanup_claim EXIT INT TERM

{
    echo "host=$(hostname -f)"
    echo "pid=$$"
    echo "started_utc=$(date -u +%FT%TZ)"
    echo "slurm_job_id=${SLURM_JOB_ID:-none}"
    echo "slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-none}"
    echo "aerosol_filter=${aerosol_filter}"
    echo "state_range=${state_range}"
    echo "force=0"
} > "${claim_dir}/owner.txt"

log_file="${restart_root}/logs/$(hostname -s)_${aerosol_filter}_${first_state}_${last_state}_$(date -u +%Y%m%dT%H%M%SZ).log"

export CUDA_DEVICE="${logical_gpu}"
export TRUTH_OUT="${output_root}"
export TRUTH_STATE_FILE="${state_table}"
export FIRST_STATE="${first_state}"
export LAST_STATE="${last_state}"
export AEROSOL_CASE_FILTER="${aerosol_filter}"
export SIF_CASE_FILTER=on
export O2_CHUNK_POINTS="${chunk_points}"
export FORCE=0
export RRS_XCO2_FLOAT_TYPE=Float32
export NLAYERS=16
export AEROSOL_NSTREAMS=9
export TRUTH_SZA_DEG=30
export TRUTH_VZA_DEG=0
export TRUTH_RELATIVE_AZIMUTH_DEG=0
export RAMAN_SHOULDER_CM=234

if [[ -z "${JULIA_DEPOT_PATH:-}" ]]; then
    export JULIA_DEPOT_PATH="${repo_root}/RRS_XCO2/.julia_depot_$(hostname -s):${HOME}/.julia"
fi

cd "${repo_root}"
{
    echo "$(date -u +%FT%TZ) corrected_SIF_truth_start filter=${aerosol_filter} states=${state_range} chunks=${chunk_points}"
    echo "host=$(hostname -f) slurm_job=${SLURM_JOB_ID:-none} CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
    "${julia_bin}" --project="${repo_root}" --startup-file=no "${program}"
    echo "$(date -u +%FT%TZ) corrected_SIF_truth_complete filter=${aerosol_filter} states=${state_range}"
} 2>&1 | tee "${log_file}"
