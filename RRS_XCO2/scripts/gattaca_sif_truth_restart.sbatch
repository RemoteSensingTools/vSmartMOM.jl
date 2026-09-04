#!/usr/bin/env bash

# Directly regenerate the corrected SIF-on O2 A-band truth on Gattaca2.
#
# This is intentionally a truth-only array, not a retrieval launcher.  The
# eight tasks own disjoint (aerosol selection, surface block) partitions, and
# the array throttle permits at most two simultaneous GPUs.  Each task calls
# run_sif_truth_restart_partition.sh, which resumes versioned chunk checkpoints
# with FORCE=0 and holds an atomic claim for its partition.

#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-7%2
#SBATCH --job-name=rrs_sif_truth_v2
# Account, QOS, and absolute log destinations are intentionally supplied by
# the private host runbook at submission time, not embedded in tracked source.

set -euo pipefail
umask 077

task_id="${SLURM_ARRAY_TASK_ID:-}"
[[ "${task_id}" =~ ^[0-7]$ ]] || {
    echo "SLURM_ARRAY_TASK_ID must be an integer from 0 through 7" >&2
    exit 2
}

aerosol_filters=(none none none none aerosol aerosol aerosol aerosol)
state_ranges=(1:16 17:32 33:48 49:64 1:16 17:32 33:48 49:64)
aerosol_filter="${aerosol_filters[task_id]}"
state_range="${state_ranges[task_id]}"

# This mode validates the complete array mapping without touching files,
# GPUs, or Slurm.  It is used by the source-level regression check documented
# in the Gattaca runbook; production must never set it.
if [[ "${GATTACA_TRUTH_PLAN_ONLY:-0}" == 1 ]]; then
    printf 'task=%s aerosol_filter=%s state_range=%s\n' \
        "${task_id}" "${aerosol_filter}" "${state_range}"
    exit 0
fi

: "${SLURM_JOB_ID:?this launcher must run inside a Slurm allocation}"
: "${CAMPAIGN_CHECKPOINT_SHA:?export the approved full 40-character source SHA at submission}"
[[ "${CAMPAIGN_CHECKPOINT_SHA}" =~ ^[0-9a-f]{40}$ ]] || {
    echo "CAMPAIGN_CHECKPOINT_SHA must be a lowercase 40-character Git SHA" >&2
    exit 2
}

repo_root="${RRS_REPO:-${HOME}/code/uni_vSmartMOM}"
private_root="${RRS_PRIVATE_ROOT:-${HOME}/RRS_XCO2_private}"
restart_root="${SIF_RESTART_ROOT:-${private_root}/results/.sif_v2_restart}"
julia_bin="${JULIA_BIN:-${HOME}/software/julia-1.12.5/bin/julia}"
depot_root="${RRS_DEPOT:-${private_root}/julia_depot_1.12.5}"
input_root="${RRS_INPUTS:-${private_root}/rrs_inputs}"
launcher="${repo_root}/RRS_XCO2/scripts/run_sif_truth_restart_partition.sh"

for directory in "${repo_root}" "${private_root}" "${restart_root}" \
                 "${depot_root}" "${input_root}"; do
    [[ -d "${directory}" ]] || {
        echo "missing required directory: ${directory}" >&2
        exit 1
    }
done
repo_root="$(cd "${repo_root}" && pwd -P)"
private_root="$(cd "${private_root}" && pwd -P)"
restart_root="$(cd "${restart_root}" && pwd -P)"

# Truth output must remain both outside Git and visibly marked as staging.
case "${restart_root}/" in
    "${private_root}/"*) ;;
    *) echo "SIF_RESTART_ROOT must be below RRS_PRIVATE_ROOT" >&2; exit 1 ;;
esac
case "${restart_root}/" in
    "${repo_root}/"*)
        echo "refusing to write truth products inside the Git checkout" >&2
        exit 1
        ;;
esac
[[ "$(basename "${restart_root}")" == .* ]] || {
    echo "SIF_RESTART_ROOT must have a hidden basename (for example .sif_v2_restart)" >&2
    exit 1
}

[[ -x "${julia_bin}" ]] || {
    echo "missing Julia 1.12.5 executable: ${julia_bin}" >&2
    exit 1
}
[[ "$("${julia_bin}" --version)" == "julia version 1.12.5" ]] || {
    echo "the direct truth restart requires exactly Julia 1.12.5" >&2
    exit 1
}
[[ -x "${launcher}" ]] || {
    echo "missing executable partition launcher: ${launcher}" >&2
    exit 1
}

source_sha="$(git -C "${repo_root}" rev-parse HEAD)"
[[ "${source_sha}" == "${CAMPAIGN_CHECKPOINT_SHA}" ]] || {
    echo "source checkout ${source_sha} does not match approved checkpoint ${CAMPAIGN_CHECKPOINT_SHA}" >&2
    exit 1
}
if git -C "${repo_root}" symbolic-ref -q HEAD >/dev/null; then
    echo "Gattaca production checkout must be detached at the approved SHA" >&2
    exit 1
fi
git -C "${repo_root}" diff --quiet -- || {
    echo "tracked source files have uncommitted changes" >&2
    exit 1
}
git -C "${repo_root}" diff --cached --quiet -- || {
    echo "tracked source files have staged changes" >&2
    exit 1
}

[[ -n "${CUDA_VISIBLE_DEVICES:-}" && "${CUDA_VISIBLE_DEVICES}" != *,* ]] || {
    echo "Slurm must expose exactly one scheduler-assigned CUDA device" >&2
    exit 1
}
assigned_gpu_memory="$({
    nvidia-smi -i "${CUDA_VISIBLE_DEVICES}" \
        --query-gpu=memory.total --format=csv,noheader,nounits
} | tr -d '[:space:]')"
[[ "${assigned_gpu_memory}" =~ ^[0-9]+$ ]] || {
    echo "could not determine memory of the scheduler-assigned GPU" >&2
    exit 1
}
(( assigned_gpu_memory >= 49152 )) || {
    echo "assigned GPU has only ${assigned_gpu_memory} MiB; at least 49152 MiB is required" >&2
    exit 1
}

required_files=(
    "${restart_root}/true_states.dat"
    "${private_root}/external_inputs.sha256"
    "${repo_root}/Manifest.toml"
    "${repo_root}/RRS_XCO2/surface_albedos/lambertian_legendre_inputs.dat"
    "${input_root}/absco_v5.2/o2_v52_v2.jld2"
    "${input_root}/absco_v5.2/wh2o_v52.jld2"
    "${input_root}/absco_v5.2/wco2_v52.jld2"
    "${input_root}/absco_v5.2/sh2o_v52.jld2"
    "${input_root}/absco_v5.2/sco2_v52.jld2"
    "${input_root}/hitran/H2O.jld2"
    "${input_root}/solar/solar.out"
    "${repo_root}/src/SIF_emission/sif-spectra.csv"
)
for required in "${required_files[@]}"; do
    [[ -f "${required}" ]] || {
        echo "missing required direct-truth input: ${required}" >&2
        exit 1
    }
done

# The transfer manifest is itself recorded below, but hashing only the
# manifest would not detect a post-transfer mutation of one of its inputs.
# Its paths are relative to private_root; verify every listed byte before any
# producer process is launched.
if ! (
    cd "${private_root}"
    sha256sum --check --quiet external_inputs.sha256
); then
    echo "one or more transferred external inputs failed SHA-256 validation" >&2
    exit 1
fi

# A transferred restart contains four SIF-on scenes for each array task.  The
# Julia producer performs the authoritative variable/metadata checks, while
# this early guard catches incomplete transfers before a GPU is initialized.
first_state="${state_range%%:*}"
if [[ "${aerosol_filter}" == none ]]; then
    first_scene=$((first_state + 4))
    scene_dir="${restart_root}/clear"
else
    first_scene=$((first_state + 12))
    scene_dir="${restart_root}/aerosol"
fi
for offset in 0 1 2 3; do
    printf -v scene '%s/hiressim_%03d.nc' "${scene_dir}" "$((first_scene + offset))"
    [[ -f "${scene}" ]] || {
        echo "missing seeded SIF-on scene: ${scene}" >&2
        exit 1
    }
done

mkdir -p "${restart_root}/logs"
provenance="${restart_root}/logs/gattaca_${aerosol_filter}_${state_range//:/_}_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${task_id}.provenance.txt"
{
    printf 'started_utc=%s\n' "$(date -u +%FT%TZ)"
    printf 'source_sha=%s\n' "${source_sha}"
    printf 'host=%s\n' "$(hostname -f)"
    printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID}"
    printf 'slurm_array_job_id=%s\n' "${SLURM_ARRAY_JOB_ID:-none}"
    printf 'slurm_array_task_id=%s\n' "${task_id}"
    printf 'aerosol_filter=%s\nstate_range=%s\nforce=0\n' \
        "${aerosol_filter}" "${state_range}"
    printf 'restart_root=%s\nCUDA_VISIBLE_DEVICES=%s\n' \
        "${restart_root}" "${CUDA_VISIBLE_DEVICES}"
    "${julia_bin}" --version
    sha256sum "${private_root}/external_inputs.sha256" \
              "${repo_root}/Manifest.toml" \
              "${repo_root}/RRS_XCO2/surface_albedos/lambertian_legendre_inputs.dat" \
              "${restart_root}/true_states.dat" \
              "${repo_root}/src/SIF_emission/sif-spectra.csv"
    nvidia-smi -i "${CUDA_VISIBLE_DEVICES}" \
        --query-gpu=name,memory.total,driver_version --format=csv,noheader
} > "${provenance}"

export REPO_ROOT="${repo_root}"
export SIF_RESTART_ROOT="${restart_root}"
export JULIA_BIN="${julia_bin}"
export JULIA_DEPOT_PATH="${depot_root}"
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export JULIA_CPU_TARGET="${JULIA_CPU_TARGET:-native}"
export VSMARTMOM_ABSCO_DIR="${input_root}/absco_v5.2"
export O2_H2O_LUT="${input_root}/hitran/H2O.jld2"
export SOLAR_OUT="${input_root}/solar/solar.out"

"${launcher}" "${aerosol_filter}" "${state_range}"
printf 'completed_utc=%s\n' "$(date -u +%FT%TZ)" >> "${provenance}"
