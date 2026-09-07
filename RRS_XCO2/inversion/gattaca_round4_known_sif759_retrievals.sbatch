#!/usr/bin/env bash

# One Gattaca2 GPU task owns one of the 40 bottom-layer SIF-on states. Submit
# state 018/task 7 as a smoke test first; production is additionally protected
# by both an afterok dependency and a validated smoke receipt.

#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-39%2
#SBATCH --job-name=rrs_r4_sif759
# Account, QOS, and absolute output/error paths are supplied at submission.

set -euo pipefail
umask 077

SIF_STATES=(
    006 007 008 009 010
    016 017 018 019 020
    026 027 028 029 030
    036 037 038 039 040
    046 047 048 049 050
    056 057 058 059 060
    066 067 068 069 070
    076 077 078 079 080
)

task_id="${SLURM_ARRAY_TASK_ID:-}"
[[ "${task_id}" =~ ^[0-9]+$ ]] &&
    (( task_id >= 0 && task_id < ${#SIF_STATES[@]} )) || {
        echo "SLURM_ARRAY_TASK_ID must be an integer from 0 through 39" >&2
        exit 2
    }
state="${SIF_STATES[task_id]}"

campaign_id="bottom_layer_round4_known_sif759_sif_on_acos_mapped_tapered_vertical_correlation_v1"
prior_filename="apriori_states_round4_known_sif759_on_acos_mapped_tapered_vertical_correlation.nc"
summary_filename="${prior_filename%.nc}.dat"

if [[ "${GATTACA_ROUND4_PLAN_ONLY:-0}" == 1 ]]; then
    printf 'task=%s state=%s sif_filter=on active=29 campaign=%s\n' \
        "${task_id}" "${state}" "${campaign_id}"
    exit 0
fi

: "${SLURM_JOB_ID:?this launcher must run inside a Slurm allocation}"
: "${ROUND4_CODE_CHECKPOINT_SHA:?export the approved 40-character checkpoint}"
: "${ROUND4_PRIOR_SHA256:?export the approved round-4 prior SHA-256}"
: "${ROUND4_PRIOR_SUMMARY_SHA256:?export the approved prior-summary SHA-256}"
: "${ROUND4_SOURCE_PRIOR_SHA256:?export the approved source-prior SHA-256}"
: "${RRS_REPO:?export the separate round-4 source checkout}"
: "${BOTTOM_XCO2_CAMPAIGN_ROOT:?export the bottom-layer campaign root in the existing data-bearing checkout}"
: "${FULL_COLUMN_TRUTH_ROOT:?export the full-column truth root in the existing data-bearing checkout}"
[[ "${ROUND4_CODE_CHECKPOINT_SHA}" =~ ^[0-9a-f]{40}$ ]] || exit 2
for hash in "${ROUND4_PRIOR_SHA256}" "${ROUND4_PRIOR_SUMMARY_SHA256}" \
            "${ROUND4_SOURCE_PRIOR_SHA256}"; do
    [[ "${hash}" =~ ^[0-9a-f]{64}$ ]] || exit 2
done

phase="${GATTACA_ROUND4_PHASE:-}"
[[ "${phase}" == smoke || "${phase}" == production ]] || {
    echo "GATTACA_ROUND4_PHASE must be smoke or production" >&2
    exit 2
}
if [[ "${phase}" == smoke && "${task_id}" != 7 ]]; then
    echo "round-4 smoke requires array task 7 (state 018)" >&2
    exit 2
fi

repo_root="${RRS_REPO}"
private_root="${RRS_PRIVATE_ROOT:-${HOME}/RRS_XCO2_private}"
input_root="${RRS_INPUTS:-${private_root}/rrs_inputs}"
depot_root="${RRS_DEPOT:-${private_root}/julia_depot_1.12.5}"
restart_root="${SIF_RESTART_ROOT:-${private_root}/results/.sif_v2_restart}"
julia_bin="${JULIA_BIN:-${HOME}/software/julia-1.12.5/bin/julia}"
bottom_campaign="${BOTTOM_XCO2_CAMPAIGN_ROOT}"
full_truth_root="${FULL_COLUMN_TRUTH_ROOT}"
bottom_archive_root="${private_root}/archive/sif_wavelength_integral_0p5_20260904/bottom_layer_XCO2_retrievals"
for directory in "${repo_root}" "${private_root}" "${input_root}" \
                 "${depot_root}" "${restart_root}" "${bottom_campaign}" \
                 "${full_truth_root}" "${bottom_archive_root}"; do
    [[ -d "${directory}" ]] || {
        echo "missing round-4 directory: ${directory}" >&2
        exit 1
    }
done

repo_root="$(cd "${repo_root}" && pwd -P)"
private_root="$(cd "${private_root}" && pwd -P)"
input_root="$(cd "${input_root}" && pwd -P)"
depot_root="$(cd "${depot_root}" && pwd -P)"
restart_root="$(cd "${restart_root}" && pwd -P)"
bottom_campaign="$(cd "${bottom_campaign}" && pwd -P)"
full_truth_root="$(cd "${full_truth_root}" && pwd -P)"
bottom_archive_root="$(cd "${bottom_archive_root}" && pwd -P)"

round4_git_root="$(git -C "${repo_root}" rev-parse --show-toplevel)"
round4_git_root="$(cd "${round4_git_root}" && pwd -P)"
[[ "${repo_root}" == "${round4_git_root}" ]] || {
    echo "RRS_REPO must name the canonical round-4 Git top level" >&2
    exit 1
}
legacy_git_root="$(git -C "${bottom_campaign}" rev-parse --show-toplevel)"
legacy_git_root="$(cd "${legacy_git_root}" && pwd -P)"
path_is_within() {
    [[ "$1" == "$2" || "$1" == "$2/"* ]]
}
if path_is_within "${repo_root}" "${legacy_git_root}" || \
   path_is_within "${legacy_git_root}" "${repo_root}"; then
    echo "round 4 requires separate, non-nested code and input checkouts" >&2
    exit 1
fi
[[ "${bottom_campaign}" == "${legacy_git_root}/RRS_XCO2/bottom_layer_XCO2_retrievals" ]] || {
    echo "BOTTOM_XCO2_CAMPAIGN_ROOT is not the canonical legacy campaign path" >&2
    exit 1
}
[[ "${full_truth_root}" == "${legacy_git_root}/RRS_XCO2/truth_map" ]] || {
    echo "FULL_COLUMN_TRUTH_ROOT must come from the same data-bearing checkout" >&2
    exit 1
}
if path_is_within "${private_root}" "${repo_root}" || \
   path_is_within "${repo_root}" "${private_root}" || \
   path_is_within "${private_root}" "${legacy_git_root}" || \
   path_is_within "${legacy_git_root}" "${private_root}"; then
    echo "RRS_PRIVATE_ROOT must be disjoint from both Git checkouts" >&2
    exit 1
fi

truth_root="${bottom_campaign}/truth"
measurement_root="${truth_root}/OCO_radiances"
noise_root="${measurement_root}/noise_covariances"
result_root="${private_root}/results/${campaign_id}"
output_root="${result_root}/retrievals"
prior_path="${result_root}/retrieval_setup/${prior_filename}"
summary_path="${result_root}/retrieval_setup/${summary_filename}"
source_prior_path="${result_root}/retrieval_setup/source_apriori_states_acos_mapped_tapered_vertical_correlation.nc"
stokes_coefficient_path="${legacy_git_root}/RRS_XCO2/inversion/instrument/representative_stokes_coefficients.nc"
scene_components_path="${bottom_campaign}/truth/scene_components.dat"
sif_template_path="${legacy_git_root}/src/SIF_emission/sif-spectra.csv"
readiness="${repo_root}/RRS_XCO2/inversion/GattacaRound4SIFReadiness.jl"
runner="${repo_root}/RRS_XCO2/inversion/run_round4_known_sif_retrievals.jl"

for required in "${readiness}" "${runner}" "${repo_root}/Manifest.toml" \
                "${truth_root}/true_states.dat" \
                "${full_truth_root}/true_states.dat" \
                "${prior_path}" "${summary_path}" "${source_prior_path}" \
                "${stokes_coefficient_path}" "${scene_components_path}" \
                "${sif_template_path}" \
                "${private_root}/external_inputs.sha256"; do
    [[ -f "${required}" ]] || {
        echo "missing round-4 input: ${required}" >&2
        exit 1
    }
done
for directory in "${measurement_root}" "${noise_root}"; do
    [[ -d "${directory}" ]] || {
        echo "missing round-4 input directory: ${directory}" >&2
        exit 1
    }
done

expected_output="${private_root}/results/${campaign_id}/retrievals"
[[ "${output_root}" == "${expected_output}" ]] || exit 1
case "${output_root}/" in
    "${repo_root}/"*|"${legacy_git_root}/"*)
        echo "refusing to write round-4 output inside either Git checkout" >&2
        exit 1
        ;;
esac

[[ -x "${julia_bin}" ]] || {
    echo "missing Julia executable: ${julia_bin}" >&2
    exit 1
}
[[ "$("${julia_bin}" --version)" == "julia version 1.12.5" ]] || {
    echo "round 4 requires Julia 1.12.5" >&2
    exit 1
}
source_sha="$(git -C "${repo_root}" rev-parse HEAD)"
[[ "${source_sha}" == "${ROUND4_CODE_CHECKPOINT_SHA}" ]] || {
    echo "checkout does not match the approved checkpoint" >&2
    exit 1
}
if git -C "${repo_root}" symbolic-ref -q HEAD >/dev/null; then
    echo "Gattaca production checkout must be detached" >&2
    exit 1
fi
git -C "${repo_root}" diff --quiet -- || {
    echo "tracked checkout has unstaged changes" >&2
    exit 1
}
git -C "${repo_root}" diff --cached --quiet -- || {
    echo "tracked checkout has staged changes" >&2
    exit 1
}
[[ -z "$(git -C "${repo_root}" status --porcelain --untracked-files=all)" ]] || {
    echo "round-4 source checkout must be completely clean" >&2
    exit 1
}

[[ -n "${CUDA_VISIBLE_DEVICES:-}" && "${CUDA_VISIBLE_DEVICES}" != *,* ]] || {
    echo "Slurm must expose exactly one CUDA device" >&2
    exit 1
}
gpu_memory="$(nvidia-smi -i "${CUDA_VISIBLE_DEVICES}" \
    --query-gpu=memory.total --format=csv,noheader,nounits | tr -d '[:space:]')"
[[ "${gpu_memory}" =~ ^[0-9]+$ ]] && (( gpu_memory >= 49152 )) || {
    echo "the assigned GPU has insufficient memory" >&2
    exit 1
}

(
    cd "${private_root}"
    sha256sum --check --quiet external_inputs.sha256
)

mkdir -p "${result_root}/logs" "${result_root}/state_claims" \
         "${output_root}"

export RRS_REPO="${repo_root}"
export RRS_PRIVATE_ROOT="${private_root}"
export RRS_INPUTS="${input_root}"
export RRS_DEPOT="${depot_root}"
export SIF_RESTART_ROOT="${restart_root}"
export BOTTOM_XCO2_CAMPAIGN_ROOT="${bottom_campaign}"
export BOTTOM_SIF_V1_ARCHIVE_ROOT="${bottom_archive_root}"
export FULL_COLUMN_TRUTH_ROOT="${full_truth_root}"
export RETRIEVAL_TRUTH_TABLE="${truth_root}/true_states.dat"
export RETRIEVAL_MEASUREMENT_DIR="${measurement_root}"
export RETRIEVAL_NOISE_DIR="${noise_root}"
export ROUND4_SOURCE_PRIOR_PATH="${source_prior_path}"
export RETRIEVAL_PRIOR_PATH="${prior_path}"
export ROUND4_PRIOR_SUMMARY_PATH="${summary_path}"
export RETRIEVAL_OUTPUT_ROOT="${output_root}"
export RETRIEVAL_STOKES_COEFFICIENT_PATH="${stokes_coefficient_path}"
export RETRIEVAL_SCENE_COMPONENTS_PATH="${scene_components_path}"
export RRS_XCO2_SIF_TEMPLATE_PATH="${sif_template_path}"
export JULIA_DEPOT_PATH="${depot_root}"
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export JULIA_CPU_TARGET="${JULIA_CPU_TARGET:-native}"
export VSMARTMOM_ABSCO_DIR="${input_root}/absco_v5.2"
export O2_H2O_LUT="${input_root}/hitran/H2O.jld2"
export SOLAR_OUT="${input_root}/solar/solar.out"

readiness_log="${result_root}/logs/readiness_job${SLURM_JOB_ID}_task${task_id}.log"
GATTACA_ROUND4_ACTION=prepare \
    "${julia_bin}" --project="${repo_root}" --startup-file=no \
    "${readiness}" > "${readiness_log}"

identity_env="${result_root}/retrieval_setup/campaign_identity.env"
[[ -f "${identity_env}" ]] || {
    echo "readiness did not publish campaign_identity.env" >&2
    exit 1
}
# The file is generated solely from validated hexadecimal values.
set -a
source "${identity_env}"
set +a
identity="${output_root}/.control/campaign_identity.dat"
identity_sha="$(sha256sum "${identity}" | awk '{print $1}')"
smoke_receipt="${result_root}/.control/state018_perturbation11_smoke.complete"
if [[ "${phase}" == production ]]; then
    [[ -f "${smoke_receipt}" ]] || {
        echo "production is blocked until the state-018 smoke passes" >&2
        exit 1
    }
    grep -qxF "campaign_identity_sha256=${identity_sha}" \
        "${smoke_receipt}" || {
        echo "smoke receipt belongs to a different campaign identity" >&2
        exit 1
    }
fi

claim="${result_root}/state_claims/state${state}.claim"
if ! mkdir "${claim}" 2>/dev/null; then
    echo "state ${state} is already claimed" >&2
    exit 1
fi
cleanup_claim() {
    rm -f "${claim}/owner.txt"
    rmdir "${claim}" 2>/dev/null || true
}
trap cleanup_claim EXIT INT TERM
{
    printf 'host=%s\npid=%s\nstate=%s\nphase=%s\n' \
        "$(hostname -f)" "$$" "${state}" "${phase}"
    printf 'slurm_job_id=%s\ncheckpoint=%s\nprior_sha256=%s\n' \
        "${SLURM_JOB_ID}" "${source_sha}" "${ROUND4_PRIOR_SHA256}"
    printf 'campaign_identity_sha256=%s\n' \
        "${ROUND4_CAMPAIGN_IDENTITY_SHA256}"
} > "${claim}/owner.txt"

first_perturbation=1
last_perturbation=11
write_manifest=0
if [[ "${phase}" == smoke ]]; then
    first_perturbation=11
    write_manifest=1
fi

env \
    CUDA_DEVICE=0 \
    ROUND4_SIF_MODE=on \
    ROUND4_KNOWN_WAVELENGTH_NM=759 \
    SIF_CASE_FILTER=on \
    AEROSOL_CASE_FILTER=all \
    FIRST_STATE="${state#0}" \
    LAST_STATE="${state#0}" \
    FIRST_PERTURBATION="${first_perturbation}" \
    LAST_PERTURBATION="${last_perturbation}" \
    RETRIEVAL_CLASS=paired \
    RETRIEVAL_ARCH=GPU \
    RETRIEVAL_FLOAT_TYPE=Float32 \
    RETRIEVAL_WRITE_MANIFEST="${write_manifest}" \
    FORCE=0 \
    FAIL_FAST=1 \
    "${julia_bin}" --project="${repo_root}" --startup-file=no "${runner}"

GATTACA_ROUND4_ACTION=validate-output \
GATTACA_ROUND4_PHASE="${phase}" \
EXPECTED_STATE="${state#0}" \
    "${julia_bin}" --project="${repo_root}" --startup-file=no "${readiness}"

if [[ "${phase}" == smoke ]]; then
    mkdir -p "$(dirname "${smoke_receipt}")"
    temporary="${smoke_receipt}.tmp.$$"
    {
        printf 'campaign_id=%s\nstate=018\nperturbation=11\n' "${campaign_id}"
        printf 'source_sha=%s\nprior_sha256=%s\n' \
            "${source_sha}" "${ROUND4_PRIOR_SHA256}"
        printf 'campaign_identity_sha256=%s\nvalidated_utc=%s\n' \
            "${identity_sha}" "$(date -u +%FT%TZ)"
    } > "${temporary}"
    mv "${temporary}" "${smoke_receipt}"
fi
