#!/usr/bin/env bash

# Submit the required state-018 smoke task and the dependent 40-state %2
# production array. This script never generates or changes truth products.

set -euo pipefail
umask 077

: "${GATTACA_SLURM_ACCOUNT:?set the Gattaca Slurm account}"
: "${GATTACA_SLURM_QOS:?set the Gattaca Slurm QOS}"
: "${ROUND4_CODE_CHECKPOINT_SHA:?set the approved checkpoint SHA}"
: "${ROUND4_PRIOR_SHA256:?set the round-4 prior SHA-256}"
: "${ROUND4_PRIOR_SUMMARY_SHA256:?set the prior-summary SHA-256}"
: "${ROUND4_SOURCE_PRIOR_SHA256:?set the source-prior SHA-256}"
: "${RRS_REPO:?set the separate, detached round-4 source checkout}"
: "${BOTTOM_XCO2_CAMPAIGN_ROOT:?set the existing data-bearing checkout's bottom-layer campaign root}"
: "${FULL_COLUMN_TRUTH_ROOT:?set the existing data-bearing checkout's full-column truth root}"

repo_root="${RRS_REPO}"
private_root="${RRS_PRIVATE_ROOT:-${HOME}/RRS_XCO2_private}"
bottom_campaign="${BOTTOM_XCO2_CAMPAIGN_ROOT}"
full_truth_root="${FULL_COLUMN_TRUTH_ROOT}"
for directory in "${repo_root}" "${private_root}" "${bottom_campaign}" \
                 "${full_truth_root}"; do
    [[ -d "${directory}" ]] || {
        echo "missing required round-4 directory: ${directory}" >&2
        exit 1
    }
done
repo_root="$(cd "${repo_root}" && pwd -P)"
private_root="$(cd "${private_root}" && pwd -P)"
bottom_campaign="$(cd "${bottom_campaign}" && pwd -P)"
full_truth_root="$(cd "${full_truth_root}" && pwd -P)"
round4_git_root="$(git -C "${repo_root}" rev-parse --show-toplevel)"
round4_git_root="$(cd "${round4_git_root}" && pwd -P)"
legacy_git_root="$(git -C "${bottom_campaign}" rev-parse --show-toplevel)"
legacy_git_root="$(cd "${legacy_git_root}" && pwd -P)"
[[ "${repo_root}" == "${round4_git_root}" ]] || {
    echo "RRS_REPO must name the canonical round-4 Git top level" >&2
    exit 1
}
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
stokes_coefficient_path="${legacy_git_root}/RRS_XCO2/inversion/instrument/representative_stokes_coefficients.nc"
scene_components_path="${bottom_campaign}/truth/scene_components.dat"
sif_template_path="${legacy_git_root}/src/SIF_emission/sif-spectra.csv"
for required in "${bottom_campaign}/truth/true_states.dat" \
                "${full_truth_root}/true_states.dat" \
                "${stokes_coefficient_path}" "${scene_components_path}" \
                "${sif_template_path}"; do
    [[ -f "${required}" ]] || {
        echo "missing data-checkout input: ${required}" >&2
        exit 1
    }
done
if path_is_within "${private_root}" "${repo_root}" || \
   path_is_within "${repo_root}" "${private_root}" || \
   path_is_within "${private_root}" "${legacy_git_root}" || \
   path_is_within "${legacy_git_root}" "${private_root}"; then
    echo "RRS_PRIVATE_ROOT must be disjoint from both Git checkouts" >&2
    exit 1
fi
[[ "${ROUND4_CODE_CHECKPOINT_SHA}" =~ ^[0-9a-f]{40}$ ]] || {
    echo "ROUND4_CODE_CHECKPOINT_SHA must be 40 lowercase hexadecimal characters" >&2
    exit 2
}
for hash in "${ROUND4_PRIOR_SHA256}" "${ROUND4_PRIOR_SUMMARY_SHA256}" \
            "${ROUND4_SOURCE_PRIOR_SHA256}"; do
    [[ "${hash}" =~ ^[0-9a-f]{64}$ ]] || {
        echo "round-4 input hashes must be 64 lowercase hexadecimal characters" >&2
        exit 2
    }
done
[[ "$(git -C "${repo_root}" rev-parse HEAD)" == \
   "${ROUND4_CODE_CHECKPOINT_SHA}" ]] || {
    echo "round-4 checkout does not match the approved checkpoint" >&2
    exit 1
}
if git -C "${repo_root}" symbolic-ref -q HEAD >/dev/null; then
    echo "round-4 source checkout must be detached before submission" >&2
    exit 1
fi
[[ -z "$(git -C "${repo_root}" status --porcelain --untracked-files=all)" ]] || {
    echo "round-4 source checkout must be completely clean" >&2
    exit 1
}

campaign_id="bottom_layer_round4_known_sif759_sif_on_acos_mapped_tapered_vertical_correlation_v1"
result_root="${private_root}/results/${campaign_id}"
log_root="${result_root}/slurm"
launcher="${repo_root}/RRS_XCO2/inversion/gattaca_round4_known_sif759_retrievals.sbatch"
jobs_environment="${result_root}/retrieval_setup/retrieval_jobs.env"

[[ -f "${launcher}" ]] || {
    echo "missing round-4 launcher: ${launcher}" >&2
    exit 1
}
[[ -f "${repo_root}/Manifest.toml" ]] || {
    echo "missing ignored runtime manifest: ${repo_root}/Manifest.toml" >&2
    echo "copy the validated Manifest.toml from the legacy checkout before submission" >&2
    exit 1
}
[[ ! -e "${jobs_environment}" ]] || {
    echo "refusing duplicate submission; inspect ${jobs_environment}" >&2
    exit 1
}

validate_export_value() {
    local name="$1"
    local value="$2"
    [[ -n "${value}" && "${value}" != *','* && "${value}" != *$'\n'* ]] || {
        echo "${name} cannot be empty or contain commas/newlines for sbatch --export" >&2
        exit 2
    }
}
for pair in \
        "RRS_REPO=${repo_root}" \
        "RRS_PRIVATE_ROOT=${private_root}" \
        "BOTTOM_XCO2_CAMPAIGN_ROOT=${bottom_campaign}" \
        "FULL_COLUMN_TRUTH_ROOT=${full_truth_root}" \
        "RETRIEVAL_STOKES_COEFFICIENT_PATH=${stokes_coefficient_path}" \
        "RETRIEVAL_SCENE_COMPONENTS_PATH=${scene_components_path}" \
        "RRS_XCO2_SIF_TEMPLATE_PATH=${sif_template_path}"; do
    validate_export_value "${pair%%=*}" "${pair#*=}"
done

mkdir -p "${log_root}" "$(dirname "${jobs_environment}")"

path_exports="RRS_REPO=${repo_root},RRS_PRIVATE_ROOT=${private_root},BOTTOM_XCO2_CAMPAIGN_ROOT=${bottom_campaign},FULL_COLUMN_TRUTH_ROOT=${full_truth_root},RETRIEVAL_STOKES_COEFFICIENT_PATH=${stokes_coefficient_path},RETRIEVAL_SCENE_COMPONENTS_PATH=${scene_components_path},RRS_XCO2_SIF_TEMPLATE_PATH=${sif_template_path}"
exports="ALL,${path_exports},GATTACA_ROUND4_PHASE=smoke,ROUND4_CODE_CHECKPOINT_SHA=${ROUND4_CODE_CHECKPOINT_SHA},ROUND4_PRIOR_SHA256=${ROUND4_PRIOR_SHA256},ROUND4_PRIOR_SUMMARY_SHA256=${ROUND4_PRIOR_SUMMARY_SHA256},ROUND4_SOURCE_PRIOR_SHA256=${ROUND4_SOURCE_PRIOR_SHA256}"
smoke_job="$(sbatch --parsable \
    --account="${GATTACA_SLURM_ACCOUNT}" \
    --qos="${GATTACA_SLURM_QOS}" \
    --array=7 \
    --output="${log_root}/round4-sif759-%A_%a.out" \
    --error="${log_root}/round4-sif759-%A_%a.err" \
    --export="${exports}" \
    "${launcher}")"

exports="ALL,${path_exports},GATTACA_ROUND4_PHASE=production,ROUND4_CODE_CHECKPOINT_SHA=${ROUND4_CODE_CHECKPOINT_SHA},ROUND4_PRIOR_SHA256=${ROUND4_PRIOR_SHA256},ROUND4_PRIOR_SUMMARY_SHA256=${ROUND4_PRIOR_SUMMARY_SHA256},ROUND4_SOURCE_PRIOR_SHA256=${ROUND4_SOURCE_PRIOR_SHA256}"
production_job="$(sbatch --parsable \
    --account="${GATTACA_SLURM_ACCOUNT}" \
    --qos="${GATTACA_SLURM_QOS}" \
    --dependency="afterok:${smoke_job}" \
    --array=0-39%2 \
    --output="${log_root}/round4-sif759-%A_%a.out" \
    --error="${log_root}/round4-sif759-%A_%a.err" \
    --export="${exports}" \
    "${launcher}")"

temporary="${jobs_environment}.tmp.$$"
{
    printf 'SMOKE_JOB=%s\nPRODUCTION_JOB=%s\n' \
        "${smoke_job}" "${production_job}"
    printf 'LOG_ROOT=%s\nCAMPAIGN_ID=%s\n' \
        "${log_root}" "${campaign_id}"
    printf 'CHECKPOINT=%s\nPRIOR_SHA256=%s\n' \
        "${ROUND4_CODE_CHECKPOINT_SHA}" "${ROUND4_PRIOR_SHA256}"
    printf 'PRIOR_SUMMARY_SHA256=%s\nSOURCE_PRIOR_SHA256=%s\n' \
        "${ROUND4_PRIOR_SUMMARY_SHA256}" "${ROUND4_SOURCE_PRIOR_SHA256}"
    printf 'RRS_REPO=%q\nLEGACY_INPUT_REPO=%q\n' \
        "${repo_root}" "${legacy_git_root}"
    printf 'BOTTOM_XCO2_CAMPAIGN_ROOT=%q\nFULL_COLUMN_TRUTH_ROOT=%q\n' \
        "${bottom_campaign}" "${full_truth_root}"
    printf 'RETRIEVAL_STOKES_COEFFICIENT_PATH=%q\n' \
        "${stokes_coefficient_path}"
    printf 'RETRIEVAL_SCENE_COMPONENTS_PATH=%q\n' \
        "${scene_components_path}"
    printf 'RRS_XCO2_SIF_TEMPLATE_PATH=%q\n' "${sif_template_path}"
} > "${temporary}"
mv "${temporary}" "${jobs_environment}"

printf 'Submitted smoke job: %s\n' "${smoke_job}"
printf 'Submitted dependent production array: %s\n' "${production_job}"
printf 'Job record: %s\n' "${jobs_environment}"
