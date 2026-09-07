#!/usr/bin/env bash

# Dedicated local launcher for the round-4, known-SIF-at-759-nm, SIF-off
# retrieval campaign.  The three accepted worker names have immutable,
# disjoint state assignments. Curry GPU 0 is deliberately not an option.
# This script is fail-closed and only validates unless
# ROUND4_NOSIF_EXECUTE=1 is supplied with caller-approved input/code hashes.

set -euo pipefail

WORKER="${1:-}"
[[ $# -le 1 ]] || {
    echo "round-4 assignments are fixed; do not pass a state range" >&2
    exit 2
}

case "${WORKER}" in
    curry0)
        echo "curry0 is reserved and is rejected by the round-4 launcher" >&2
        exit 2
        ;;
    curry1)
        EXPECTED_HOST=curry
        PHYSICAL_GPU=1
        SCENE_CLASS=none
        DEPOT_NAME=.julia_depot_curry
        STATE_BLOCKS=("1:5" "21:25" "41:45" "61:65")
        ;;
    wurst0)
        EXPECTED_HOST=wurst
        PHYSICAL_GPU=0
        SCENE_CLASS=aerosol
        DEPOT_NAME=.julia_depot_wurst
        STATE_BLOCKS=("11:15" "51:55")
        ;;
    wurst1)
        EXPECTED_HOST=wurst
        PHYSICAL_GPU=1
        SCENE_CLASS=aerosol
        DEPOT_NAME=.julia_depot_wurst
        STATE_BLOCKS=("31:35" "71:75")
        ;;
    *)
        echo "usage: $0 curry1|wurst0|wurst1" >&2
        echo "curry0 is intentionally reserved" >&2
        exit 2
        ;;
esac

SCRIPT_DIRECTORY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="${REPO_ROOT:-$(git -C "${SCRIPT_DIRECTORY}" rev-parse --show-toplevel)}"
CAMPAIGN_ROOT="${BOTTOM_RETRIEVAL_CAMPAIGN_ROOT:-${REPO_ROOT}/RRS_XCO2/bottom_layer_XCO2_retrievals}"
ROUND4_ROOT="${CAMPAIGN_ROOT}/round4_known_sif759"
TRUTH_ROOT="${CAMPAIGN_ROOT}/truth"
OCO_ROOT="${TRUTH_ROOT}/OCO_radiances"
NOISE_ROOT="${OCO_ROOT}/noise_covariances"
MODEL_NAME="acos_mapped_tapered_vertical_correlation"
PRIOR_PATH="${ROUND4_NOSIF_PRIOR_PATH:-${ROUND4_ROOT}/retrieval_setup/apriori_states_round4_known_sif759_off_${MODEL_NAME}.nc}"
SOURCE_PRIOR_PATH="${ROUND4_SOURCE_PRIOR_PATH:-${CAMPAIGN_ROOT}/retrieval_setup/apriori_states_${MODEL_NAME}.nc}"
OUTPUT_ROOT="${ROUND4_NOSIF_OUTPUT_ROOT:-${ROUND4_ROOT}/retrievals_nosif}"
RUNNER="${ROUND4_RUNNER_PATH:-${REPO_ROOT}/RRS_XCO2/inversion/run_round4_known_sif_retrievals.jl}"
PREFLIGHT="${ROUND4_PREFLIGHT_PATH:-${REPO_ROOT}/RRS_XCO2/inversion/preflight_round4_known_sif759_no_sif_retrievals.jl}"
JULIA_BIN="${JULIA_BIN:-julia}"
JULIA_DEPOT_PATH_VALUE="${JULIA_DEPOT_PATH_VALUE:-${REPO_ROOT}/RRS_XCO2/${DEPOT_NAME}:${HOME}/.julia}"
EXECUTE="${ROUND4_NOSIF_EXECUTE:-0}"
PLAN_ONLY="${ROUND4_NOSIF_PLAN_ONLY:-0}"
PRINT_HASHES="${ROUND4_NOSIF_PRINT_HASHES:-0}"
REQUIRE_ROUND3_COMPLETE="${ROUND4_REQUIRE_ROUND3_COMPLETE:-1}"
CLAIM_ROOT="${OUTPUT_ROOT}/.state_claims"

for flag_name in EXECUTE PLAN_ONLY PRINT_HASHES REQUIRE_ROUND3_COMPLETE; do
    flag_value="${!flag_name}"
    [[ "${flag_value}" == 0 || "${flag_value}" == 1 ]] || {
        echo "${flag_name} must be 0 or 1" >&2
        exit 2
    }
done

declare -a OWNED_STATES=()
for block in "${STATE_BLOCKS[@]}"; do
    first_state="${block%%:*}"
    last_state="${block##*:}"
    for ((state=first_state; state<=last_state; state++)); do
        OWNED_STATES+=("${state}")
    done
done

NORMALIZED_STATE_SPEC="$({
    separator=""
    for state in "${OWNED_STATES[@]}"; do
        printf '%s%s' "${separator}" "${state}"
        separator=,
    done
})"

timestamp() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

print_plan() {
    echo "campaign=bottom_layer_round4_known_sif759_nosif_v1"
    echo "worker=${WORKER}"
    echo "expected_host=${EXPECTED_HOST}"
    echo "physical_gpu=${PHYSICAL_GPU}"
    echo "scene_class=${SCENE_CLASS}"
    echo "blocks=$(IFS=,; echo "${STATE_BLOCKS[*]}")"
    echo "states=${NORMALIZED_STATE_SPEC}"
    echo "measurement_classes=corrected,uncorrected"
    echo "perturbation_order=11,1:10"
    echo "sif_case_filter=off"
    echo "round4_sif_mode=off"
    echo "known_sif_wavelength_nm=759"
    echo "prior=$(realpath -m "${PRIOR_PATH}")"
    echo "source_prior=$(realpath -m "${SOURCE_PRIOR_PATH}")"
    echo "output=$(realpath -m "${OUTPUT_ROOT}")"
    echo "execute=${EXECUTE}"
}

if [[ "${PLAN_ONLY}" == 1 ]]; then
    print_plan
    exit 0
fi

for required in "${TRUTH_ROOT}/true_states.dat" "${PRIOR_PATH}" \
        "${SOURCE_PRIOR_PATH}" "${RUNNER}" "${PREFLIGHT}"; do
    [[ -f "${required}" ]] || {
        echo "missing required round-4 input or program: ${required}" >&2
        exit 1
    }
done

host="$(hostname -s)"
[[ "${host}" == "${EXPECTED_HOST}"* ]] || {
    echo "${WORKER} is reserved for ${EXPECTED_HOST}; current host is ${host}" >&2
    exit 1
}

preflight_env() {
    env \
        ROUND4_EXPECTED_STATES="${NORMALIZED_STATE_SPEC}" \
        ROUND4_EXPECTED_SCENE_CLASS="${SCENE_CLASS}" \
        BOTTOM_RETRIEVAL_CAMPAIGN_ROOT="${CAMPAIGN_ROOT}" \
        RETRIEVAL_PRIOR_PATH="${PRIOR_PATH}" \
        ROUND4_SOURCE_PRIOR_PATH="${SOURCE_PRIOR_PATH}" \
        RETRIEVAL_OUTPUT_ROOT="${OUTPUT_ROOT}" \
        ROUND4_REPO_ROOT="${REPO_ROOT}" \
        ROUND4_REQUIRE_ROUND3_COMPLETE="${REQUIRE_ROUND3_COMPLETE}" \
        ROUND4_PRIOR_SHA256="${ROUND4_PRIOR_SHA256:-}" \
        ROUND4_SOURCE_PRIOR_SHA256="${ROUND4_SOURCE_PRIOR_SHA256:-}" \
        ROUND4_INPUT_SET_SHA256="${ROUND4_INPUT_SET_SHA256:-}" \
        ROUND4_CODE_CHECKPOINT_SHA="${ROUND4_CODE_CHECKPOINT_SHA:-}" \
        ROUND4_CODESET_SHA256="${ROUND4_CODESET_SHA256:-}" \
        ROUND4_NOSIF_INITIALIZE_IDENTITY="${ROUND4_NOSIF_INITIALIZE_IDENTITY:-0}" \
        ROUND4_NOSIF_PRINT_CANDIDATE_HASHES="${ROUND4_NOSIF_PRINT_CANDIDATE_HASHES:-0}" \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        "${JULIA_BIN}" --project="${REPO_ROOT}" --startup-file=no "${PREFLIGHT}"
}

if [[ "${PRINT_HASHES}" == 1 ]]; then
    print_plan
    ROUND4_NOSIF_PRINT_CANDIDATE_HASHES=1 \
    ROUND4_NOSIF_INITIALIZE_IDENTITY=0 \
    preflight_env
    echo "Review and export the five values above before enabling execution."
    exit 0
fi

: "${ROUND4_PRIOR_SHA256:?export the caller-approved round-4 prior SHA-256}"
: "${ROUND4_SOURCE_PRIOR_SHA256:?export the caller-approved source-prior SHA-256}"
: "${ROUND4_INPUT_SET_SHA256:?export the caller-approved round-4 input-set SHA-256}"
: "${ROUND4_CODE_CHECKPOINT_SHA:?export the caller-approved 40-character Git checkpoint}"
: "${ROUND4_CODESET_SHA256:?export the caller-approved round-4 code-set SHA-256}"
for hash_name in ROUND4_PRIOR_SHA256 ROUND4_SOURCE_PRIOR_SHA256 \
        ROUND4_INPUT_SET_SHA256 ROUND4_CODESET_SHA256; do
    [[ "${!hash_name}" =~ ^[0-9a-f]{64}$ ]] || {
        echo "${hash_name} must be a lowercase 64-character SHA-256" >&2
        exit 2
    }
done
[[ "${ROUND4_CODE_CHECKPOINT_SHA}" =~ ^[0-9a-f]{40}$ ]] || {
    echo "ROUND4_CODE_CHECKPOINT_SHA must be a lowercase 40-character Git SHA" >&2
    exit 2
}

echo "$(timestamp) validating round-4 inputs and immutable identity"
ROUND4_NOSIF_INITIALIZE_IDENTITY="${EXECUTE}" preflight_env

if [[ "${EXECUTE}" != 1 ]]; then
    print_plan
    echo "$(timestamp) PRECHECK ONLY: no retrieval was started."
    echo "Set ROUND4_NOSIF_EXECUTE=1 only after approving this exact plan."
    exit 0
fi

IDENTITY="${OUTPUT_ROOT}/.control/campaign_identity.dat"
[[ -f "${IDENTITY}" ]] || {
    echo "preflight did not initialize the round-4 campaign identity" >&2
    exit 1
}
IDENTITY_SHA256="$(sha256sum "${IDENTITY}" | awk '{print $1}')"

assigned_gpu_uuid="$(nvidia-smi -i "${PHYSICAL_GPU}" --query-gpu=uuid \
    --format=csv,noheader 2>/dev/null)" || {
        echo "cannot resolve ${host} physical GPU ${PHYSICAL_GPU}" >&2
        exit 1
    }
assigned_gpu_uuid="${assigned_gpu_uuid//[[:space:]]/}"

require_assigned_gpu_idle() {
    local current_uuid applications row app_uuid
    current_uuid="$(nvidia-smi -i "${PHYSICAL_GPU}" --query-gpu=uuid \
        --format=csv,noheader 2>/dev/null)" || {
            echo "cannot re-resolve physical GPU ${PHYSICAL_GPU}" >&2
            return 1
        }
    current_uuid="${current_uuid//[[:space:]]/}"
    [[ "${current_uuid}" == "${assigned_gpu_uuid}" ]] || {
        echo "physical GPU ${PHYSICAL_GPU} identity changed" >&2
        return 1
    }
    applications="$(nvidia-smi \
        --query-compute-apps=gpu_uuid,pid,process_name \
        --format=csv,noheader,nounits 2>/dev/null)" || {
            echo "cannot inspect GPU compute processes" >&2
            return 1
        }
    while IFS= read -r row; do
        [[ -n "${row//[[:space:]]/}" ]] || continue
        app_uuid="${row%%,*}"
        app_uuid="${app_uuid//[[:space:]]/}"
        if [[ "${app_uuid}" == "${assigned_gpu_uuid}" ]]; then
            echo "physical GPU ${PHYSICAL_GPU} is busy: ${row}" >&2
            return 1
        fi
    done <<< "${applications}"
    return 0
}

require_assigned_gpu_idle

mkdir -p "${CLAIM_ROOT}"
declare -a ACQUIRED_CLAIMS=()
cleanup_claims() {
    local claim
    for claim in "${ACQUIRED_CLAIMS[@]}"; do
        rm -f "${claim}/owner.dat"
        rmdir "${claim}" 2>/dev/null || true
    done
}
trap cleanup_claims EXIT INT TERM

for state in "${OWNED_STATES[@]}"; do
    claim="${CLAIM_ROOT}/state$(printf '%03d' "${state}").claim"
    if ! mkdir "${claim}" 2>/dev/null; then
        echo "state ${state} is already claimed at ${claim}" >&2
        [[ -f "${claim}/owner.dat" ]] && sed -n '1,40p' \
            "${claim}/owner.dat" >&2
        exit 1
    fi
    ACQUIRED_CLAIMS+=("${claim}")
    {
        echo "campaign=${ROUND4_CAMPAIGN_ID:-bottom_layer_round4_known_sif759_nosif_v1}"
        echo "host=${host}"
        echo "pid=$$"
        echo "started_utc=$(timestamp)"
        echo "worker=${WORKER}"
        echo "physical_gpu=${PHYSICAL_GPU}"
        echo "state=${state}"
        echo "prior=$(realpath -m "${PRIOR_PATH}")"
        echo "prior_sha256=${ROUND4_PRIOR_SHA256}"
        echo "source_prior_sha256=${ROUND4_SOURCE_PRIOR_SHA256}"
        echo "input_set_sha256=${ROUND4_INPUT_SET_SHA256}"
        echo "code_checkpoint=${ROUND4_CODE_CHECKPOINT_SHA}"
        echo "codeset_sha256=${ROUND4_CODESET_SHA256}"
        echo "campaign_identity_sha256=${IDENTITY_SHA256}"
        echo "output=$(realpath -m "${OUTPUT_ROOT}")"
        echo "sif_case=off"
    } > "${claim}/owner.dat"
done

verify_immutable_inputs() {
    local current_identity
    printf '%s  %s\n' "${ROUND4_PRIOR_SHA256}" "${PRIOR_PATH}" |
        sha256sum --check --strict --quiet
    printf '%s  %s\n' "${ROUND4_SOURCE_PRIOR_SHA256}" "${SOURCE_PRIOR_PATH}" |
        sha256sum --check --strict --quiet
    current_identity="$(sha256sum "${IDENTITY}" | awk '{print $1}')"
    [[ "${current_identity}" == "${IDENTITY_SHA256}" ]] || {
        echo "round-4 campaign identity changed after preflight" >&2
        return 1
    }
    ROUND4_REQUIRE_ROUND3_COMPLETE=0 \
    ROUND4_NOSIF_INITIALIZE_IDENTITY=0 \
    preflight_env
}

run_block() {
    local first_state="$1"
    local last_state="$2"
    env \
        CUDA_VISIBLE_DEVICES="${PHYSICAL_GPU}" \
        CUDA_DEVICE=0 \
        RETRIEVAL_CLASS=paired \
        RETRIEVAL_ARCH=GPU \
        RETRIEVAL_FLOAT_TYPE=Float32 \
        RETRIEVAL_TRUTH_TABLE="${TRUTH_ROOT}/true_states.dat" \
        RETRIEVAL_MEASUREMENT_DIR="${OCO_ROOT}" \
        RETRIEVAL_NOISE_DIR="${NOISE_ROOT}" \
        RETRIEVAL_PRIOR_PATH="${PRIOR_PATH}" \
        RETRIEVAL_OUTPUT_ROOT="${OUTPUT_ROOT}" \
        RETRIEVAL_WRITE_MANIFEST=0 \
        ROUND4_SIF_MODE=off \
        ROUND4_KNOWN_WAVELENGTH_NM=759 \
        ROUND4_CODE_CHECKPOINT_SHA="${ROUND4_CODE_CHECKPOINT_SHA}" \
        ROUND4_CODESET_SHA256="${ROUND4_CODESET_SHA256}" \
        ROUND4_INPUT_SET_SHA256="${ROUND4_INPUT_SET_SHA256}" \
        ROUND4_CAMPAIGN_IDENTITY_SHA256="${IDENTITY_SHA256}" \
        SIF_CASE_FILTER=off \
        AEROSOL_CASE_FILTER="${SCENE_CLASS}" \
        FIRST_STATE="${first_state}" \
        LAST_STATE="${last_state}" \
        FIRST_PERTURBATION=1 \
        LAST_PERTURBATION=11 \
        FORCE=0 \
        FAIL_FAST=1 \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        JULIA_CPU_TARGET=native \
        "${JULIA_BIN}" --project="${REPO_ROOT}" --startup-file=no "${RUNNER}"
}

print_plan
for block in "${STATE_BLOCKS[@]}"; do
    first_state="${block%%:*}"
    last_state="${block##*:}"
    echo "$(timestamp) block_preflight worker=${WORKER} states=${first_state}-${last_state}"
    require_assigned_gpu_idle
    verify_immutable_inputs
    echo "$(timestamp) block_start worker=${WORKER} states=${first_state}-${last_state}"
    run_block "${first_state}" "${last_state}"
    verify_immutable_inputs
    echo "$(timestamp) block_pass worker=${WORKER} states=${first_state}-${last_state}"
done

echo "$(timestamp) round-4 no-SIF partition complete worker=${WORKER} states=${NORMALIZED_STATE_SPEC}"
