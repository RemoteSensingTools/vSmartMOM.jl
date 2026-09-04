#!/usr/bin/env bash

# Run one explicitly assigned Curry/Wurst partition of the isolated
# bottom-layer retrieval round that uses the tapered mapped-ACOS CO2 prior.
# This launcher is permanently restricted to SIF-off truth states. It is
# fail-closed: without TAPERED_NOSIF_EXECUTE=1 it performs validation only.
# Every non-plan invocation also requires caller-approved
# TAPERED_PRIOR_SHA256 and CURRENT_PRIOR_SHA256 values; these become part of
# the immutable campaign identity shared by all Curry/Wurst workers.

set -euo pipefail

WORKER="${1:-}"
STATE_SPEC="${2:-${TAPERED_NOSIF_STATE_BLOCKS:-}}"

case "${WORKER}" in
    curry0)
        EXPECTED_HOST=curry
        PHYSICAL_GPU=0
        DEPOT_NAME=.julia_depot_curry
        ;;
    curry1)
        EXPECTED_HOST=curry
        PHYSICAL_GPU=1
        DEPOT_NAME=.julia_depot_curry
        ;;
    wurst0)
        EXPECTED_HOST=wurst
        PHYSICAL_GPU=0
        DEPOT_NAME=.julia_depot_wurst
        ;;
    wurst1)
        EXPECTED_HOST=wurst
        PHYSICAL_GPU=1
        DEPOT_NAME=.julia_depot_wurst
        ;;
    *)
        echo "usage: $0 curry0|curry1|wurst0|wurst1 STATE-RANGES" >&2
        echo "example state ranges: 1-5,21-25" >&2
        exit 2
        ;;
esac

[[ -n "${STATE_SPEC}" ]] || {
    echo "A comma-separated no-SIF state assignment is required." >&2
    echo "Example: $0 ${WORKER} 1-5,21-25" >&2
    exit 2
}

SCRIPT_DIRECTORY="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="${REPO_ROOT:-$(git -C "${SCRIPT_DIRECTORY}" rev-parse --show-toplevel)}"
DEFAULT_DEPOT="${REPO_ROOT}/RRS_XCO2/${DEPOT_NAME}:${HOME}/.julia"
CAMPAIGN_ROOT="${BOTTOM_RETRIEVAL_CAMPAIGN_ROOT:-${REPO_ROOT}/RRS_XCO2/bottom_layer_XCO2_retrievals}"
TRUTH_ROOT="${CAMPAIGN_ROOT}/truth"
OCO_ROOT="${TRUTH_ROOT}/OCO_radiances"
NOISE_ROOT="${OCO_ROOT}/noise_covariances"
MODEL_NAME="acos_mapped_tapered_vertical_correlation"
PRIOR_PATH="${TAPERED_NOSIF_PRIOR_PATH:-${CAMPAIGN_ROOT}/retrieval_setup/apriori_states_${MODEL_NAME}.nc}"
OUTPUT_ROOT="${TAPERED_NOSIF_OUTPUT_ROOT:-${CAMPAIGN_ROOT}/retrievals_${MODEL_NAME}_nosif}"
ACTIVE_PRIOR="${CAMPAIGN_ROOT}/retrieval_setup/apriori_states.nc"
ACTIVE_OUTPUT="${CAMPAIGN_ROOT}/retrievals"
RUNNER="${REPO_ROOT}/RRS_XCO2/inversion/run_retrievals.jl"
PREFLIGHT="${REPO_ROOT}/RRS_XCO2/inversion/preflight_tapered_no_sif_retrievals.jl"
VALIDATOR="${REPO_ROOT}/RRS_XCO2/inversion/validate_bottom_layer_retrievals.jl"
JULIA_BIN="${JULIA_BIN:-julia}"
JULIA_DEPOT_PATH_VALUE="${JULIA_DEPOT_PATH_VALUE:-${DEFAULT_DEPOT}}"
EXECUTE="${TAPERED_NOSIF_EXECUTE:-0}"
PLAN_ONLY="${TAPERED_NOSIF_PLAN_ONLY:-0}"
REQUIRE_CURRENT_COMPLETE="${TAPERED_NOSIF_REQUIRE_CURRENT_COMPLETE:-1}"
CLAIM_ROOT="${OUTPUT_ROOT}/.state_claims"
ALL_NOSIF_SPEC="1-5,11-15,21-25,31-35,41-45,51-55,61-65,71-75"

for flag_name in EXECUTE PLAN_ONLY REQUIRE_CURRENT_COMPLETE; do
    flag_value="${!flag_name}"
    [[ "${flag_value}" == 0 || "${flag_value}" == 1 ]] || {
        echo "${flag_name} must be 0 or 1" >&2
        exit 2
    }
done

declare -a STATE_BLOCKS=()
declare -a OWNED_STATES=()
declare -A SEEN_STATES=()
IFS=',' read -r -a RAW_BLOCKS <<< "${STATE_SPEC}"
for raw_block in "${RAW_BLOCKS[@]}"; do
    block="${raw_block//[[:space:]]/}"
    if [[ "${block}" =~ ^([0-9]+)(-([0-9]+))?$ ]]; then
        first_state="${BASH_REMATCH[1]}"
        last_state="${BASH_REMATCH[3]:-${first_state}}"
    else
        echo "invalid state range: ${raw_block}" >&2
        exit 2
    fi
    (( first_state <= last_state )) || {
        echo "descending state range is not allowed: ${block}" >&2
        exit 2
    }
    for ((state=first_state; state<=last_state; state++)); do
        (( state >= 1 && state <= 80 )) || {
            echo "state ${state} is outside the bottom-layer campaign" >&2
            exit 2
        }
        decade_offset=$(( (state - 1) % 10 + 1 ))
        (( decade_offset <= 5 )) || {
            echo "state ${state} is SIF-on; this launcher accepts only no-SIF states" >&2
            exit 2
        }
        [[ -z "${SEEN_STATES[${state}]:-}" ]] || {
            echo "state ${state} is assigned more than once" >&2
            exit 2
        }
        SEEN_STATES[${state}]=1
        OWNED_STATES+=("${state}")
    done
    STATE_BLOCKS+=("${first_state}:${last_state}")
done

NORMALIZED_STATE_SPEC="$({
    separator=""
    for state in "${OWNED_STATES[@]}"; do
        printf '%s%s' "${separator}" "${state}"
        separator=,
    done
})"

canonical_active_prior="$(realpath -m "${ACTIVE_PRIOR}")"
canonical_selected_prior="$(realpath -m "${PRIOR_PATH}")"
canonical_active_output="$(realpath -m "${ACTIVE_OUTPUT}")"
canonical_selected_output="$(realpath -m "${OUTPUT_ROOT}")"
[[ "${canonical_selected_prior}" != "${canonical_active_prior}" ]] || {
    echo "refusing to use the active campaign prior: ${ACTIVE_PRIOR}" >&2
    exit 1
}
[[ "${canonical_selected_output}" != "${canonical_active_output}" &&
   "${canonical_selected_output}" != "${canonical_active_output}/"* ]] || {
    echo "refusing to write in or beneath the active output: ${ACTIVE_OUTPUT}" >&2
    exit 1
}

timestamp() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

print_plan() {
    echo "worker=${WORKER}"
    echo "expected_host=${EXPECTED_HOST}"
    echo "physical_gpu=${PHYSICAL_GPU}"
    echo "states=${NORMALIZED_STATE_SPEC}"
    echo "sif_case_filter=off"
    echo "measurement_classes=corrected,uncorrected"
    echo "perturbation_order=11,1:10"
    echo "prior=${canonical_selected_prior}"
    echo "output=${canonical_selected_output}"
    echo "execute=${EXECUTE}"
}

if [[ "${PLAN_ONLY}" == 1 ]]; then
    print_plan
    exit 0
fi

: "${TAPERED_PRIOR_SHA256:?export the caller-approved tapered-prior SHA-256}"
: "${CURRENT_PRIOR_SHA256:?export the caller-approved current-prior SHA-256}"
[[ "${TAPERED_PRIOR_SHA256}" =~ ^[0-9a-f]{64}$ ]] || {
    echo "TAPERED_PRIOR_SHA256 must be a lowercase 64-character SHA-256" >&2
    exit 2
}
[[ "${CURRENT_PRIOR_SHA256}" =~ ^[0-9a-f]{64}$ ]] || {
    echo "CURRENT_PRIOR_SHA256 must be a lowercase 64-character SHA-256" >&2
    exit 2
}

host="$(hostname -s)"
[[ "${host}" == "${EXPECTED_HOST}"* ||
   "${TAPERED_NOSIF_ALLOW_OTHER_HOST:-0}" == 1 ]] || {
    echo "${WORKER} is reserved for ${EXPECTED_HOST}; current host is ${host}" >&2
    exit 1
}

for required in "${TRUTH_ROOT}/true_states.dat" "${PRIOR_PATH}" \
        "${RUNNER}" "${PREFLIGHT}" "${VALIDATOR}"; do
    [[ -f "${required}" ]] || {
        echo "missing required input or program: ${required}" >&2
        exit 1
    }
done

if [[ "${REQUIRE_CURRENT_COMPLETE}" == 1 ]]; then
    active_claim_root="${ACTIVE_OUTPUT}/.partition_claims"
    if [[ -d "${active_claim_root}" ]] &&
            find "${active_claim_root}" -mindepth 1 -maxdepth 1 -type d \
                -print -quit | grep -q .; then
        echo "active campaign still owns one or more partition claims:" >&2
        find "${active_claim_root}" -mindepth 1 -maxdepth 2 -type f \
            -name owner.txt -print >&2
        exit 1
    fi
    echo "$(timestamp) validating completion of the current no-SIF round"
    env \
        EXPECTED_STATES="${ALL_NOSIF_SPEC}" \
        EXPECTED_PERTURBATIONS=1-11 \
        VALIDATION_MODE=products \
        BOTTOM_RETRIEVAL_SCENE_CLASS=all \
        BOTTOM_RETRIEVAL_CAMPAIGN_ROOT="${CAMPAIGN_ROOT}" \
        RETRIEVAL_OUTPUT_ROOT="${ACTIVE_OUTPUT}" \
        RETRIEVAL_PRIOR_PATH="${ACTIVE_PRIOR}" \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        "${JULIA_BIN}" --project="${REPO_ROOT}" --startup-file=no "${VALIDATOR}"
fi

echo "$(timestamp) validating tapered no-SIF inputs"
env \
    EXPECTED_STATES="${NORMALIZED_STATE_SPEC}" \
    BOTTOM_RETRIEVAL_CAMPAIGN_ROOT="${CAMPAIGN_ROOT}" \
    RETRIEVAL_PRIOR_PATH="${PRIOR_PATH}" \
    CURRENT_PRIOR_PATH="${ACTIVE_PRIOR}" \
    TAPERED_PRIOR_SHA256="${TAPERED_PRIOR_SHA256}" \
    CURRENT_PRIOR_SHA256="${CURRENT_PRIOR_SHA256}" \
    RETRIEVAL_OUTPUT_ROOT="${OUTPUT_ROOT}" \
    TAPERED_NOSIF_INITIALIZE_IDENTITY="${EXECUTE}" \
    JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
    "${JULIA_BIN}" --project="${REPO_ROOT}" --startup-file=no "${PREFLIGHT}"

if [[ "${EXECUTE}" != 1 ]]; then
    print_plan
    echo "$(timestamp) PRECHECK ONLY: no retrieval was started."
    echo "Set TAPERED_NOSIF_EXECUTE=1 only after approving this assignment."
    exit 0
fi

IDENTITY="${OUTPUT_ROOT}/.control/campaign_identity.dat"
[[ -f "${IDENTITY}" ]] || {
    echo "preflight did not initialize the tapered no-SIF campaign identity" >&2
    exit 1
}
IDENTITY_SHA256="$(sha256sum "${IDENTITY}" | awk '{print $1}')"

gpu_uuid="$(nvidia-smi -i "${PHYSICAL_GPU}" --query-gpu=uuid \
    --format=csv,noheader 2>/dev/null)" || {
        echo "cannot resolve ${host} physical GPU ${PHYSICAL_GPU}" >&2
        exit 1
    }

gpu_is_busy() {
    local applications row app_uuid process_name
    applications="$(nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name \
        --format=csv,noheader,nounits 2>/dev/null)" || return 0
    while IFS= read -r row; do
        [[ -n "${row//[[:space:]]/}" ]] || continue
        app_uuid="${row%%,*}"
        [[ "${app_uuid//[[:space:]]/}" == "${gpu_uuid//[[:space:]]/}" ]] ||
            continue
        process_name="${row##*,}"
        if [[ "${WORKER}" == curry0 &&
              "${process_name}" == *"postgres: GPU0 memory keeper"* ]]; then
            # Curry GPU0 has a small permanent PostgreSQL memory keeper. It
            # is not a retrieval/RT worker and must not reserve the device.
            continue
        fi
        return 0
    done <<< "${applications}"
    return 1
}

if gpu_is_busy; then
    echo "physical GPU ${PHYSICAL_GPU} is busy; refusing the tapered launch" >&2
    exit 1
fi

mkdir -p "${CLAIM_ROOT}"
declare -a ACQUIRED_CLAIMS=()
cleanup_claims() {
    local claim
    for claim in "${ACQUIRED_CLAIMS[@]}"; do
        rm -f "${claim}/owner.txt"
        rmdir "${claim}" 2>/dev/null || true
    done
}
trap cleanup_claims EXIT INT TERM

for state in "${OWNED_STATES[@]}"; do
    claim="${CLAIM_ROOT}/state$(printf '%03d' "${state}").claim"
    if ! mkdir "${claim}" 2>/dev/null; then
        echo "state ${state} is already claimed at ${claim}" >&2
        [[ -f "${claim}/owner.txt" ]] && sed -n '1,30p' \
            "${claim}/owner.txt" >&2
        exit 1
    fi
    ACQUIRED_CLAIMS+=("${claim}")
    {
        echo "host=${host}"
        echo "pid=$$"
        echo "started_utc=$(timestamp)"
        echo "worker=${WORKER}"
        echo "physical_gpu=${PHYSICAL_GPU}"
        echo "state=${state}"
        echo "prior=${canonical_selected_prior}"
        echo "prior_sha256=${TAPERED_PRIOR_SHA256}"
        echo "reference_prior_sha256=${CURRENT_PRIOR_SHA256}"
        echo "campaign_identity_sha256=${IDENTITY_SHA256}"
        echo "output=${canonical_selected_output}"
        echo "sif_case_filter=off"
    } > "${claim}/owner.txt"
done

run_block() {
    local first_state="$1"
    local last_state="$2"
    local current_identity_sha256
    printf '%s  %s\n' "${TAPERED_PRIOR_SHA256}" "${PRIOR_PATH}" |
        sha256sum --check --strict --quiet
    printf '%s  %s\n' "${CURRENT_PRIOR_SHA256}" "${ACTIVE_PRIOR}" |
        sha256sum --check --strict --quiet
    current_identity_sha256="$(sha256sum "${IDENTITY}" | awk '{print $1}')"
    [[ "${current_identity_sha256}" == "${IDENTITY_SHA256}" ]] || {
        echo "tapered no-SIF campaign identity changed after preflight" >&2
        exit 1
    }
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
        SIF_CASE_FILTER=off \
        AEROSOL_CASE_FILTER=all \
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

validate_block() {
    local first_state="$1"
    local last_state="$2"
    env \
        EXPECTED_STATES="${first_state}-${last_state}" \
        EXPECTED_PERTURBATIONS=1-11 \
        VALIDATION_MODE=products \
        BOTTOM_RETRIEVAL_SCENE_CLASS=all \
        BOTTOM_RETRIEVAL_CAMPAIGN_ROOT="${CAMPAIGN_ROOT}" \
        RETRIEVAL_OUTPUT_ROOT="${OUTPUT_ROOT}" \
        RETRIEVAL_PRIOR_PATH="${PRIOR_PATH}" \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        "${JULIA_BIN}" --project="${REPO_ROOT}" --startup-file=no "${VALIDATOR}"
}

print_plan
for block in "${STATE_BLOCKS[@]}"; do
    first_state="${block%%:*}"
    last_state="${block##*:}"
    echo "$(timestamp) block_start worker=${WORKER} states=${first_state}-${last_state}"
    run_block "${first_state}" "${last_state}"
    validate_block "${first_state}" "${last_state}"
    echo "$(timestamp) block_pass worker=${WORKER} states=${first_state}-${last_state}"
done

echo "$(timestamp) tapered no-SIF partition complete worker=${WORKER} states=${NORMALIZED_STATE_SPEC}"
