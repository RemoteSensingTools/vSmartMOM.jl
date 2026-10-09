#!/usr/bin/env bash

# Run one disjoint Curry-GPU partition of the bottom-layer CO2 clear-sky
# retrieval campaign. Both measurement classes and all eleven perturbation
# indices remain together on the same GPU for each truth state. Clear SIF-on
# work is released only after all clear and aerosol SIF-off states validate.

set -euo pipefail

PARTITION="${1:-}"
case "${PARTITION}" in
    curry0)
        PHYSICAL_GPU=0
        WRITE_MANIFEST=1
        NOSIF_BLOCKS=("1:5" "21:25")
        SIF_BLOCKS=("6:10" "26:30")
        OWNED_NOSIF_SPEC="1-5,21-25"
        OWNED_SIF_SPEC="6-10,26-30"
        EXPECTED_STATES=(
            001 002 003 004 005 006 007 008 009 010
            021 022 023 024 025 026 027 028 029 030)
        ;;
    curry1)
        PHYSICAL_GPU=1
        WRITE_MANIFEST=0
        NOSIF_BLOCKS=("41:45" "61:65")
        SIF_BLOCKS=("46:50" "66:70")
        OWNED_NOSIF_SPEC="41-45,61-65"
        OWNED_SIF_SPEC="46-50,66-70"
        EXPECTED_STATES=(
            041 042 043 044 045 046 047 048 049 050
            061 062 063 064 065 066 067 068 069 070)
        ;;
    *)
        echo "usage: $0 curry0|curry1" >&2
        exit 2
        ;;
esac

REPO_ROOT="${REPO_ROOT:-/home/sanghavi/code/github/uni_vSmartMOM}"
CAMPAIGN_ROOT="${REPO_ROOT}/RRS_XCO2/bottom_layer_XCO2_retrievals"
TRUTH_ROOT="${CAMPAIGN_ROOT}/truth"
OCO_ROOT="${TRUTH_ROOT}/OCO_radiances"
NOISE_ROOT="${OCO_ROOT}/noise_covariances"
PRIOR_PATH="${CAMPAIGN_ROOT}/retrieval_setup/apriori_states.nc"
OUTPUT_ROOT="${CAMPAIGN_ROOT}/retrievals"
RUNNER="${REPO_ROOT}/RRS_XCO2/inversion/run_retrievals.jl"
VALIDATOR="${REPO_ROOT}/RRS_XCO2/inversion/validate_bottom_layer_retrievals.jl"
POLL_SECONDS="${BOTTOM_RETRIEVAL_POLL_SECONDS:-600}"
WAIT_FOR_GPU="${BOTTOM_RETRIEVAL_WAIT_FOR_GPU:-1}"
JULIA_DEPOT_PATH_VALUE="${JULIA_DEPOT_PATH_VALUE:-${REPO_ROOT}/RRS_XCO2/.julia_depot_curry:/home/sanghavi/.julia}"
LOG_ROOT="${CAMPAIGN_ROOT}/logs"
CLAIM_ROOT="${OUTPUT_ROOT}/.partition_claims"
CLAIM_DIR="${CLAIM_ROOT}/${PARTITION}.claim"
ALL_NOSIF_SPEC="1-5,11-15,21-25,31-35,41-45,51-55,61-65,71-75"
ALL_SIF_SPEC="6-10,16-20,26-30,36-40,46-50,56-60,66-70,76-80"

timestamp() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

if [[ "$(hostname -s)" != curry* && "${BOTTOM_RETRIEVAL_ALLOW_OTHER_HOST:-0}" != 1 ]]; then
    echo "${PARTITION} is reserved for Curry; current host is $(hostname -s)" >&2
    exit 1
fi

for required in "${TRUTH_ROOT}/true_states.dat" "${PRIOR_PATH}" "${RUNNER}"; do
    [[ -f "${required}" ]] || {
        echo "missing required input: ${required}" >&2
        exit 1
    }
done

for state in "${EXPECTED_STATES[@]}"; do
    [[ -f "${TRUTH_ROOT}/hiressim_${state}.nc" ]] || {
        echo "missing truth scene ${state}" >&2
        exit 1
    }
    [[ -f "${OCO_ROOT}/OCO2sims_${state}.nc" ]] || {
        echo "missing OCO measurement ${state}" >&2
        exit 1
    }
    [[ -f "${NOISE_ROOT}/OCO2noise_${state}.nc" ]] || {
        echo "missing noise covariance ${state}" >&2
        exit 1
    }
done

gpu_is_busy() {
    local rows
    rows="$(nvidia-smi -i "${PHYSICAL_GPU}" \
        --query-compute-apps=pid,process_name \
        --format=csv,noheader,nounits 2>/dev/null || true)"
    [[ -z "${rows//[[:space:]]/}" ]] && return 1
    while IFS= read -r row; do
        [[ -z "${row//[[:space:]]/}" ]] && continue
        if [[ "${PHYSICAL_GPU}" == 0 && "${row}" == *"postgres: GPU0 memory keeper"* ]]; then
            continue
        fi
        return 0
    done <<< "${rows}"
    return 1
}

while gpu_is_busy; do
    [[ "${WAIT_FOR_GPU}" == 1 ]] || {
        echo "$(timestamp) physical GPU ${PHYSICAL_GPU} is busy; refusing launch" >&2
        exit 1
    }
    echo "$(timestamp) waiting: Curry physical GPU ${PHYSICAL_GPU} is busy"
    sleep "${POLL_SECONDS}"
done

mkdir -p "${LOG_ROOT}" "${CLAIM_ROOT}"
if ! mkdir "${CLAIM_DIR}" 2>/dev/null; then
    echo "partition claim already exists: ${CLAIM_DIR}" >&2
    echo "Inspect its owner record before deciding whether it is stale." >&2
    exit 1
fi
{
    echo "host=$(hostname -s)"
    echo "pid=$$"
    echo "started_utc=$(timestamp)"
    echo "partition=${PARTITION}"
    echo "physical_gpu=${PHYSICAL_GPU}"
    echo "owned_nosif_states=${OWNED_NOSIF_SPEC}"
    echo "owned_sif_states=${OWNED_SIF_SPEC}"
    echo "global_nosif_barrier=${ALL_NOSIF_SPEC}"
    echo "global_sif_barrier=${ALL_SIF_SPEC}"
    echo "perturbation_order=11,1:10"
    echo "measurement_classes=corrected,uncorrected"
} > "${CLAIM_DIR}/owner.txt"

cleanup_claim() {
    rm -f "${CLAIM_DIR}/owner.txt"
    rmdir "${CLAIM_DIR}" 2>/dev/null || true
}
trap cleanup_claim EXIT

run_subset() {
    local first_state="$1"
    local last_state="$2"
    local sif_filter="$3"
    local first_perturbation="$4"
    local last_perturbation="$5"
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
        RETRIEVAL_WRITE_MANIFEST="${WRITE_MANIFEST}" \
        SIF_CASE_FILTER="${sif_filter}" \
        AEROSOL_CASE_FILTER=none \
        FIRST_STATE="${first_state}" \
        LAST_STATE="${last_state}" \
        FIRST_PERTURBATION="${first_perturbation}" \
        LAST_PERTURBATION="${last_perturbation}" \
        FORCE=0 \
        FAIL_FAST=1 \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        JULIA_CPU_TARGET=native \
        julia --project="${REPO_ROOT}" "${RUNNER}"
}

validate_retrieval_products() {
    local state_spec="$1"
    local perturbation_spec="$2"
    local mode="$3"
    local scene_class="${4:-all}"
    env \
        EXPECTED_STATES="${state_spec}" \
        EXPECTED_PERTURBATIONS="${perturbation_spec}" \
        VALIDATION_MODE="${mode}" \
        BOTTOM_RETRIEVAL_SCENE_CLASS="${scene_class}" \
        BOTTOM_RETRIEVAL_CAMPAIGN_ROOT="${CAMPAIGN_ROOT}" \
        RETRIEVAL_OUTPUT_ROOT="${OUTPUT_ROOT}" \
        RETRIEVAL_PRIOR_PATH="${PRIOR_PATH}" \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        julia --project="${REPO_ROOT}" "${VALIDATOR}"
}

run_owned_blocks() {
    local phase="$1"
    local sif_filter="$2"
    shift 2
    local state_range first_state last_state
    for state_range in "$@"; do
        first_state="${state_range%%:*}"
        last_state="${state_range##*:}"
        echo "$(timestamp) block_start partition=${PARTITION} phase=${phase} states=${state_range}"
        run_subset "${first_state}" "${last_state}" "${sif_filter}" 1 11
        validate_retrieval_products "${first_state}-${last_state}" 1-11 products all
        echo "$(timestamp) block_pass partition=${PARTITION} phase=${phase} states=${state_range}"
    done
}

wait_for_global_phase() {
    local phase="$1"
    local state_spec="$2"
    echo "$(timestamp) global_phase_barrier_wait phase=${phase} states=${state_spec}"
    until validate_retrieval_products "${state_spec}" 1-11 products all \
            >/dev/null 2>&1; do
        sleep "${POLL_SECONDS}"
    done
    echo "$(timestamp) global_phase_barrier_pass phase=${phase}"
}

cd "${REPO_ROOT}"
echo "$(timestamp) launch partition=${PARTITION} physical_gpu=${PHYSICAL_GPU} " \
     "owned_nosif=${OWNED_NOSIF_SPEC} owned_sif=${OWNED_SIF_SPEC}"

# The first device performs a paired noiseless desert 400-ppm control smoke
# test. State 043 lies numerically in Curry1's range, so Curry1 waits for this
# validated result and then skips the two already-complete products. This
# prevents the cross-partition smoke solve from being computed twice.
if [[ "${PARTITION}" == curry0 && "${BOTTOM_RETRIEVAL_SKIP_SMOKE:-0}" != 1 ]]; then
    echo "$(timestamp) smoke_test state=043 perturbation=11 classes=paired"
    run_subset 43 43 off 11 11
    env \
        EXPECTED_STATES=43 \
        EXPECTED_PERTURBATIONS=11 \
        VALIDATION_MODE=smoke \
        RETRIEVAL_OUTPUT_ROOT="${OUTPUT_ROOT}" \
        RETRIEVAL_PRIOR_PATH="${PRIOR_PATH}" \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        julia --project="${REPO_ROOT}" "${VALIDATOR}"
fi

if [[ "${PARTITION}" == curry1 ]]; then
    echo "$(timestamp) waiting_for_validated_desert_smoke state=043 perturbation=11"
    until env \
        EXPECTED_STATES=43 \
        EXPECTED_PERTURBATIONS=11 \
        VALIDATION_MODE=smoke \
        RETRIEVAL_OUTPUT_ROOT="${OUTPUT_ROOT}" \
        RETRIEVAL_PRIOR_PATH="${PRIOR_PATH}" \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        julia --project="${REPO_ROOT}" "${VALIDATOR}" >/dev/null 2>&1; do
        sleep "${POLL_SECONDS}"
    done
    echo "$(timestamp) validated_desert_smoke_ready"
fi

echo "$(timestamp) production partition=${PARTITION} phase=nosif"
run_owned_blocks nosif off "${NOSIF_BLOCKS[@]}"
wait_for_global_phase nosif "${ALL_NOSIF_SPEC}"

echo "$(timestamp) production partition=${PARTITION} phase=sif"
run_owned_blocks sif on "${SIF_BLOCKS[@]}"
wait_for_global_phase sif "${ALL_SIF_SPEC}"

echo "$(timestamp) complete partition=${PARTITION} global_bottom_layer_retrievals=1760"
