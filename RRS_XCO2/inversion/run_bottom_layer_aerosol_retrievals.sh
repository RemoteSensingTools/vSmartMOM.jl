#!/usr/bin/env bash

# Guarded, statically partitioned Wurst queue for every bottom-layer CO2
# aerosol retrieval. All 40 clear and aerosol SIF-off states must validate
# across both hosts before either host may start a SIF-on state. Corrected and
# uncorrected retrievals and their shared perturbation draw stay adjacent on
# one worker.

set -euo pipefail

PARTITION="${1:-}"
case "${PARTITION}" in
    wurst0)
        PHYSICAL_GPU=0
        WRITE_MANIFEST=1
        NOSIF_BLOCKS=("11:15" "31:35" "51:51")
        SIF_BLOCKS=("16:20" "36:40")
        OWNED_NOSIF_SPEC="11-15,31-35,51"
        OWNED_SIF_SPEC="16-20,36-40"
        ;;
    wurst1)
        PHYSICAL_GPU=1
        WRITE_MANIFEST=0
        NOSIF_BLOCKS=("52:55" "71:75")
        SIF_BLOCKS=("56:60" "76:80")
        OWNED_NOSIF_SPEC="52-55,71-75"
        OWNED_SIF_SPEC="56-60,76-80"
        ;;
    *)
        echo "usage: $0 wurst0|wurst1" >&2
        exit 2
        ;;
esac

REPO_ROOT="${REPO_ROOT:-/home/sanghavi/code/github/uni_vSmartMOM}"
CAMPAIGN_ROOT="${BOTTOM_RETRIEVAL_CAMPAIGN_ROOT:-${REPO_ROOT}/RRS_XCO2/bottom_layer_XCO2_retrievals}"
TRUTH_ROOT="${CAMPAIGN_ROOT}/truth"
OCO_ROOT="${TRUTH_ROOT}/OCO_radiances"
NOISE_ROOT="${OCO_ROOT}/noise_covariances"
PRIOR_PATH="${CAMPAIGN_ROOT}/retrieval_setup/apriori_states.nc"
OUTPUT_ROOT="${CAMPAIGN_ROOT}/retrievals"
RUNNER="${REPO_ROOT}/RRS_XCO2/inversion/run_retrievals.jl"
PREFLIGHT="${REPO_ROOT}/RRS_XCO2/inversion/preflight_bottom_layer_aerosol_retrievals.jl"
RETRIEVAL_VALIDATOR="${REPO_ROOT}/RRS_XCO2/inversion/validate_bottom_layer_retrievals.jl"
RADIANCE_VALIDATOR="${REPO_ROOT}/RRS_XCO2/inversion/instrument/validate_oco_radiances.jl"
NOISE_VALIDATOR="${REPO_ROOT}/RRS_XCO2/inversion/instrument/validate_noise_covariances.jl"
POLL_SECONDS="${BOTTOM_RETRIEVAL_POLL_SECONDS:-600}"
WAIT_FOR_GPU="${BOTTOM_RETRIEVAL_WAIT_FOR_GPU:-1}"
JULIA_DEPOT_PATH_VALUE="${JULIA_DEPOT_PATH_VALUE:-${REPO_ROOT}/RRS_XCO2/.julia_depot_wurst:/home/sanghavi/.julia}"
CLAIM_ROOT="${OUTPUT_ROOT}/.partition_claims"
CLAIM_DIR="${CLAIM_ROOT}/${PARTITION}_aerosol_all_sif.claim"
GATE_ROOT="${OUTPUT_ROOT}/.aerosol_release_gates"
MANIFEST_PATH="${OUTPUT_ROOT}/retrieval_manifest_wurst0_aerosol_all.dat"

ALL_NOSIF_SPEC="1-5,11-15,21-25,31-35,41-45,51-55,61-65,71-75"
ALL_SIF_SPEC="6-10,16-20,26-30,36-40,46-50,56-60,66-70,76-80"

timestamp() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

if [[ "$(hostname -s)" != wurst* && "${BOTTOM_RETRIEVAL_ALLOW_OTHER_HOST:-0}" != 1 ]]; then
    echo "${PARTITION} is reserved for Wurst; current host is $(hostname -s)." >&2
    exit 1
fi

for required in "${TRUTH_ROOT}/true_states.dat" "${PRIOR_PATH}" \
        "${RUNNER}" "${PREFLIGHT}" "${RETRIEVAL_VALIDATOR}" \
        "${RADIANCE_VALIDATOR}" "${NOISE_VALIDATOR}"; do
    [[ -f "${required}" ]] || {
        echo "missing required input or program: ${required}" >&2
        exit 1
    }
done

validate_instrument_products() {
    echo "$(timestamp) validate_complete_80_state_oco_dataset"
    env \
        EXPECTED_STATES=1-80 \
        SYNTHETIC_OCO_OUT="${OCO_ROOT}" \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        julia --project="${REPO_ROOT}" "${RADIANCE_VALIDATOR}"
    echo "$(timestamp) validate_complete_80_state_noise_dataset"
    env \
        EXPECTED_STATES=1-80 \
        SYNTHETIC_OCO_DIR="${OCO_ROOT}" \
        OCO_NOISE_OUT="${NOISE_ROOT}" \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        julia --project="${REPO_ROOT}" "${NOISE_VALIDATOR}"
    echo "$(timestamp) validate_aerosol_retrieval_selection"
    env \
        BOTTOM_RETRIEVAL_CAMPAIGN_ROOT="${CAMPAIGN_ROOT}" \
        BOTTOM_AEROSOL_MANIFEST_PATH="${MANIFEST_PATH}" \
        BOTTOM_AEROSOL_WRITE_MANIFEST="${WRITE_MANIFEST}" \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        julia --project="${REPO_ROOT}" "${PREFLIGHT}"
}

cd "${REPO_ROOT}"
validate_instrument_products
if [[ "${BOTTOM_RETRIEVAL_PREFLIGHT_ONLY:-0}" == 1 ]]; then
    echo "$(timestamp) preflight_only_complete partition=${PARTITION} " \
         "owned_nosif=${OWNED_NOSIF_SPEC} owned_sif=${OWNED_SIF_SPEC}"
    exit 0
fi

gpu_uuid="$(nvidia-smi -i "${PHYSICAL_GPU}" --query-gpu=uuid \
    --format=csv,noheader 2>/dev/null)" || {
        echo "cannot resolve Wurst physical GPU ${PHYSICAL_GPU}" >&2
        exit 1
    }

gpu_is_busy() {
    local apps
    apps="$(nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name \
        --format=csv,noheader,nounits 2>/dev/null)" || return 0
    [[ -n "$(awk -F, -v uuid="${gpu_uuid}" '$1 == uuid { print }' <<< "${apps}")" ]]
}

mkdir -p "${CLAIM_ROOT}" "${GATE_ROOT}"
if ! mkdir "${CLAIM_DIR}" 2>/dev/null; then
    echo "partition claim already exists: ${CLAIM_DIR}" >&2
    [[ -f "${CLAIM_DIR}/owner.txt" ]] && sed -n '1,30p' "${CLAIM_DIR}/owner.txt" >&2
    exit 1
fi
{
    echo "host=$(hostname -s)"
    echo "pid=$$"
    echo "started_utc=$(timestamp)"
    echo "partition=${PARTITION}"
    echo "physical_gpu=${PHYSICAL_GPU}"
    echo "physical_gpu_uuid=${gpu_uuid}"
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
trap cleanup_claim EXIT INT TERM

while gpu_is_busy; do
    [[ "${WAIT_FOR_GPU}" == 1 ]] || {
        echo "$(timestamp) Wurst physical GPU ${PHYSICAL_GPU} is busy; refusing launch" >&2
        exit 1
    }
    echo "$(timestamp) waiting: Wurst physical GPU ${PHYSICAL_GPU} is busy"
    sleep "${POLL_SECONDS}"
done

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
        RETRIEVAL_WRITE_MANIFEST=0 \
        SIF_CASE_FILTER="${sif_filter}" \
        AEROSOL_CASE_FILTER=aerosol \
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
    local scene_class="${4:-aerosol_all_sif}"
    env \
        EXPECTED_STATES="${state_spec}" \
        EXPECTED_PERTURBATIONS="${perturbation_spec}" \
        VALIDATION_MODE="${mode}" \
        BOTTOM_RETRIEVAL_SCENE_CLASS="${scene_class}" \
        BOTTOM_RETRIEVAL_CAMPAIGN_ROOT="${CAMPAIGN_ROOT}" \
        RETRIEVAL_OUTPUT_ROOT="${OUTPUT_ROOT}" \
        RETRIEVAL_PRIOR_PATH="${PRIOR_PATH}" \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        julia --project="${REPO_ROOT}" "${RETRIEVAL_VALIDATOR}"
}

wait_for_validated_smoke() {
    local phase="$1"
    local state="$2"
    local gate="${GATE_ROOT}/${phase}_smoke.complete"
    echo "$(timestamp) waiting_for_smoke phase=${phase} state=${state}"
    until [[ -f "${gate}" ]] && \
            validate_retrieval_products "${state}" 11 smoke >/dev/null 2>&1; do
        sleep "${POLL_SECONDS}"
    done
    echo "$(timestamp) smoke_from_wurst0_validated phase=${phase} state=${state}"
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
        validate_retrieval_products "${first_state}-${last_state}" 1-11 products
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

publish_smoke_gate() {
    local phase="$1"
    local state="$2"
    local gate="${GATE_ROOT}/${phase}_smoke.complete"
    local temporary="${gate}.tmp.$$"
    {
        echo "validated_utc=$(timestamp)"
        echo "state=${state}"
        echo "prior_sha256=$(sha256sum "${PRIOR_PATH}" | awk '{print $1}')"
    } > "${temporary}"
    mv "${temporary}" "${gate}"
}

echo "$(timestamp) launch partition=${PARTITION} physical_gpu=${PHYSICAL_GPU} " \
     "owned_nosif=${OWNED_NOSIF_SPEC} owned_sif=${OWNED_SIF_SPEC}"

if [[ "${PARTITION}" == wurst0 ]]; then
    echo "$(timestamp) smoke_start phase=nosif state=013 perturbation=11"
    run_subset 13 13 off 11 11
    validate_retrieval_products 13 11 smoke
    publish_smoke_gate nosif 13
    echo "$(timestamp) smoke_pass phase=nosif state=013"
else
    wait_for_validated_smoke nosif 13
fi

run_owned_blocks nosif off "${NOSIF_BLOCKS[@]}"
wait_for_global_phase nosif "${ALL_NOSIF_SPEC}"

# The global no-SIF barrier above is deliberately before the SIF smoke gate.
# Consequently no worker can start state 018 (or any other SIF-on state) while
# even one of the 40 clear or aerosol no-SIF states remains incomplete or invalid.
if [[ "${PARTITION}" == wurst0 ]]; then
    echo "$(timestamp) smoke_start phase=sif state=018 perturbation=11"
    run_subset 18 18 on 11 11
    validate_retrieval_products 18 11 smoke
    publish_smoke_gate sif 18
    echo "$(timestamp) smoke_pass phase=sif state=018"
else
    wait_for_validated_smoke sif 18
fi

run_owned_blocks sif on "${SIF_BLOCKS[@]}"
wait_for_global_phase sif "${ALL_SIF_SPEC}"

echo "$(timestamp) complete partition=${PARTITION} global_aerosol_retrievals=880"
