#!/usr/bin/env bash

# Isolated noiseless retrieval test for sigma(ln(AOD760)) = 0.75.
#
# This script never reads from or writes to the production retrieval output
# directories.  It temporarily claims Wurst GPU0's production partition so a
# normal production launcher cannot race this sensitivity test, but it does
# not publish either production smoke-release gate.

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/sanghavi/code/github/uni_vSmartMOM}"
CAMPAIGN_ROOT="${REPO_ROOT}/RRS_XCO2/bottom_layer_XCO2_retrievals"
TRUTH_ROOT="${CAMPAIGN_ROOT}/truth"
OCO_ROOT="${TRUTH_ROOT}/OCO_radiances"
NOISE_ROOT="${OCO_ROOT}/noise_covariances"
TEST_ROOT="${CAMPAIGN_ROOT}/sensitivity_tests/ln_aod_sigma_0p75"
PRIOR_PATH="${TEST_ROOT}/retrieval_setup/apriori_states.nc"
OUTPUT_ROOT="${TEST_ROOT}/retrievals"
RUNNER="${REPO_ROOT}/RRS_XCO2/inversion/run_retrievals.jl"
VALIDATOR="${REPO_ROOT}/RRS_XCO2/inversion/validate_bottom_layer_retrievals.jl"
PHYSICAL_GPU=0
JULIA_DEPOT_PATH_VALUE="${JULIA_DEPOT_PATH_VALUE:-${REPO_ROOT}/RRS_XCO2/.julia_depot_wurst:/home/sanghavi/.julia}"
PRODUCTION_CLAIM="${CAMPAIGN_ROOT}/retrievals/.partition_claims/wurst0_aerosol_all_sif.claim"
TEST_CLAIM="${TEST_ROOT}/.gpu0_test.claim"

timestamp() {
    date -u +%Y-%m-%dT%H:%M:%SZ
}

if [[ "$(hostname -s)" != wurst* && "${BOTTOM_RETRIEVAL_ALLOW_OTHER_HOST:-0}" != 1 ]]; then
    echo "This sensitivity test is reserved for Wurst; current host is $(hostname -s)." >&2
    exit 1
fi

for required in "${TRUTH_ROOT}/true_states.dat" "${PRIOR_PATH}" \
        "${OCO_ROOT}/OCO2sims_001.nc" "${OCO_ROOT}/OCO2sims_013.nc" \
        "${NOISE_ROOT}/OCO2noise_001.nc" "${NOISE_ROOT}/OCO2noise_013.nc" \
        "${RUNNER}" "${VALIDATOR}"; do
    [[ -f "${required}" ]] || {
        echo "missing required test input: ${required}" >&2
        exit 1
    }
done

mkdir -p "$(dirname "${PRODUCTION_CLAIM}")" "${TEST_ROOT}"
if ! mkdir "${TEST_CLAIM}" 2>/dev/null; then
    echo "sensitivity-test claim already exists: ${TEST_CLAIM}" >&2
    [[ -f "${TEST_CLAIM}/owner.txt" ]] && sed -n '1,30p' "${TEST_CLAIM}/owner.txt" >&2
    exit 1
fi
OWN_TEST_CLAIM=1
OWN_PRODUCTION_CLAIM=0

cleanup_claims() {
    if [[ "${OWN_TEST_CLAIM}" == 1 ]]; then
        rm -f "${TEST_CLAIM}/owner.txt"
        rmdir "${TEST_CLAIM}" 2>/dev/null || true
    fi
    if [[ "${OWN_PRODUCTION_CLAIM}" == 1 ]]; then
        rm -f "${PRODUCTION_CLAIM}/owner.txt"
        rmdir "${PRODUCTION_CLAIM}" 2>/dev/null || true
    fi
}
trap cleanup_claims EXIT INT TERM

if ! mkdir "${PRODUCTION_CLAIM}" 2>/dev/null; then
    echo "Wurst GPU0 production partition is already claimed:" >&2
    [[ -f "${PRODUCTION_CLAIM}/owner.txt" ]] && \
        sed -n '1,30p' "${PRODUCTION_CLAIM}/owner.txt" >&2
    exit 1
fi
OWN_PRODUCTION_CLAIM=1

gpu_uuid="$(nvidia-smi -i "${PHYSICAL_GPU}" --query-gpu=uuid \
    --format=csv,noheader 2>/dev/null)" || {
        echo "cannot resolve Wurst physical GPU ${PHYSICAL_GPU}" >&2
        exit 1
    }
apps="$(nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name \
    --format=csv,noheader,nounits 2>/dev/null || true)"
if [[ -n "$(awk -F, -v uuid="${gpu_uuid}" '$1 == uuid { print }' <<< "${apps}")" ]]; then
    echo "Wurst physical GPU ${PHYSICAL_GPU} became busy; refusing test launch." >&2
    exit 1
fi

for claim in "${TEST_CLAIM}" "${PRODUCTION_CLAIM}"; do
    {
        echo "host=$(hostname -s)"
        echo "pid=$$"
        echo "started_utc=$(timestamp)"
        echo "purpose=ln_aod_sigma_0p75_sensitivity"
        echo "physical_gpu=${PHYSICAL_GPU}"
        echo "physical_gpu_uuid=${gpu_uuid}"
        echo "states=001,013"
        echo "perturbations=11"
        echo "measurement_classes=corrected,uncorrected"
        echo "production_outputs_modified=0"
        echo "production_release_gates_modified=0"
    } > "${claim}/owner.txt"
done

echo "$(timestamp) begin ln_aod_sigma_0p75 sensitivity test"
echo "prior=${PRIOR_PATH}"
echo "prior_sha256=$(sha256sum "${PRIOR_PATH}" | awk '{print $1}')"
echo "output_root=${OUTPUT_ROOT}"

run_state() {
    local state="$1"
    local aerosol_filter="$2"
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
        AEROSOL_CASE_FILTER="${aerosol_filter}" \
        FIRST_STATE="${state}" \
        LAST_STATE="${state}" \
        FIRST_PERTURBATION=11 \
        LAST_PERTURBATION=11 \
        FORCE=0 \
        FAIL_FAST=1 \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        JULIA_CPU_TARGET=native \
        julia --project="${REPO_ROOT}" "${RUNNER}"
}

validate_state() {
    local state="$1"
    local scene_class="$2"
    env \
        EXPECTED_STATES="${state}" \
        EXPECTED_PERTURBATIONS=11 \
        VALIDATION_MODE=products \
        BOTTOM_RETRIEVAL_SCENE_CLASS="${scene_class}" \
        BOTTOM_RETRIEVAL_CAMPAIGN_ROOT="${CAMPAIGN_ROOT}" \
        RETRIEVAL_OUTPUT_ROOT="${OUTPUT_ROOT}" \
        RETRIEVAL_PRIOR_PATH="${PRIOR_PATH}" \
        JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH_VALUE}" \
        julia --project="${REPO_ROOT}" "${VALIDATOR}"
}

echo "$(timestamp) run state=001 aerosol=none sif=off perturbation=11"
run_state 1 none
validate_state 1 clear_nosif

echo "$(timestamp) run state=013 aerosol=aod760_0p28 sif=off perturbation=11"
run_state 13 aerosol
validate_state 13 aerosol_all_sif

echo "$(timestamp) sensitivity test complete"
