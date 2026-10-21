#!/usr/bin/env bash

# Gattaca2-only bottom-layer SIF-on retrievals with the tapered ACOS CO2 prior.
#
# This launcher cannot select no-SIF truth and cannot write to the historical
# Curry/Wurst result tree.  Before model allocation it runs the independent
# Julia readiness gate, which requires the complete corrected-SIF publication
# chain and caller-pinned tapered/reference prior hashes.  Use
# GATTACA_TAPERED_PHASE=smoke
# with `sbatch --array=7` first; production refuses to start until that paired
# noiseless state-018 smoke result has passed validation.
# The fixed private `.../slurm` directory must exist before `sbatch`, because
# Slurm opens stdout/stderr before this script can create directories.

#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=0-39%2
#SBATCH --job-name=rrs_sif_tapered
# Account, QOS, and absolute log destinations are intentionally supplied by
# the private host runbook at submission time, not embedded in tracked source.

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

# Safe mapping inspection for the source-level regression test and operator
# review.  This exits before any filesystem, Slurm, or CUDA requirement.
if [[ "${GATTACA_TAPERED_PLAN_ONLY:-0}" == 1 ]]; then
    printf 'task=%s state=%s sif_filter=on prior_model=%s campaign=%s\n' \
        "${task_id}" "${state}" \
        "acos_mapped_tapered_vertical_correlation" \
        "bottom_layer_sif_acos_mapped_tapered_vertical_correlation_v1"
    exit 0
fi

: "${SLURM_JOB_ID:?this launcher must run inside a Slurm allocation}"
: "${CAMPAIGN_CHECKPOINT_SHA:?export the approved full 40-character source SHA}"
: "${TAPERED_PRIOR_SHA256:?export the approved tapered-prior SHA-256}"
: "${CURRENT_PRIOR_SHA256:?export the approved current-prior SHA-256}"
[[ "${CAMPAIGN_CHECKPOINT_SHA}" =~ ^[0-9a-f]{40}$ ]] || {
    echo "CAMPAIGN_CHECKPOINT_SHA must be a lowercase 40-character Git SHA" >&2
    exit 2
}
[[ "${TAPERED_PRIOR_SHA256}" =~ ^[0-9a-f]{64}$ ]] || {
    echo "TAPERED_PRIOR_SHA256 must be a lowercase 64-character SHA-256" >&2
    exit 2
}
[[ "${CURRENT_PRIOR_SHA256}" =~ ^[0-9a-f]{64}$ ]] || {
    echo "CURRENT_PRIOR_SHA256 must be a lowercase 64-character SHA-256" >&2
    exit 2
}

phase="${GATTACA_TAPERED_PHASE:-}"
[[ "${phase}" == smoke || "${phase}" == production ]] || {
    echo "GATTACA_TAPERED_PHASE must be explicitly set to smoke or production" >&2
    exit 2
}
if [[ "${phase}" == smoke && "${task_id}" != 7 ]]; then
    echo "the smoke phase requires only array task 7 (bottom-layer state 018)" >&2
    exit 2
fi

repo_root="${RRS_REPO:-${HOME}/code/uni_vSmartMOM}"
private_root="${RRS_PRIVATE_ROOT:-${HOME}/RRS_XCO2_private}"
input_root="${RRS_INPUTS:-${private_root}/rrs_inputs}"
depot_root="${RRS_DEPOT:-${private_root}/julia_depot_1.12.5}"
restart_root="${SIF_RESTART_ROOT:-${private_root}/results/.sif_v2_restart}"
julia_bin="${JULIA_BIN:-${HOME}/software/julia-1.12.5/bin/julia}"
bottom_campaign="${repo_root}/RRS_XCO2/bottom_layer_XCO2_retrievals"
bottom_archive_root="${private_root}/archive/sif_wavelength_integral_0p5_20260904/bottom_layer_XCO2_retrievals"
truth_root="${bottom_campaign}/truth"
measurement_root="${truth_root}/OCO_radiances"
noise_root="${measurement_root}/noise_covariances"
campaign_id="bottom_layer_sif_acos_mapped_tapered_vertical_correlation_v1"
result_root="${private_root}/results/${campaign_id}"
output_root="${result_root}/retrievals"
prior_path="${result_root}/retrieval_setup/apriori_states_acos_mapped_tapered_vertical_correlation.nc"
current_prior_path="${bottom_campaign}/retrieval_setup/apriori_states.nc"
readiness_script="${repo_root}/RRS_XCO2/inversion/GattacaTaperedSIFReadiness.jl"
runner="${repo_root}/RRS_XCO2/inversion/run_retrievals.jl"
validator="${repo_root}/RRS_XCO2/inversion/validate_bottom_layer_retrievals.jl"

for directory in "${repo_root}" "${private_root}" "${input_root}" \
                 "${depot_root}" "${restart_root}" "${bottom_campaign}" \
                 "${bottom_archive_root}"; do
    [[ -d "${directory}" ]] || {
        echo "missing required directory: ${directory}" >&2
        exit 1
    }
done
for required in "${readiness_script}" "${runner}" "${validator}" \
                "${repo_root}/Manifest.toml" "${prior_path}" \
                "${current_prior_path}" \
                "${private_root}/external_inputs.sha256"; do
    [[ -f "${required}" ]] || {
        echo "missing required tapered-SIF retrieval input: ${required}" >&2
        exit 1
    }
done

repo_root="$(cd "${repo_root}" && pwd -P)"
private_root="$(cd "${private_root}" && pwd -P)"
result_root="${private_root}/results/${campaign_id}"
output_root="${result_root}/retrievals"
case "${output_root}/" in
    "${private_root}/results/${campaign_id}/retrievals/") ;;
    *) echo "internal output-isolation error: ${output_root}" >&2; exit 1 ;;
esac
case "${output_root}/" in
    "${repo_root}/"*)
        echo "refusing to write retrieval products inside the Git checkout" >&2
        exit 1
        ;;
esac

[[ -x "${julia_bin}" ]] || {
    echo "missing Julia 1.12.5 executable: ${julia_bin}" >&2
    exit 1
}
[[ "$("${julia_bin}" --version)" == "julia version 1.12.5" ]] || {
    echo "the tapered-SIF campaign requires exactly Julia 1.12.5" >&2
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
    echo "could not determine scheduler-assigned GPU memory" >&2
    exit 1
}
(( assigned_gpu_memory >= 49152 )) || {
    echo "assigned GPU has ${assigned_gpu_memory} MiB; at least 49152 MiB is required" >&2
    exit 1
}

# These hashes protect the spectroscopy inputs independently of the SIF-chain
# checks performed below.  The transferred manifest is rooted at private_root.
(
    cd "${private_root}"
    sha256sum --check --quiet external_inputs.sha256
)

mkdir -p "${result_root}/slurm" "${result_root}/logs" \
         "${result_root}/state_claims" "${output_root}"

export RRS_REPO="${repo_root}"
export RRS_PRIVATE_ROOT="${private_root}"
export RRS_INPUTS="${input_root}"
export RRS_DEPOT="${depot_root}"
export SIF_RESTART_ROOT="${restart_root}"
export BOTTOM_XCO2_CAMPAIGN_ROOT="${bottom_campaign}"
export BOTTOM_SIF_V1_ARCHIVE_ROOT="${bottom_archive_root}"
export FULL_COLUMN_TRUTH_ROOT="${repo_root}/RRS_XCO2/truth_map"
export RETRIEVAL_TRUTH_TABLE="${truth_root}/true_states.dat"
export RETRIEVAL_MEASUREMENT_DIR="${measurement_root}"
export RETRIEVAL_NOISE_DIR="${noise_root}"
export RETRIEVAL_PRIOR_PATH="${prior_path}"
export CURRENT_PRIOR_PATH="${current_prior_path}"
export RETRIEVAL_OUTPUT_ROOT="${output_root}"
export JULIA_DEPOT_PATH="${depot_root}"
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export JULIA_CPU_TARGET="${JULIA_CPU_TARGET:-native}"
export VSMARTMOM_ABSCO_DIR="${input_root}/absco_v5.2"
export O2_H2O_LUT="${input_root}/hitran/H2O.jld2"
export SOLAR_OUT="${input_root}/solar/solar.out"

readiness_log="${result_root}/logs/readiness_job${SLURM_JOB_ID}_task${task_id}.log"
"${julia_bin}" --project="${repo_root}" --startup-file=no \
    "${readiness_script}" > "${readiness_log}"

identity="${output_root}/.control/campaign_identity.dat"
[[ -f "${identity}" ]] || {
    echo "readiness gate did not publish the campaign identity" >&2
    exit 1
}
identity_sha256="$(sha256sum "${identity}" | awk '{print $1}')"
smoke_receipt="${result_root}/.control/state018_perturbation11_smoke.complete"
if [[ "${phase}" == production ]]; then
    [[ -f "${smoke_receipt}" ]] || {
        echo "production is blocked until the tapered-prior state-018 smoke gate passes" >&2
        exit 1
    }
    grep -qxF "campaign_identity_sha256=${identity_sha256}" "${smoke_receipt}" || {
        echo "smoke receipt does not match the current campaign identity" >&2
        exit 1
    }
fi

claim="${result_root}/state_claims/state${state}.claim"
if ! mkdir "${claim}" 2>/dev/null; then
    echo "state ${state} is already claimed: ${claim}" >&2
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
    printf 'slurm_job_id=%s\nsource_sha=%s\nprior_sha256=%s\n' \
        "${SLURM_JOB_ID}" "${source_sha}" "${TAPERED_PRIOR_SHA256}"
    printf 'campaign_identity_sha256=%s\n' "${identity_sha256}"
} > "${claim}/owner.txt"

provenance="${result_root}/logs/state${state}_job${SLURM_JOB_ID}_${phase}.provenance.txt"
{
    printf 'started_utc=%s\nhost=%s\nstate=%s\nphase=%s\n' \
        "$(date -u +%FT%TZ)" "$(hostname -f)" "${state}" "${phase}"
    printf 'slurm_job_id=%s\nsource_sha=%s\n' \
        "${SLURM_JOB_ID}" "${source_sha}"
    printf 'prior_model=%s\nprior_path=%s\nprior_sha256=%s\n' \
        "acos_mapped_tapered_vertical_correlation" "${prior_path}" \
        "${TAPERED_PRIOR_SHA256}"
    printf 'current_prior_path=%s\ncurrent_prior_sha256=%s\n' \
        "${current_prior_path}" "${CURRENT_PRIOR_SHA256}"
    printf 'output_root=%s\ncampaign_identity_sha256=%s\n' \
        "${output_root}" "${identity_sha256}"
    sha256sum "${private_root}/external_inputs.sha256" \
              "${repo_root}/Manifest.toml" "${identity}"
    nvidia-smi -i "${CUDA_VISIBLE_DEVICES}" \
        --query-gpu=name,memory.total,driver_version --format=csv,noheader
    "${julia_bin}" --version
} > "${provenance}"

first_perturbation=1
last_perturbation=11
write_manifest=0
if [[ "${phase}" == smoke ]]; then
    first_perturbation=11
    last_perturbation=11
    write_manifest=1
fi

# SIF_CASE_FILTER is intentionally literal, not inherited from the submitter.
# This source file therefore cannot launch the later no-SIF tapered-prior run.
env \
    CUDA_DEVICE=0 \
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

validation_mode=products
scene_class=all
expected_perturbations=1-11
if [[ "${phase}" == smoke ]]; then
    validation_mode=smoke
    scene_class=aerosol_all_sif
    expected_perturbations=11
fi
env \
    EXPECTED_STATES="${state#0}" \
    EXPECTED_PERTURBATIONS="${expected_perturbations}" \
    VALIDATION_MODE="${validation_mode}" \
    BOTTOM_RETRIEVAL_SCENE_CLASS="${scene_class}" \
    BOTTOM_RETRIEVAL_CAMPAIGN_ROOT="${bottom_campaign}" \
    RETRIEVAL_OUTPUT_ROOT="${output_root}" \
    RETRIEVAL_PRIOR_PATH="${prior_path}" \
    "${julia_bin}" --project="${repo_root}" --startup-file=no "${validator}"

if [[ "${phase}" == smoke ]]; then
    mkdir -p "$(dirname "${smoke_receipt}")"
    temporary="${smoke_receipt}.tmp.$$"
    {
        printf 'campaign_id=%s\nstate=018\nperturbation=11\n' "${campaign_id}"
        printf 'source_sha=%s\nprior_sha256=%s\n' \
            "${source_sha}" "${TAPERED_PRIOR_SHA256}"
        printf 'campaign_identity_sha256=%s\n' "${identity_sha256}"
        printf 'validated_utc=%s\n' "$(date -u +%FT%TZ)"
    } > "${temporary}"
    mv "${temporary}" "${smoke_receipt}"
fi

printf 'completed_utc=%s\n' "$(date -u +%FT%TZ)" >> "${provenance}"
