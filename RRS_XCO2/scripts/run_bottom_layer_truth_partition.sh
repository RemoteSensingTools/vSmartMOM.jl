#!/usr/bin/env bash

set -euo pipefail

usage() {
    echo "Usage: $0 {none|aerosol} {nosif|sif}" >&2
    echo "  none    : Curry partition (20 clear scenes)" >&2
    echo "  aerosol : Wurst partition (20 aerosol scenes)" >&2
    echo "  nosif must finish on both hosts before either sif stage is started." >&2
}

partition="${1:-}"
stage="${2:-}"
if [[ "$partition" != "none" && "$partition" != "aerosol" ]]; then
    usage
    exit 2
fi
if [[ "$stage" != "nosif" && "$stage" != "sif" ]]; then
    usage
    exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"
campaign_root="${BOTTOM_XCO2_CAMPAIGN_ROOT:-$project_root/RRS_XCO2/bottom_layer_XCO2_retrievals}"
truth_root="$campaign_root/truth"
log_root="$campaign_root/logs"
host="$(hostname -s)"
physical_gpu="${BOTTOM_TRUTH_PHYSICAL_GPU:-1}"
if [[ ! "$physical_gpu" =~ ^[0-9]+$ ]]; then
    echo "BOTTOM_TRUTH_PHYSICAL_GPU must be a non-negative physical GPU index." >&2
    exit 2
fi
expected_host="curry"
if [[ "$partition" == "aerosol" ]]; then
    expected_host="wurst"
fi

if [[ "$host" != *"$expected_host"* && "${ALLOW_BOTTOM_TRUTH_HOST_OVERRIDE:-0}" != "1" ]]; then
    echo "Refusing $partition production on $host; this partition is assigned to $expected_host." >&2
    echo "Set ALLOW_BOTTOM_TRUTH_HOST_OVERRIDE=1 only after explicitly reassigning the work." >&2
    exit 1
fi

wait_for_idle="${WAIT_FOR_RETRIEVALS:-1}"
poll_seconds="${IDLE_POLL_SECONDS:-600}"
retrieval_pattern='[j]ulia.*RRS_XCO2/inversion/run_retrievals[.]jl'

# Curry keeps a small, permanent PostgreSQL CUDA context on GPU 0. It is
# infrastructure rather than a compute workload. Ignore only that explicitly
# named keeper; every other process on the selected GPU remains blocking.
selected_gpu_busy_apps() {
    local selected_uuid="$1"
    local all_apps="$2"
    local line
    while IFS= read -r line; do
        [[ "$line" == "$selected_uuid,"* ]] || continue
        if [[ "$line" =~ postgres:\ GPU[0-9]+\ memory\ keeper ]]; then
            continue
        fi
        printf '%s\n' "$line"
    done <<< "$all_apps"
    return 0
}

if [[ "${BOTTOM_TRUTH_LAUNCH_DRY_RUN:-0}" == "1" ]]; then
    retrieval_status="idle"
    pgrep -af "$retrieval_pattern" >/dev/null && retrieval_status="retrieval-active"
    gpu_status="not-needed"
    if [[ "$stage" == "nosif" ]]; then
        gpu_status="unavailable"
        if gpu_uuid="$(nvidia-smi -i "$physical_gpu" --query-gpu=uuid --format=csv,noheader 2>/dev/null)"; then
            if gpu_apps="$(nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name --format=csv,noheader 2>/dev/null)"; then
                gpu_busy_apps="$(selected_gpu_busy_apps "$gpu_uuid" "$gpu_apps")"
                gpu_status="idle"
                [[ -n "$gpu_busy_apps" ]] && gpu_status="compute-active"
            else
                gpu_status="query-failed"
            fi
        fi
    fi
    echo "Validated launch: host=$host partition=$partition stage=$stage physical_gpu=$physical_gpu process_status=$retrieval_status gpu_status=$gpu_status"
    echo "Project: $project_root"
    echo "Campaign: $campaign_root"
    exit 0
fi
if [[ "${BOTTOM_TRUTH_ALLOW_RETRIEVAL_OVERLAP:-0}" != "1" ]]; then
    if [[ "$wait_for_idle" == "1" ]]; then
        while pgrep -af "$retrieval_pattern" >/dev/null; do
            echo "$(date -u +%FT%TZ) waiting: a retrieval worker is still active on $host"
            sleep "$poll_seconds"
        done
    elif pgrep -af "$retrieval_pattern" >/dev/null; then
        echo "Refusing to overlap truth production with an active retrieval worker on $host." >&2
        exit 1
    fi
elif pgrep -af "$retrieval_pattern" >/dev/null; then
    echo "$(date -u +%FT%TZ) allowing the active retrieval worker because BOTTOM_TRUTH_ALLOW_RETRIEVAL_OVERLAP=1; selected-GPU isolation is still enforced"
fi

if [[ "$stage" == "nosif" ]]; then
    if ! gpu_uuid="$(nvidia-smi -i "$physical_gpu" --query-gpu=uuid --format=csv,noheader 2>/dev/null)"; then
        echo "Cannot verify physical GPU $physical_gpu; refusing to launch production." >&2
        exit 1
    fi
    while true; do
        if ! gpu_apps="$(nvidia-smi --query-compute-apps=gpu_uuid,pid,process_name --format=csv,noheader 2>/dev/null)"; then
            echo "Cannot query active GPU processes; refusing to assume GPU $physical_gpu is idle." >&2
            exit 1
        fi
        gpu_busy_apps="$(selected_gpu_busy_apps "$gpu_uuid" "$gpu_apps")"
        if [[ -z "$gpu_busy_apps" ]]; then
            break
        fi
        if [[ "$wait_for_idle" != "1" ]]; then
            echo "Refusing to use physical GPU $physical_gpu while another compute process is active." >&2
            exit 1
        fi
        echo "$(date -u +%FT%TZ) waiting: physical GPU $physical_gpu still has a non-infrastructure compute process"
        printf '%s\n' "$gpu_busy_apps"
        sleep "$poll_seconds"
    done
fi

# Hide every device except the explicitly selected physical GPU so Julia
# addresses that device as its local device 0 unambiguously.
export CUDA_VISIBLE_DEVICES="$physical_gpu"
export CUDA_DEVICE="0"
export AEROSOL_CASE_FILTER="$partition"
export BOTTOM_TRUTH_STAGE="$stage"
export BOTTOM_XCO2_CAMPAIGN_ROOT="$campaign_root"
export RRS_XCO2_FLOAT_TYPE="Float32"
export NLAYERS="16"
export AEROSOL_NSTREAMS="9"
export TRUTH_SZA_DEG="30"
export TRUTH_VZA_DEG="0"
export TRUTH_RELATIVE_AZIMUTH_DEG="0"
if [[ -z "${JULIA_DEPOT_PATH:-}" ]]; then
    if [[ "$partition" == "aerosol" ]]; then
        export JULIA_DEPOT_PATH="$project_root/RRS_XCO2/.julia_depot_wurst:/home/sanghavi/.julia"
    else
        export JULIA_DEPOT_PATH="$project_root/RRS_XCO2/.julia_depot_curry:/home/sanghavi/.julia"
    fi
fi
export JULIA_CPU_TARGET="${JULIA_CPU_TARGET:-native}"

mkdir -p "$truth_root/.claims" "$log_root"
launcher_claim="$truth_root/.claims/launcher_${partition}_${stage}.claim"
if ! mkdir "$launcher_claim"; then
    echo "Partition already claimed: $launcher_claim" >&2
    [[ -f "$launcher_claim/owner.txt" ]] && sed -n '1,20p' "$launcher_claim/owner.txt" >&2
    exit 1
fi

cleanup_launcher_claim() {
    rm -f "$launcher_claim/owner.txt"
    rmdir "$launcher_claim" 2>/dev/null || true
}
trap cleanup_launcher_claim EXIT INT TERM

{
    echo "host=$host"
    echo "pid=$$"
    echo "started_utc=$(date -u +%FT%TZ)"
    echo "partition=$partition"
    echo "stage=$stage"
    echo "physical_gpu=$physical_gpu"
    echo "physical_gpu_uuid=${gpu_uuid:-not-needed}"
    echo "logical_cuda_device=0"
    echo "retrieval_overlap_opt_in=${BOTTOM_TRUTH_ALLOW_RETRIEVAL_OVERLAP:-0}"
    echo "ignored_gpu_infrastructure=postgres GPU[0-9]+ memory keeper"
} > "$launcher_claim/owner.txt"

cd "$project_root"
julia --project=. RRS_XCO2/scripts/generate_bottom_layer_truth_states.jl

log_file="$log_root/${host}_gpu${physical_gpu}_${partition}_${stage}_$(date -u +%Y%m%dT%H%M%SZ).log"
echo "Writing production log to $log_file"
julia --project=. RRS_XCO2/scripts/generate_bottom_layer_truth.jl 2>&1 | tee "$log_file"
