#!/usr/bin/env bash
# No queue, polling, background launch, or GPU takeover. One host owns both modes.

round7_gpu_uuid() {
    local uuid
    uuid="$(nvidia-smi -i 1 --query-gpu=uuid --format=csv,noheader,nounits)" || {
        echo 'STOP: cannot query physical GPU1 through NVML' >&2; return 1;
    }
    [[ "$uuid" =~ ^GPU-[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$ ]] || {
        echo 'STOP: ambiguous/invalid GPU1 UUID' >&2; return 1;
    }
    printf '%s\n' "$uuid"
}

round7_require_gpu_idle() {
    local uuid="$1" applications gpu pid extra busy=0
    applications="$(nvidia-smi --query-compute-apps=gpu_uuid,pid --format=csv,noheader,nounits)" || {
        echo 'STOP: NVML process query failed; vacancy is unknown' >&2; return 1;
    }
    [[ -n "$applications" ]] || return 0
    while IFS=, read -r gpu pid extra; do
        gpu="${gpu//[[:space:]]/}"; pid="${pid//[[:space:]]/}"
        [[ "$gpu" =~ ^GPU-[0-9a-fA-F-]{36}$ && "$pid" =~ ^[0-9]+$ && -z "$extra" ]] || {
            echo 'STOP: malformed/unsupported NVML process output; vacancy is unknown' >&2; return 1;
        }
        if [[ "$gpu" == "$uuid" ]]; then
            echo "STOP: physical GPU1 is occupied by PID $pid; no queue or launch" >&2
            busy=1
        fi
    done <<< "$applications"
    [[ "$busy" == 0 ]]
}

round7_selection() {
    local mode="$1" worker="$2" smoke="$3"
    export FIRST_STATE=1 LAST_STATE=80 FIRST_PERTURBATION=1 LAST_PERTURBATION=11 ROUND7_SMOKE=0
    if [[ "$smoke" == 1 ]]; then
        case "$worker:$mode" in
            curry1:off) FIRST_STATE=1 ;; curry1:on) FIRST_STATE=6 ;;
            wurst1:off) FIRST_STATE=11 ;; wurst1:on) FIRST_STATE=16 ;;
            *) return 1 ;;
        esac
        export LAST_STATE="$FIRST_STATE" FIRST_PERTURBATION=11 LAST_PERTURBATION=11 ROUND7_SMOKE=1
    fi
}

round7_main() (
    set -euo pipefail
    umask 077
    worker="${1:?usage: run_round7_partition.sh curry1|wurst1 [--check-only] [--smoke-only]}"
    shift
    check_only="${ROUND7_CHECK_ONLY:-0}"
    smoke_only="${ROUND7_SMOKE_ONLY:-${ROUND7_SMOKE:-0}}"
    for option in "$@"; do
        case "$option" in
            --check-only) check_only=1 ;; --smoke-only) smoke_only=1 ;;
            *) echo "STOP: unknown option: $option" >&2; exit 2 ;;
        esac
    done
    [[ "$check_only" =~ ^[01]$ && "$smoke_only" =~ ^[01]$ ]] || exit 2
    [[ "${ROUND7_AUTO_QUEUE:-0}" == 0 && "${ROUND7_WAIT_FOR_GPU:-0}" == 0 ]] || {
        echo 'STOP: automatic queue/wait is not authorized' >&2; exit 2;
    }
    case "$worker" in
        curry1) expected_host=curry; aerosol=none ;;
        wurst1) expected_host=wurst; aerosol=aerosol ;;
        *) echo 'STOP: only curry1 and wurst1 are authorized' >&2; exit 2 ;;
    esac
    [[ "$(hostname -s)" == "$expected_host" ]] || { echo 'STOP: worker/host mismatch' >&2; exit 1; }
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
    source "$script_dir/round7_environment.sh"
    configure_round7 off local
    export JULIA_DEPOT_PATH="${ROUND7_DEPOT:-$RRS_XCO2_DATA_ROOT/.julia_depot_${expected_host}:$HOME/.julia:}"
    export AEROSOL_CASE_FILTER="$aerosol" CUDA_DEVICE=0
    runner="$script_dir/inversion/run_round7_imperfect_correction_retrievals.jl"
    if [[ "$check_only" == 0 ]]; then
        mkdir -p "$ROUND7_ROOT/.control"
        # Append mode preserves lock files; never unlink a lock inode. Hold FD9
        # across preflight and both Julia children, including their descendants.
        exec 9>>"$ROUND7_ROOT/.control/$worker.lock"
        flock -n 9 || { echo "STOP: $worker already owns this partition" >&2; exit 1; }
        gpu_uuid="$(round7_gpu_uuid)"
        round7_require_gpu_idle "$gpu_uuid"
    fi
    # Validate BOTH modes before any solve. Check-only creates no locks, outputs,
    # identity files or GPU contexts and can run while another user owns GPU1.
    for mode in off on; do
        configure_round7 "$mode" local
        round7_selection "$mode" "$worker" "$smoke_only"
        ROUND7_PREFLIGHT_ONLY=1 CUDA_VISIBLE_DEVICES='' \
            "$JULIA_BIN" --project="$ROUND7_SOURCE_ROOT" --startup-file=no "$runner"
    done
    [[ "$check_only" == 0 ]] || exit 0
    for mode in off on; do
        configure_round7 "$mode" local
        round7_selection "$mode" "$worker" "$smoke_only"
        # Recheck before each child. The shared lock protects our workers; NVML
        # protects against already-running jobs outside this campaign.
        [[ "$(round7_gpu_uuid)" == "$gpu_uuid" ]] || { echo 'STOP: GPU1 identity changed' >&2; exit 1; }
        round7_require_gpu_idle "$gpu_uuid"
        identity_file="$ROUND7_ROOT/.control/${mode}_campaign_identity.sha256"
        exec 8>>"$ROUND7_ROOT/.control/identity.lock"
        flock -n 8 || { echo 'STOP: another worker is sealing campaign identity; retry explicitly' >&2; exit 1; }
        if [[ -e "$identity_file" || -L "$identity_file" ]]; then
            [[ -f "$identity_file" && ! -L "$identity_file" && "$(<"$identity_file")" == "$ROUND7_CAMPAIGN_IDENTITY_SHA256" ]] || {
                echo 'STOP: existing output campaign identity differs' >&2; exit 1;
            }
        else
            (set -o noclobber; printf '%s\n' "$ROUND7_CAMPAIGN_IDENTITY_SHA256" > "$identity_file")
        fi
        flock -u 8
        exec 8>&-
        export CUDA_VISIBLE_DEVICES="$gpu_uuid" CUDA_DEVICE=0
        echo "started_utc=$(date -u +%FT%TZ) worker=$worker physical_gpu=1 uuid=$gpu_uuid SIF=$mode smoke=$smoke_only"
        "$JULIA_BIN" --project="$ROUND7_SOURCE_ROOT" --startup-file=no "$runner"
        echo "finished_utc=$(date -u +%FT%TZ) worker=$worker SIF=$mode"
    done
)

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    round7_main "$@"
fi
