#!/usr/bin/env bash
# Round 7_altgrd pilot on Curry physical GPU 1: for every scene, noisy member
# 01 and noiseless member 11 (160 solves). No queue, polling, or GPU takeover.
#   run_round7_altgrd_pilot.sh [--check-only] [--smoke-only]
# --smoke-only: state 001 (no aerosol, SIF off) and state 036 (aerosol, SIF
# on), member 11 only, with strict convergence and fit-quality gates.

r7a_gpu_uuid() {
    local uuid
    uuid="$(nvidia-smi -i 1 --query-gpu=uuid --format=csv,noheader,nounits)" || {
        echo 'STOP: cannot query physical GPU1 through NVML' >&2; return 1; }
    [[ "$uuid" =~ ^GPU-[0-9a-fA-F-]{36}$ ]] || { echo 'STOP: invalid GPU1 UUID' >&2; return 1; }
    printf '%s\n' "$uuid"
}

r7a_require_gpu_idle() {
    local uuid="$1" applications gpu pid extra busy=0
    applications="$(nvidia-smi --query-compute-apps=gpu_uuid,pid --format=csv,noheader,nounits)" || {
        echo 'STOP: NVML process query failed; vacancy is unknown' >&2; return 1; }
    [[ -n "$applications" ]] || return 0
    while IFS=, read -r gpu pid extra; do
        gpu="${gpu//[[:space:]]/}"; pid="${pid//[[:space:]]/}"
        if [[ "$gpu" == "$uuid" ]]; then
            echo "STOP: physical GPU1 is occupied by PID $pid; no queue or launch" >&2
            busy=1
        fi
    done <<< "$applications"
    [[ "$busy" == 0 ]]
}

r7a_selection() {
    local mode="$1" smoke="$2"
    export FIRST_STATE=1 LAST_STATE=80 ROUND7_ALTGRD_PERTURBATIONS=1,11 ROUND7_ALTGRD_SMOKE=0
    export AEROSOL_CASE_FILTER=all
    if [[ "$smoke" == 1 ]]; then
        case "$mode" in off) FIRST_STATE=1 ;; on) FIRST_STATE=36 ;; esac
        export LAST_STATE="$FIRST_STATE" ROUND7_ALTGRD_PERTURBATIONS=11 ROUND7_ALTGRD_SMOKE=1
    fi
}

r7a_main() (
    set -euo pipefail
    umask 022
    local check_only=0 smoke_only=0 option mode gpu_uuid identity_file runner
    for option in "$@"; do
        case "$option" in
            --check-only) check_only=1 ;; --smoke-only) smoke_only=1 ;;
            *) echo "STOP: unknown option: $option" >&2; exit 2 ;;
        esac
    done
    [[ "$(hostname -s)" == curry ]] || { echo 'STOP: the pilot runs on curry only' >&2; exit 1; }
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
    source "$script_dir/round7_altgrd_environment.sh"
    configure_round7_altgrd off
    runner="$script_dir/run_round7_altgrd_retrievals.jl"
    if [[ "$check_only" == 0 ]]; then
        mkdir -p "$ROUND7_ALTGRD_ROOT/.control"
        exec 9>>"$ROUND7_ALTGRD_ROOT/.control/curry1.lock"
        flock -n 9 || { echo 'STOP: the Round 7_altgrd worker already runs' >&2; exit 1; }
        gpu_uuid="$(r7a_gpu_uuid)"
        r7a_require_gpu_idle "$gpu_uuid"
    fi
    # Validate both modes before any solve, without a GPU context.
    for mode in off on; do
        configure_round7_altgrd "$mode"
        r7a_selection "$mode" "$smoke_only"
        ROUND7_ALTGRD_PREFLIGHT_ONLY=1 CUDA_VISIBLE_DEVICES='' \
            "$JULIA_BIN" --project="$ROUND7_ALTGRD_SOURCE_ROOT" --startup-file=no "$runner"
    done
    [[ "$check_only" == 0 ]] || exit 0
    for mode in off on; do
        configure_round7_altgrd "$mode"
        r7a_selection "$mode" "$smoke_only"
        [[ "$(r7a_gpu_uuid)" == "$gpu_uuid" ]] || { echo 'STOP: GPU1 identity changed' >&2; exit 1; }
        r7a_require_gpu_idle "$gpu_uuid"
        identity_file="$ROUND7_ALTGRD_ROOT/.control/${mode}_campaign_identity.sha256"
        if [[ -e "$identity_file" || -L "$identity_file" ]]; then
            [[ -f "$identity_file" && ! -L "$identity_file" && "$(<"$identity_file")" == "$ROUND7_ALTGRD_CAMPAIGN_IDENTITY_SHA256" ]] || {
                echo 'STOP: existing output campaign identity differs' >&2; exit 1; }
        else
            (set -o noclobber; printf '%s\n' "$ROUND7_ALTGRD_CAMPAIGN_IDENTITY_SHA256" > "$identity_file")
        fi
        # Physical GPU 1 only, addressed by UUID; process-local ordinal 0.
        export CUDA_VISIBLE_DEVICES="$gpu_uuid" CUDA_DEVICE=0
        echo "started_utc=$(date -u +%FT%TZ) worker=curry1 physical_gpu=1 uuid=$gpu_uuid SIF=$mode smoke=$smoke_only"
        "$JULIA_BIN" --project="$ROUND7_ALTGRD_SOURCE_ROOT" --startup-file=no "$runner"
        echo "finished_utc=$(date -u +%FT%TZ) worker=curry1 SIF=$mode"
    done
)

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    r7a_main "$@"
fi
