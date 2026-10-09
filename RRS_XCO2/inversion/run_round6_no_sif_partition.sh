#!/usr/bin/env bash
set -euo pipefail
umask 077
worker="${1:?usage: run_round6_no_sif_partition.sh curry1|wurst1}"
[[ $# == 1 ]] || exit 2
case "$worker" in
    curry1) expected_host=curry; aerosol=none ;;
    wurst1) expected_host=wurst; aerosol=aerosol ;;
    *) echo "Only curry1 (no aerosol) and wurst1 (with aerosol) are authorized" >&2; exit 2 ;;
esac
[[ "$(hostname -s)" == "$expected_host" ]]
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "$script_dir/round6_environment.sh"
export JULIA_DEPOT_PATH="${ROUND6_DEPOT:-$HOME/code/github/uni_vSmartMOM/RRS_XCO2/.julia_depot_${expected_host}:$HOME/.julia:}"
configure_round6 off local
export RETRIEVAL_OUTPUT_ROOT="${ROUND6_OUTPUT_ROOT:-$RRS_XCO2_DATA_ROOT/bottom_layer_XCO2_retrievals/round6_fixed_sif/retrievals_nosif}"
export AEROSOL_CASE_FILTER="$aerosol" FIRST_STATE=1 LAST_STATE=80 FIRST_PERTURBATION=1 LAST_PERTURBATION=11
export ROUND6_INITIAL_PAIR_GATE=1
if [[ "${ROUND6_CHECK_ONLY:-0}" == 1 ]]; then
    export ROUND6_PREFLIGHT_ONLY=1
    "$JULIA_BIN" --project="$ROUND6_SOURCE_ROOT" --startup-file=no "$script_dir/inversion/run_round6_fixed_sif_retrievals.jl"
    exit
fi
mkdir -p "$RETRIEVAL_OUTPUT_ROOT/.control"
exec 8>"$RETRIEVAL_OUTPUT_ROOT/.control/identity.lock"
flock 8
identity_file="$RETRIEVAL_OUTPUT_ROOT/.control/campaign_identity.sha256"
if [[ -f "$identity_file" ]]; then
    [[ "$(<"$identity_file")" == "$ROUND6_CAMPAIGN_IDENTITY_SHA256" ]] || {
        echo "STOP: existing Round-6 outputs have a different campaign identity" >&2; exit 1;
    }
else
    printf '%s\n' "$ROUND6_CAMPAIGN_IDENTITY_SHA256" > "$identity_file"
fi
flock -u 8
exec 9>"$RETRIEVAL_OUTPUT_ROOT/.control/$worker.lock"
flock -n 9 || { echo "$worker already owns this Round-6 partition" >&2; exit 1; }
gpu_uuid="$(nvidia-smi -i 1 --query-gpu=uuid --format=csv,noheader)"
applications="$(nvidia-smi --query-compute-apps=gpu_uuid,pid --format=csv,noheader)"
if printf '%s\n' "$applications" | awk -F, -v gpu="$gpu_uuid" '$1==gpu {found=1} END {exit !found}'; then
    echo "STOP: physical GPU 1 already has a compute process" >&2
    exit 1
fi
# Restrict visibility to physical device 1; its process-local CUDA ordinal is 0.
export CUDA_VISIBLE_DEVICES="$gpu_uuid" CUDA_DEVICE=0
echo "worker=$worker physical_gpu=1 uuid=$gpu_uuid aerosol_filter=$aerosol"
echo "started_utc=$(date -u +%FT%TZ) output=$RETRIEVAL_OUTPUT_ROOT"
"$JULIA_BIN" --project="$ROUND6_SOURCE_ROOT" --startup-file=no "$script_dir/inversion/run_round6_fixed_sif_retrievals.jl"
echo "finished_utc=$(date -u +%FT%TZ) worker=$worker"
