#!/usr/bin/env bash
#SBATCH --account=cmml
#SBATCH --qos=normal
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=rrs_r6_sif
set -euo pipefail
umask 077
bundle="${ROUND6_BUNDLE:?submit with submit_round6_sif_on_gattaca.sh}"
source "$bundle/round6_environment.sh"
configure_round6 on gattaca
SIF_STATES=(6 7 8 9 10 16 17 18 19 20 26 27 28 29 30 36 37 38 39 40
            46 47 48 49 50 56 57 58 59 60 66 67 68 69 70 76 77 78 79 80)
task="${SLURM_ARRAY_TASK_ID:?missing array task index}"
[[ "$task" =~ ^[0-9]+$ ]] && (( task >= 0 && task < 40 )) || exit 2
state="${SIF_STATES[task]}"
phase="${ROUND6_PHASE:?missing phase}"
[[ "$phase" == smoke || "$phase" == production ]] || exit 2
export RETRIEVAL_OUTPUT_ROOT="${ROUND6_PRIVATE_ROOT:-$HOME/RRS_XCO2_private}/results/bottom_layer_round6_fixed_sif_on_standard_utls_acos_mapped_tapered_vertical_correlation_v1/retrievals"
export FIRST_STATE="$state" LAST_STATE="$state" AEROSOL_CASE_FILTER=all
export FIRST_PERTURBATION=1 LAST_PERTURBATION=11
if [[ "$phase" == smoke ]]; then
    [[ "$task" == 7 ]] || exit 2
    export FIRST_PERTURBATION=11 ROUND6_SMOKE=1
fi
[[ -n "${CUDA_VISIBLE_DEVICES:-}" && "$CUDA_VISIBLE_DEVICES" != *,* ]]
export CUDA_DEVICE=0
claim_root="$(dirname "$RETRIEVAL_OUTPUT_ROOT")/state_claims"
mkdir -p "$claim_root"
claim="$claim_root/state${state}.claim"
mkdir "$claim" || { echo "state $state is already owned by another task" >&2; exit 1; }
trap 'rmdir "$claim"' EXIT
"$JULIA_BIN" --project="$ROUND6_SOURCE_ROOT" --startup-file=no "$bundle/inversion/run_round6_fixed_sif_retrievals.jl"
