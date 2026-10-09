#!/usr/bin/env bash
set -euo pipefail
umask 077
bundle="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "$bundle/round6_environment.sh"
configure_round6 on gattaca
result_root="${ROUND6_PRIVATE_ROOT:-$HOME/RRS_XCO2_private}/results/bottom_layer_round6_fixed_sif_on_standard_utls_acos_mapped_tapered_vertical_correlation_v1"
if [[ "${ROUND6_CHECK_ONLY:-0}" == 1 ]]; then
    echo "Round-6 static checks passed. No jobs submitted."
    exit
fi
mkdir -p "$result_root/slurm" "$result_root/retrieval_setup"
exec 9>"$result_root/retrieval_setup/submission.lock"
flock -n 9 || { echo "Another Round-6 submission is active" >&2; exit 1; }
record="$result_root/retrieval_setup/retrieval_jobs.env"
[[ ! -e "$record" ]] || { echo "STOP: submission already recorded at $record" >&2; exit 1; }
export ROUND6_BUNDLE="$bundle"
common=(--parsable --output="$result_root/slurm/round6-sif-%A_%a.out"
        --error="$result_root/slurm/round6-sif-%A_%a.err")
smoke_job="$(sbatch "${common[@]}" --array=7 --export=ALL,ROUND6_PHASE=smoke "$bundle/round6_sif_on_gattaca.sbatch")"
smoke_job="${smoke_job%%;*}"
printf 'SMOKE_JOB=%s\nCAMPAIGN_IDENTITY_SHA256=%s\n' "$smoke_job" "$ROUND6_CAMPAIGN_IDENTITY_SHA256" > "$record"
production_job="$(sbatch "${common[@]}" --array=0-39%2 --dependency="afterok:$smoke_job" \
    --export=ALL,ROUND6_PHASE=production "$bundle/round6_sif_on_gattaca.sbatch")"
production_job="${production_job%%;*}"
printf 'PRODUCTION_JOB=%s\nRESULT_ROOT=%q\n' "$production_job" "$result_root" >> "$record"
echo "Smoke (state 018, noiseless, both classes): $smoke_job"
echo "Dependent full 40-state SIF-on array (up to 2 GPUs): $production_job"
echo "Results: $result_root/retrievals"
echo "Logs: $result_root/slurm"
echo "Submission record: $record"
