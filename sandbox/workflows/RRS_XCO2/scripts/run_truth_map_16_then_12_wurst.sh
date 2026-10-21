#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TRUTH_ROOT="${REPO_ROOT}/RRS_XCO2/truth_map"
LAYERS12_OUT="${TRUTH_ROOT}/sims_12layer"

export CUDA_DEVICE="${CUDA_DEVICE:-1}"

cd "${REPO_ROOT}"
mkdir -p "${LAYERS12_OUT}"

echo "$(date --iso-8601=seconds) starting corrected 16-layer truth map on CUDA device ${CUDA_DEVICE}"
NLAYERS=16 FORCE=1 TRUTH_OUT="${TRUTH_ROOT}" \
    julia --project=. RRS_XCO2/scripts/generate_truth_map.jl

echo "$(date --iso-8601=seconds) completed 16-layer truth map; starting 12-layer truth map"
NLAYERS=12 FORCE=1 TRUTH_OUT="${LAYERS12_OUT}" \
    julia --project=. RRS_XCO2/scripts/generate_truth_map.jl

echo "$(date --iso-8601=seconds) completed corrected 12-layer truth map"
