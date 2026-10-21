#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

export CUDA_DEVICE="${CUDA_DEVICE:-1}"
export NREPS="${NREPS:-2}"
export RUN_LINEAR="${RUN_LINEAR:-1}"
export RRS_XCO2_OUTPUT_DIR="${RRS_XCO2_OUTPUT_DIR:-${REPO_ROOT}/RRS_XCO2/output/wurst_gpu1}"

cd "${REPO_ROOT}"
mkdir -p "${RRS_XCO2_OUTPUT_DIR}"
julia --project=. RRS_XCO2/scripts/benchmark_basis_gpu.jl \
  2>&1 | tee "${RRS_XCO2_OUTPUT_DIR}/basis_benchmark_gpu1.log"
python3 RRS_XCO2/scripts/plot_basis_results.py
