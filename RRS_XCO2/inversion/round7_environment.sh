#!/usr/bin/env bash
# Source only from a sealed release. configure_round7 off|on [local].

configure_round7() {
    local mode="${1:?SIF mode required}" site="${2:-local}"
    local bundle_dir data_checkout assignments key value count prior_name
    local -A identity_keys=()
    [[ "$mode" == off || "$mode" == on ]] || { echo 'STOP: mode must be off or on' >&2; return 1; }
    [[ "$site" == local ]] || { echo 'STOP: only Curry/Wurst deployment is authorized' >&2; return 1; }
    bundle_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)" || return 1
    data_checkout="$(realpath -e "${ROUND7_DATA_CHECKOUT:-$HOME/code/github/uni_vSmartMOM}")" || return 1
    export ROUND7_SOURCE_ROOT
    ROUND7_SOURCE_ROOT="$(realpath -e "${ROUND7_SOURCE_ROOT:-$HOME/code/github/uni_vSmartMOM_round5_jacobian}")" || return 1
    export ROUND7_ROOT
    ROUND7_ROOT="$(realpath -e "${ROUND7_ROOT:-$data_checkout/RRS_XCO2/bottom_layer_XCO2_retrievals/round7_fixed_sif_imperfect_correction}")" || return 1
    export ROUND7_BUNDLE="$bundle_dir" ROUND7_OBSERVATION_DIR="$ROUND7_ROOT/observations"
    export JULIA_BIN="${JULIA_BIN:-$HOME/.julia/juliaup/julia-1.12.5+0.x64.linux.gnu/bin/julia}"
    [[ -x "$JULIA_BIN" ]] || { echo "STOP: missing Julia: $JULIA_BIN" >&2; return 1; }
    [[ "$("$JULIA_BIN" --startup-file=no --version)" == 'julia version 1.12.5' ]] || {
        echo 'STOP: Round7 requires Julia 1.12.5' >&2; return 1;
    }
    export VSMARTMOM_ABSCO_DIR="${VSMARTMOM_ABSCO_DIR:-/net/fluo/data1/ABSCO_CS_Database/v5.2_final}"
    export O2_H2O_LUT="${O2_H2O_LUT:-$HOME/data/HITRAN_LUTs/H2O.jld2}"
    export SOLAR_OUT="${SOLAR_OUT:-$HOME/Raman_misc/workflows/worktrees/uni_vSmartMOM_sanghavi_2025-03-18/src/SolarModel/solar.out}"
    # Override inherited individual-table settings; all five ABSCO tables come
    # from one explicitly selected database, as in the Round-5 production run.
    export O2_ABSCO_LUT="$VSMARTMOM_ABSCO_DIR/o2_v52_v2.jld2"
    export WEAK_H2O_ABSCO_LUT="$VSMARTMOM_ABSCO_DIR/wh2o_v52.jld2"
    export WEAK_CO2_ABSCO_LUT="$VSMARTMOM_ABSCO_DIR/wco2_v52.jld2"
    export STRONG_H2O_ABSCO_LUT="$VSMARTMOM_ABSCO_DIR/sh2o_v52.jld2"
    export STRONG_CO2_ABSCO_LUT="$VSMARTMOM_ABSCO_DIR/sco2_v52.jld2"
    for value in "$O2_ABSCO_LUT" "$O2_H2O_LUT" "$WEAK_H2O_ABSCO_LUT" "$WEAK_CO2_ABSCO_LUT" "$STRONG_H2O_ABSCO_LUT" "$STRONG_CO2_ABSCO_LUT"; do
        [[ -f "$value" && -r "$value" ]] || { echo "STOP: missing spectroscopy table: $value" >&2; return 1; }
    done
    export RRS_XCO2_DATA_ROOT="$data_checkout/RRS_XCO2"
    export RRS_XCO2_CONFIG="$ROUND7_SOURCE_ROOT/sandbox/workflows/RRS_XCO2/config/oco_grass_3aerosol.yaml"
    export RETRIEVAL_TRUTH_TABLE="$bundle_dir/inputs/true_states_corrected_sif_v2.dat"
    prior_name="apriori_states_round5_fixed_sif_${mode}_tight_utls_acos_mapped_tapered_vertical_correlation.nc"
    export RETRIEVAL_PRIOR_PATH="$bundle_dir/inputs/$prior_name"
    export ROUND7_SOURCE_PRIOR="$bundle_dir/inputs/source_prior.nc"
    export RETRIEVAL_STOKES_COEFFICIENT_PATH="$bundle_dir/inputs/representative_stokes_coefficients.nc"
    export RETRIEVAL_SCENE_COMPONENTS_PATH="$bundle_dir/inputs/scene_components.dat"
    export RRS_XCO2_SIF_TEMPLATE_PATH="$bundle_dir/inputs/sif-spectra.csv"
    # The new runner consumes only ROUND7_OBSERVATION_DIR. Never inherit the
    # legacy SIF measurement/noise directories from a previous shell session.
    unset RETRIEVAL_MEASUREMENT_DIR RETRIEVAL_NOISE_DIR FULL_COLUMN_TRUTH_ROOT
    unset ROUND7_PRINT_IDENTITY ROUND7_PREFLIGHT_ONLY
    export PYTHONDONTWRITEBYTECODE=1
    assignments="$("${ROUND7_PYTHON:-python3}" "$bundle_dir/package_round7_release.py" verify \
        --bundle "$bundle_dir" --source-root "$ROUND7_SOURCE_ROOT" --round7-root "$ROUND7_ROOT" \
        --solar-out "$SOLAR_OUT" --mode "$mode")" || return 1
    count=0
    while IFS='=' read -r key value; do
        case "$key" in
            ROUND7_CODE_CHECKPOINT_SHA) [[ "$value" =~ ^[0-9a-f]{40}$ ]] || return 1 ;;
            ROUND7_CODESET_SHA256|ROUND7_INPUT_SET_SHA256|ROUND7_CAMPAIGN_IDENTITY_SHA256|ROUND7_OBSERVATION_MANIFEST_SHA256)
                [[ "$value" =~ ^[0-9a-f]{64}$ ]] || return 1 ;;
            *) echo "STOP: unexpected identity key: $key" >&2; return 1 ;;
        esac
        [[ -z "${identity_keys[$key]:-}" ]] || { echo "STOP: duplicate identity key: $key" >&2; return 1; }
        identity_keys[$key]=1
        printf -v "$key" '%s' "$value"
        export "$key"
        count=$((count + 1))
    done <<< "$assignments"
    [[ "$count" == 5 ]] || { echo 'STOP: incomplete identity' >&2; return 1; }
    export OBSERVATION_MANIFEST_SHA256="$ROUND7_OBSERVATION_MANIFEST_SHA256"
    export SIF_CASE_FILTER="$mode" ROUND7_SIF_MODE="$mode" ROUND7_KNOWN_WAVELENGTH_NM=759
    export RETRIEVAL_OUTPUT_ROOT="$ROUND7_ROOT/retrievals_nosif"
    [[ "$mode" != on ]] || export RETRIEVAL_OUTPUT_ROOT="$ROUND7_ROOT/retrievals_sif"
    export RETRIEVAL_CLASS=corrected RETRIEVAL_ARCH=GPU RETRIEVAL_FLOAT_TYPE=Float32 RETRIEVAL_NSTREAMS=9
    export RETRIEVAL_WRITE_MANIFEST=0 FORCE=0 FAIL_FAST=1
    export JULIA_CPU_TARGET=generic JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
    export JULIA_LOAD_PATH='@:@stdlib' JULIA_PKG_PRECOMPILE_AUTO=0
    export VSMARTMOM_CUDA_BATCH_INV_PREALLOC=0
    echo "Round7 SIF=$mode identity=$ROUND7_CAMPAIGN_IDENTITY_SHA256 observations=$ROUND7_OBSERVATION_DIR" >&2
}
