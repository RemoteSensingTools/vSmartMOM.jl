#!/usr/bin/env bash
# Round 7_altgrd environment. Source, then:
#   seal_round7_altgrd_release <commit>   # once: pins the clean checkout
#   configure_round7_altgrd off|on        # every launch: verify and export
#
# Every input is the sealed Round 7 release, read-only: its observations,
# priors, truth table, analyzer, aerosol components, SIF template, Julia
# Manifest and depot. Only the code checkout (this branch at the pinned
# commit) and the new output tree differ.

_r7a_fail() { echo "STOP: $*" >&2; return 1; }

_r7a_paths() {
    local script_dir
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)" || return 1
    export ROUND7_ALTGRD_SOURCE_ROOT
    ROUND7_ALTGRD_SOURCE_ROOT="$(realpath -e "${ROUND7_ALTGRD_SOURCE_ROOT:-$script_dir/../../../..}")" || return 1
    local data_checkout
    data_checkout="$(realpath -e "${ROUND7_DATA_CHECKOUT:-$HOME/code/github/uni_vSmartMOM}")" || return 1
    export RRS_XCO2_DATA_ROOT="$data_checkout/RRS_XCO2"
    local bottom="$RRS_XCO2_DATA_ROOT/bottom_layer_XCO2_retrievals"
    export ROUND7_REFERENCE_ROOT="$bottom/round7_fixed_sif_imperfect_correction"
    export ROUND7_BUNDLE="$ROUND7_REFERENCE_ROOT/deployment/round7_fixed_sif_imperfect_correction_v1"
    export ROUND7_ALTGRD_ROOT="$bottom/round7_altgrd_fixed_sif_imperfect_correction"
    export ROUND7_ALTGRD_RELEASE="$ROUND7_ALTGRD_ROOT/deployment/round7_altgrd_release.json"
    export JULIA_BIN="${JULIA_BIN:-$HOME/.julia/juliaup/julia-1.12.5+0.x64.linux.gnu/bin/julia}"
}

# Values recorded by the completed Round 7 deployment (release.json, logs).
ROUND7_RELEASE_SHA256SUMS=7014e80434d8a55761b9893fdc96a9098e31745a6fae3423ba4526db1915da79
ROUND7_OBSERVATION_SHA256SUMS=a281bedd941037177a1958fa1455ce3e5f53674ccc1854a20796dfd2aafe8f47
ROUND7_SOLAR_SHA256=9e44e424788ce9f91c654398a789cd9da80205b9d6dc0c00df5e20f500ee8644
ROUND7_CONFIG_SHA256=30fac998ea635ed0bd994bab7b89d173042b24c715ba0e8e8b22e96e38647c81

_r7a_sha() { sha256sum "$1" | awk '{print $1}'; }

# The six spectroscopy tables used by Round 7 (same defaults as round7_environment.sh).
_r7a_spectroscopy() {
    export VSMARTMOM_ABSCO_DIR="${VSMARTMOM_ABSCO_DIR:-/net/fluo/data1/ABSCO_CS_Database/v5.2_final}"
    export O2_H2O_LUT="${O2_H2O_LUT:-$HOME/data/HITRAN_LUTs/H2O.jld2}"
    export O2_ABSCO_LUT="$VSMARTMOM_ABSCO_DIR/o2_v52_v2.jld2"
    export WEAK_H2O_ABSCO_LUT="$VSMARTMOM_ABSCO_DIR/wh2o_v52.jld2"
    export WEAK_CO2_ABSCO_LUT="$VSMARTMOM_ABSCO_DIR/wco2_v52.jld2"
    export STRONG_H2O_ABSCO_LUT="$VSMARTMOM_ABSCO_DIR/sh2o_v52.jld2"
    export STRONG_CO2_ABSCO_LUT="$VSMARTMOM_ABSCO_DIR/sco2_v52.jld2"
    R7A_SPECTROSCOPY_NAMES=(O2_ABSCO_LUT O2_H2O_LUT WEAK_H2O_ABSCO_LUT WEAK_CO2_ABSCO_LUT STRONG_H2O_ABSCO_LUT STRONG_CO2_ABSCO_LUT)
    local name
    for name in "${R7A_SPECTROSCOPY_NAMES[@]}"; do
        [[ -f "${!name}" && -r "${!name}" ]] || _r7a_fail "missing spectroscopy table: ${!name}" || return 1
    done
}

_r7a_verify_round7_inputs() {
    [[ -d "$ROUND7_BUNDLE" ]] || _r7a_fail "missing Round 7 release $ROUND7_BUNDLE" || return 1
    [[ "$(_r7a_sha "$ROUND7_BUNDLE/SHA256SUMS")" == "$ROUND7_RELEASE_SHA256SUMS" ]] ||
        _r7a_fail "Round 7 release manifest differs" || return 1
    (cd "$ROUND7_BUNDLE" && sha256sum --check --strict --quiet SHA256SUMS) ||
        _r7a_fail "Round 7 release files differ from their manifest" || return 1
    [[ "$(_r7a_sha "$ROUND7_REFERENCE_ROOT/observations/SHA256SUMS")" == "$ROUND7_OBSERVATION_SHA256SUMS" ]] ||
        _r7a_fail "Round 7 observation manifest differs" || return 1
    cmp -s "$ROUND7_BUNDLE/inputs/OBSERVATION_SHA256SUMS" "$ROUND7_REFERENCE_ROOT/observations/SHA256SUMS" ||
        _r7a_fail "Round 7 sealed observation manifest copy differs" || return 1
}

_r7a_verify_checkout() {
    local commit="$1" root="$ROUND7_ALTGRD_SOURCE_ROOT"
    [[ "$commit" =~ ^[0-9a-f]{40}$ ]] || _r7a_fail "checkpoint must be a full commit SHA" || return 1
    [[ "$(git -C "$root" rev-parse HEAD)" == "$commit" ]] ||
        _r7a_fail "checkout HEAD is not the pinned checkpoint $commit" || return 1
    [[ -z "$(git -C "$root" status --porcelain --untracked-files=all)" ]] ||
        _r7a_fail "pinned checkout is dirty" || return 1
    [[ ! -e "$root/LocalPreferences.toml" ]] || _r7a_fail "unsealed LocalPreferences.toml" || return 1
    cmp -s "$root/Manifest.toml" "$ROUND7_BUNDLE/inputs/Manifest.toml" ||
        _r7a_fail "Julia Manifest differs from the Round 7 release" || return 1
    [[ "$(_r7a_sha "$root/sandbox/workflows/RRS_XCO2/config/oco_grass_3aerosol.yaml")" == "$ROUND7_CONFIG_SHA256" ]] ||
        _r7a_fail "RT configuration differs from Round 7" || return 1
    [[ -x "$JULIA_BIN" && "$("$JULIA_BIN" --startup-file=no --version)" == 'julia version 1.12.5' ]] ||
        _r7a_fail "Round 7_altgrd requires Julia 1.12.5 at $JULIA_BIN" || return 1
}

# Pin the checkout once. The release record is never overwritten.
seal_round7_altgrd_release() {
    local commit="${1:-}" tree
    [[ -n "$commit" ]] || _r7a_fail "usage: seal_round7_altgrd_release <commit>" || return 1
    _r7a_paths || return 1
    _r7a_verify_round7_inputs || return 1
    _r7a_verify_checkout "$commit" || return 1
    tree="$(git -C "$ROUND7_ALTGRD_SOURCE_ROOT" ls-tree -r --full-tree HEAD | sha256sum | awk '{print $1}')"
    _r7a_spectroscopy || return 1
    local name spectroscopy=""
    for name in "${R7A_SPECTROSCOPY_NAMES[@]}"; do
        spectroscopy+="$(printf '    "%s": "%s",' "$name" "$(_r7a_sha "${!name}")")"$'\n'
    done
    spectroscopy="${spectroscopy%,$'\n'}"
    mkdir -p "$(dirname "$ROUND7_ALTGRD_RELEASE")"
    [[ ! -e "$ROUND7_ALTGRD_RELEASE" ]] || _r7a_fail "release already sealed: $ROUND7_ALTGRD_RELEASE" || return 1
    (set -o noclobber; printf '{\n  "release": "round7_altgrd_fixed_sif_imperfect_correction_pilot_v1",\n  "checkpoint": "%s",\n  "base_tree_sha256": "%s",\n  "branch": "experiment/round7-altgrd",\n  "reference_release_sha256sums": "%s",\n  "observation_manifest_sha256": "%s",\n  "spectroscopy_sha256": {\n%s\n  },\n  "sealed_utc": "%s"\n}\n' \
        "$commit" "$tree" "$ROUND7_RELEASE_SHA256SUMS" "$ROUND7_OBSERVATION_SHA256SUMS" "$spectroscopy" "$(date -u +%FT%TZ)" \
        > "$ROUND7_ALTGRD_RELEASE") || return 1
    python3 -c 'import json,sys; json.load(open(sys.argv[1]))' "$ROUND7_ALTGRD_RELEASE" ||
        _r7a_fail "sealed release is not valid JSON" || return 1
    chmod 400 "$ROUND7_ALTGRD_RELEASE"
    echo "Sealed $ROUND7_ALTGRD_RELEASE at $commit" >&2
}

configure_round7_altgrd() {
    local mode="${1:-}" checkpoint tree sealed_tree assignments key value count=0
    local -A identity_keys=()
    [[ "$mode" == off || "$mode" == on ]] || _r7a_fail "mode must be off or on" || return 1
    _r7a_paths || return 1
    [[ -f "$ROUND7_ALTGRD_RELEASE" ]] || _r7a_fail "unsealed release; run seal_round7_altgrd_release" || return 1
    checkpoint="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["checkpoint"])' "$ROUND7_ALTGRD_RELEASE")" || return 1
    sealed_tree="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["base_tree_sha256"])' "$ROUND7_ALTGRD_RELEASE")" || return 1
    _r7a_verify_round7_inputs || return 1
    _r7a_verify_checkout "$checkpoint" || return 1
    tree="$(git -C "$ROUND7_ALTGRD_SOURCE_ROOT" ls-tree -r --full-tree HEAD | sha256sum | awk '{print $1}')"
    [[ "$tree" == "$sealed_tree" ]] || _r7a_fail "checkout tree differs from the sealed release" || return 1

    # Identical spectroscopy (contents sealed at release time), solar, and inputs.
    _r7a_spectroscopy || return 1
    local name sealed
    for name in "${R7A_SPECTROSCOPY_NAMES[@]}"; do
        sealed="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["spectroscopy_sha256"][sys.argv[2]])' "$ROUND7_ALTGRD_RELEASE" "$name")" || return 1
        [[ "$(_r7a_sha "${!name}")" == "$sealed" ]] || _r7a_fail "spectroscopy table $name differs from the sealed release" || return 1
    done
    export SOLAR_OUT="${SOLAR_OUT:-$HOME/Raman_misc/workflows/worktrees/uni_vSmartMOM_sanghavi_2025-03-18/src/SolarModel/solar.out}"
    [[ "$(_r7a_sha "$SOLAR_OUT")" == "$ROUND7_SOLAR_SHA256" ]] || _r7a_fail "solar spectrum differs" || return 1
    export RRS_XCO2_CONFIG="$ROUND7_ALTGRD_SOURCE_ROOT/sandbox/workflows/RRS_XCO2/config/oco_grass_3aerosol.yaml"
    export RETRIEVAL_TRUTH_TABLE="$ROUND7_BUNDLE/inputs/true_states_corrected_sif_v2.dat"
    export RETRIEVAL_PRIOR_PATH="$ROUND7_BUNDLE/inputs/apriori_states_round5_fixed_sif_${mode}_tight_utls_acos_mapped_tapered_vertical_correlation.nc"
    export ROUND7_SOURCE_PRIOR="$ROUND7_BUNDLE/inputs/source_prior.nc"
    export RETRIEVAL_STOKES_COEFFICIENT_PATH="$ROUND7_BUNDLE/inputs/representative_stokes_coefficients.nc"
    export RETRIEVAL_SCENE_COMPONENTS_PATH="$ROUND7_BUNDLE/inputs/scene_components.dat"
    export RRS_XCO2_SIF_TEMPLATE_PATH="$ROUND7_BUNDLE/inputs/sif-spectra.csv"
    export ROUND7_OBSERVATION_DIR="$ROUND7_REFERENCE_ROOT/observations"
    unset RETRIEVAL_MEASUREMENT_DIR RETRIEVAL_NOISE_DIR FULL_COLUMN_TRUTH_ROOT
    unset ROUND7_ALTGRD_PRINT_IDENTITY ROUND7_ALTGRD_PREFLIGHT_ONLY
    export SIF_CASE_FILTER="$mode" ROUND7_KNOWN_WAVELENGTH_NM=759
    export RETRIEVAL_OUTPUT_ROOT="$ROUND7_ALTGRD_ROOT/retrievals_nosif"
    [[ "$mode" != on ]] || export RETRIEVAL_OUTPUT_ROOT="$ROUND7_ALTGRD_ROOT/retrievals_sif"
    export RETRIEVAL_CLASS=corrected RETRIEVAL_ARCH=GPU RETRIEVAL_FLOAT_TYPE=Float32 RETRIEVAL_NSTREAMS=9
    export RETRIEVAL_WRITE_MANIFEST=0 FORCE=0 FAIL_FAST=1
    export JULIA_CPU_TARGET=generic JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
    export JULIA_LOAD_PATH='@:@stdlib' JULIA_PKG_PRECOMPILE_AUTO=0
    export VSMARTMOM_CUDA_BATCH_INV_PREALLOC=0
    # Round 7's depot is read, never written: compiled caches for this branch
    # go to a private first depot inside the new output tree.
    mkdir -p "$ROUND7_ALTGRD_ROOT/.julia_depot_altgrd_curry"
    export JULIA_DEPOT_PATH="$ROUND7_ALTGRD_ROOT/.julia_depot_altgrd_curry:$RRS_XCO2_DATA_ROOT/.julia_depot_curry:$HOME/.julia:"
    export ROUND7_ALTGRD_CODE_CHECKPOINT_SHA="$checkpoint"

    # The runner computes the identity from the pinned files; never hand-typed.
    assignments="$(ROUND7_ALTGRD_PRINT_IDENTITY=1 CUDA_VISIBLE_DEVICES='' "$JULIA_BIN" \
        --project="$ROUND7_ALTGRD_SOURCE_ROOT" --startup-file=no \
        "$ROUND7_ALTGRD_SOURCE_ROOT/sandbox/workflows/RRS_XCO2/inversion/run_round7_altgrd_retrievals.jl" |
        grep '^export ')" || _r7a_fail "identity computation failed" || return 1
    while IFS='=' read -r key value; do
        key="${key#export }"
        case "$key" in
            ROUND7_ALTGRD_CODE_CHECKPOINT_SHA)
                [[ "$value" == "$checkpoint" ]] || _r7a_fail "identity checkpoint differs" || return 1 ;;
            ROUND7_ALTGRD_CODESET_SHA256|ROUND7_ALTGRD_INPUT_SET_SHA256|ROUND7_ALTGRD_CAMPAIGN_IDENTITY_SHA256|ROUND7_ALTGRD_OBSERVATION_MANIFEST_SHA256)
                [[ "$value" =~ ^[0-9a-f]{64}$ ]] || _r7a_fail "malformed identity value for $key" || return 1 ;;
            *) _r7a_fail "unexpected identity key: $key" || return 1 ;;
        esac
        [[ -z "${identity_keys[$key]:-}" ]] || _r7a_fail "duplicate identity key: $key" || return 1
        identity_keys[$key]=1
        printf -v "$key" '%s' "$value"
        export "$key"
        count=$((count + 1))
    done <<< "$assignments"
    [[ "$count" == 5 ]] || _r7a_fail "incomplete identity" || return 1
    echo "Round7_altgrd SIF=$mode checkpoint=$checkpoint identity=$ROUND7_ALTGRD_CAMPAIGN_IDENTITY_SHA256" >&2
}
