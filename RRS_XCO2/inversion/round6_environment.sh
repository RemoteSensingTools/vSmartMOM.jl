#!/usr/bin/env bash
# Source this file from a released Round-6 bundle, then configure_round6 off|on local|gattaca.

configure_round6() {
    local mode="$1" site="$2" base_tree release_hash input_hash codeset identity
    local bundle_dir data_checkout external_inputs prior_name
    bundle_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
    export ROUND6_BUNDLE="$bundle_dir"
    case "$site" in
        local)
            export ROUND6_SOURCE_ROOT="${ROUND6_SOURCE_ROOT:-$HOME/code/github/uni_vSmartMOM_round5_jacobian}"
            data_checkout="${ROUND6_DATA_CHECKOUT:-$HOME/code/github/uni_vSmartMOM}"
            export JULIA_BIN="${JULIA_BIN:-$HOME/.julia/juliaup/julia-1.12.5+0.x64.linux.gnu/bin/julia}"
            export VSMARTMOM_ABSCO_DIR="${VSMARTMOM_ABSCO_DIR:-/net/fluo/data1/ABSCO_CS_Database/v5.2_final}"
            export O2_H2O_LUT="${O2_H2O_LUT:-$HOME/data/HITRAN_LUTs/H2O.jld2}"
            export SOLAR_OUT="${SOLAR_OUT:-$HOME/Raman_misc/workflows/worktrees/uni_vSmartMOM_sanghavi_2025-03-18/src/SolarModel/solar.out}"
            ;;
        gattaca)
            export ROUND6_SOURCE_ROOT="${ROUND6_SOURCE_ROOT:-$HOME/code/uni_vSmartMOM_round5_jacobian}"
            data_checkout="${ROUND6_DATA_CHECKOUT:-$HOME/code/uni_vSmartMOM}"
            external_inputs="${ROUND6_PRIVATE_ROOT:-$HOME/RRS_XCO2_private}/rrs_inputs"
            export JULIA_BIN="${JULIA_BIN:-$HOME/software/julia-1.12.5/bin/julia}"
            export VSMARTMOM_ABSCO_DIR="$external_inputs/absco_v5.2"
            export O2_H2O_LUT="$external_inputs/hitran/H2O.jld2"
            export SOLAR_OUT="$external_inputs/solar/solar.out"
            export JULIA_DEPOT_PATH="${ROUND6_DEPOT:-$HOME/RRS_XCO2_private/julia_depot_1.12.5:}"
            ;;
        *) echo "unknown Round-6 site: $site" >&2; return 1 ;;
    esac
    [[ "$mode" == off || "$mode" == on ]] || return 1
    [[ -x "$JULIA_BIN" ]] || { echo "missing Julia: $JULIA_BIN" >&2; return 1; }
    (cd "$bundle_dir" && sha256sum --check --strict SHA256SUMS)
    [[ "$(git -C "$ROUND6_SOURCE_ROOT" rev-parse HEAD)" == 7acab57000dae259207a6760faae156cfa1734f6 ]]
    [[ -z "$(git -C "$ROUND6_SOURCE_ROOT" status --porcelain --untracked-files=all)" ]]
    base_tree="$(git -C "$ROUND6_SOURCE_ROOT" ls-tree -r --full-tree HEAD | sha256sum | awk '{print $1}')"
    [[ "$base_tree" == 1b2c4061bd9d7b3131523d2c76269cf103f0cc4f014904d69b4e9e161cc940f3 ]]
    cmp "$bundle_dir/inputs/Manifest.toml" "$ROUND6_SOURCE_ROOT/Manifest.toml"
    export RRS_XCO2_DATA_ROOT="$data_checkout/RRS_XCO2"
    export FULL_COLUMN_TRUTH_ROOT="$RRS_XCO2_DATA_ROOT/truth_map"
    export RRS_XCO2_CONFIG="$ROUND6_SOURCE_ROOT/sandbox/workflows/RRS_XCO2/config/oco_grass_3aerosol.yaml"
    export RETRIEVAL_TRUTH_TABLE="$bundle_dir/inputs/true_states_corrected_sif_v2.dat"
    if [[ "$site" == gattaca ]]; then
        # The SIF release barrier resolves high-resolution truth relative to this table.
        export RETRIEVAL_TRUTH_TABLE="$RRS_XCO2_DATA_ROOT/bottom_layer_XCO2_retrievals/truth/true_states.dat"
        cmp "$RETRIEVAL_TRUTH_TABLE" "$bundle_dir/inputs/true_states_corrected_sif_v2.dat"
    fi
    export RETRIEVAL_MEASUREMENT_DIR="$RRS_XCO2_DATA_ROOT/bottom_layer_XCO2_retrievals/truth/OCO_radiances"
    export RETRIEVAL_NOISE_DIR="$RETRIEVAL_MEASUREMENT_DIR/noise_covariances"
    export RETRIEVAL_STOKES_COEFFICIENT_PATH="$bundle_dir/inputs/representative_stokes_coefficients.nc"
    export RETRIEVAL_SCENE_COMPONENTS_PATH="$bundle_dir/inputs/scene_components.dat"
    export RRS_XCO2_SIF_TEMPLATE_PATH="$bundle_dir/inputs/sif-spectra.csv"
    export ROUND6_SOURCE_PRIOR="$bundle_dir/inputs/source_prior.nc"
    prior_name="apriori_states_round6_fixed_sif_${mode}_standard_utls_acos_mapped_tapered_vertical_correlation.nc"
    export RETRIEVAL_PRIOR_PATH="$bundle_dir/inputs/$prior_name"
    [[ -d "$RETRIEVAL_NOISE_DIR" && -d "$VSMARTMOM_ABSCO_DIR" && -f "$O2_H2O_LUT" ]]
    [[ "$(sha256sum "$SOLAR_OUT" | awk '{print $1}')" == 9e44e424788ce9f91c654398a789cd9da80205b9d6dc0c00df5e20f500ee8644 ]]
    release_hash="$(sha256sum "$bundle_dir/SHA256SUMS" | awk '{print $1}')"
    codeset="$(printf '%s\n%s\n' "$base_tree" "$release_hash" | sha256sum | awk '{print $1}')"
    input_hash="$(sha256sum "$RETRIEVAL_TRUTH_TABLE" "$RETRIEVAL_PRIOR_PATH" \
        "$ROUND6_SOURCE_PRIOR" "$RETRIEVAL_STOKES_COEFFICIENT_PATH" \
        "$RETRIEVAL_SCENE_COMPONENTS_PATH" "$RRS_XCO2_SIF_TEMPLATE_PATH" \
        "$SOLAR_OUT" "$RRS_XCO2_CONFIG" | awk '{print $1}' | sha256sum | awk '{print $1}')"
    identity="$(printf '%s\n%s\n%s\n' round6_fixed_sif_standard_utls_v1 "$codeset" "$input_hash" | sha256sum | awk '{print $1}')"
    export ROUND6_CODE_CHECKPOINT_SHA=7acab57000dae259207a6760faae156cfa1734f6
    export ROUND6_CODESET_SHA256="$codeset" ROUND6_INPUT_SET_SHA256="$input_hash"
    export ROUND6_CAMPAIGN_IDENTITY_SHA256="$identity"
    export ROUND6_SIF_MODE="$mode" SIF_CASE_FILTER="$mode" ROUND6_KNOWN_WAVELENGTH_NM=759
    export RETRIEVAL_CLASS=paired RETRIEVAL_ARCH=GPU RETRIEVAL_FLOAT_TYPE=Float32
    export RETRIEVAL_WRITE_MANIFEST=0 FORCE=0 FAIL_FAST=1
    export JULIA_CPU_TARGET=generic JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}" OPENBLAS_NUM_THREADS=1
    echo "Round 6: SIF=$mode; fixed SIF only; original UTLS and CO2 priors"
    echo "source=$ROUND6_SOURCE_ROOT release_sha256=$release_hash campaign_identity=$identity"
}
