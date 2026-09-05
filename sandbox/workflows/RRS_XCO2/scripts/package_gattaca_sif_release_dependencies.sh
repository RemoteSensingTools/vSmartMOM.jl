#!/usr/bin/env bash

# Build or install the immutable, non-spectroscopy inputs needed to publish a
# corrected full-column SIF release and rebuild the bottom-layer SIF campaign
# on a disconnected compute host.  The payload is data-only: source code must
# arrive through the separately pinned Git checkout.
#
# Package (on the authoritative-data host):
#   EXPECTED_SIF_OWNERSHIP_SHA256=<sha256> \
#     ./package_gattaca_sif_release_dependencies.sh package /safe/output/dir
#
# Install (on the compute host):
#   RRS_REPO=/path/to/pinned/checkout \
#   RRS_PRIVATE_ROOT=/private/root \
#     ./package_gattaca_sif_release_dependencies.sh install /transfer/dir
#
# Installation is fail-closed and idempotent.  It first validates every byte
# in an isolated private staging directory, then refuses any existing target
# whose checksum differs.  Missing files are linked into place atomically; it
# never overwrites a destination.

set -euo pipefail
umask 077

readonly archive_name="sif_wavelength_integral_0p5_20260904"
readonly repo_tar_name="gattaca_sif_release_repo_dependencies.tar"
readonly archive_tar_name="gattaca_sif_release_full_legacy_archive.tar"
readonly repo_manifest_name="gattaca_sif_release_repo_dependencies.sha256"
readonly archive_manifest_name="gattaca_sif_release_full_legacy_archive.sha256"
readonly metadata_name="BUNDLE_METADATA.dat"
readonly envelope_name="TRANSFER_SHA256SUMS"
readonly packager_relative_path="RRS_XCO2/scripts/package_gattaca_sif_release_dependencies.sh"

die() {
    printf 'STOP: %s\n' "$*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  package_gattaca_sif_release_dependencies.sh package OUTPUT_DIRECTORY
  package_gattaca_sif_release_dependencies.sh install BUNDLE_DIRECTORY

Package mode requires EXPECTED_SIF_OWNERSHIP_SHA256 and requires this helper
itself to match the version committed at HEAD. Unrelated tracked edits outside
the explicitly enumerated data payload do not block packaging; every payload
byte is independently hashed and the completed archives are re-extracted and
verified. Install mode requires RRS_PRIVATE_ROOT and a clean tracked checkout.
RRS_REPO defaults to the Git top level containing this script.
EOF
}

script_repo="$({ git -C "$(dirname "$0")" rev-parse --show-toplevel; } 2>/dev/null)" ||
    die "this script must be run from a vSmartMOM Git checkout"
repo_root="$(realpath "${RRS_REPO:-$script_repo}")"
[[ -d "$repo_root/.git" ]] || die "RRS_REPO is not a Git checkout: $repo_root"
readonly packager_path="$(realpath "$0")"
readonly expected_packager_path="$(realpath -m "$repo_root/$packager_relative_path")"

require_file() {
    [[ -f "$1" ]] || die "missing required file: $1"
}

require_absent() {
    [[ ! -e "$1" ]] || die "pre-publication input unexpectedly exists: $1"
}

require_count() {
    local expected="$1"
    local directory="$2"
    local pattern="$3"
    local actual
    actual="$(find "$directory" -maxdepth 1 -type f -name "$pattern" -print | wc -l)"
    [[ "$actual" -eq "$expected" ]] ||
        die "$directory contains $actual $pattern files; expected $expected"
}

sha_of() {
    sha256sum "$1" | awk '{print $1}'
}

# Package metadata names the Git checkpoint whose source will be used on the
# disconnected host.  A dirty, unrelated tracked file in the authoritative
# data worktree must not prevent a data-only handoff, but the helper defining
# that handoff must be byte-for-byte identical to its version at that exact
# checkpoint.  Comparing Git blob IDs is independent of the index state and
# therefore catches both staged and unstaged edits to this script.
require_packager_at_head() {
    local head_blob worktree_blob object_type
    [[ "$packager_path" == "$expected_packager_path" ]] ||
        die "run the canonical packager from $expected_packager_path"
    head_blob="$(git -C "$repo_root" rev-parse --verify \
        "HEAD:$packager_relative_path")" ||
        die "packager is absent from the current Git checkpoint"
    object_type="$(git -C "$repo_root" cat-file -t "$head_blob")" ||
        die "cannot inspect the committed packager object"
    [[ "$object_type" == blob ]] || die "committed packager object is not a file"
    worktree_blob="$(git -C "$repo_root" hash-object -- "$packager_path")" ||
        die "cannot hash the working packager"
    [[ "$worktree_blob" == "$head_blob" ]] ||
        die "packager differs from $packager_relative_path at HEAD"
    sha_of "$packager_path"
}

validate_legacy_full_column_archive() {
    local full_truth="$repo_root/RRS_XCO2/truth_map"
    local legacy_root="$repo_root/RRS_XCO2/obsolete/$archive_name"
    local manifest="$legacy_root/full_column_input_manifest.dat"
    local expected_table_hash
    local row_count=0
    local off_index on_index expected_hash name archive_file canonical_file

    require_file "$manifest"
    require_file "$legacy_root/truth_map/true_states.dat"
    expected_table_hash="$(awk '$1 == "#" && $2 == "legacy_true_states_sha256" {print $3}' "$manifest")"
    [[ "$expected_table_hash" =~ ^[0-9a-f]{64}$ ]] ||
        die "legacy archive manifest lacks a valid state-table checksum"
    [[ "$(sha_of "$legacy_root/truth_map/true_states.dat")" == "$expected_table_hash" ]] ||
        die "legacy full-column state table differs from its archive manifest"
    [[ "$(sha_of "$full_truth/true_states.dat")" == "$expected_table_hash" ]] ||
        die "canonical full-column table is not the archived pre-publication table"

    while read -r off_index on_index expected_hash _; do
        [[ -n "${off_index:-}" ]] || continue
        [[ "$off_index" == "#" ]] && continue
        [[ "$off_index" =~ ^[0-9]{3}$ && "$on_index" =~ ^[0-9]{3}$ &&
           "$expected_hash" =~ ^[0-9a-f]{64}$ ]] ||
            die "malformed full-column archive manifest row: $off_index $on_index $expected_hash"
        name="hiressim_${on_index}.nc"
        if [[ -f "$legacy_root/truth_map/$name" &&
              -f "$legacy_root/truth_map/aerosol_chunked/$name" ]]; then
            die "legacy SIF state $on_index occurs in both scene classes"
        elif [[ -f "$legacy_root/truth_map/$name" ]]; then
            archive_file="$legacy_root/truth_map/$name"
            canonical_file="$full_truth/$name"
        elif [[ -f "$legacy_root/truth_map/aerosol_chunked/$name" ]]; then
            archive_file="$legacy_root/truth_map/aerosol_chunked/$name"
            canonical_file="$full_truth/aerosol_chunked/$name"
        else
            die "legacy archive lacks full-column SIF state $on_index"
        fi
        require_file "$canonical_file"
        [[ "$(sha_of "$archive_file")" == "$expected_hash" ]] ||
            die "archived full-column state $on_index differs from its manifest"
        [[ "$(sha_of "$canonical_file")" == "$expected_hash" ]] ||
            die "canonical full-column state $on_index is not the archived legacy input"
        row_count=$((row_count + 1))
    done < "$manifest"
    [[ "$row_count" -eq 32 ]] ||
        die "full-column archive manifest contains $row_count scene rows; expected 32"
}

# Build the legacy payload from the validated manifest, not from a recursive
# directory sweep.  This makes an unexpected note, credential, or intermediate
# file a hard failure instead of silently adding it to a transfer bundle.
collect_archive_files() {
    local obsolete_root="$repo_root/RRS_XCO2/obsolete"
    local legacy_root="$obsolete_root/$archive_name"
    local manifest="$legacy_root/full_column_input_manifest.dat"
    local off_index on_index expected_hash name relative
    local -a actual_files

    ARCHIVE_FILES=(
        "$archive_name/README.md"
        "$archive_name/full_column_input_manifest.dat"
        "$archive_name/truth_map/true_states.dat"
    )
    while read -r off_index on_index expected_hash _; do
        [[ -n "${off_index:-}" ]] || continue
        [[ "$off_index" == "#" ]] && continue
        name="hiressim_${on_index}.nc"
        if [[ -f "$legacy_root/truth_map/$name" ]]; then
            relative="$archive_name/truth_map/$name"
        elif [[ -f "$legacy_root/truth_map/aerosol_chunked/$name" ]]; then
            relative="$archive_name/truth_map/aerosol_chunked/$name"
        else
            die "validated legacy archive lost SIF state $on_index"
        fi
        ARCHIVE_FILES+=("$relative")
    done < "$manifest"
    mapfile -t ARCHIVE_FILES < <(
        printf '%s\n' "${ARCHIVE_FILES[@]}" | LC_ALL=C sort -u)
    [[ "${#ARCHIVE_FILES[@]}" -eq 35 ]] ||
        die "legacy archive allowlist has ${#ARCHIVE_FILES[@]} files; expected 35"

    mapfile -t actual_files < <(
        cd "$obsolete_root"
        find "$archive_name" -type f -print | LC_ALL=C sort
    )
    [[ "${#actual_files[@]}" -eq "${#ARCHIVE_FILES[@]}" ]] ||
        die "legacy archive contains files outside its 35-member allowlist"
    for index in "${!ARCHIVE_FILES[@]}"; do
        [[ "${actual_files[$index]}" == "${ARCHIVE_FILES[$index]}" ]] ||
            die "unexpected legacy archive file: ${actual_files[$index]}"
    done
}

collect_repo_files() {
    local full_truth="$repo_root/RRS_XCO2/truth_map"
    local bottom_root="$repo_root/RRS_XCO2/bottom_layer_XCO2_retrievals"
    local bottom_truth="$bottom_root/truth"
    local measurement_root="$bottom_truth/OCO_radiances"
    local noise_root="$measurement_root/noise_covariances"
    local path

    REPO_FILES=(
        "RRS_XCO2/truth_map/true_states.dat"
        "RRS_XCO2/truth_map/scene_components.dat"
        "RRS_XCO2/truth_map/sim_wavelength.nc"
        "RRS_XCO2/truth_map/aerosol_chunked/sim_wavelength.nc"
        "RRS_XCO2/bottom_layer_XCO2_retrievals/truth/true_states.dat"
        "RRS_XCO2/bottom_layer_XCO2_retrievals/truth/control_reuse_map.dat"
        "RRS_XCO2/bottom_layer_XCO2_retrievals/truth/scene_components.dat"
        "RRS_XCO2/bottom_layer_XCO2_retrievals/truth/sim_wavelength.nc"
        "RRS_XCO2/bottom_layer_XCO2_retrievals/truth/aerosol_chunked/sim_wavelength.nc"
        "RRS_XCO2/bottom_layer_XCO2_retrievals/truth/OCO_radiances/noise_covariances/noise_covariance_manifest.dat"
        "RRS_XCO2/bottom_layer_XCO2_retrievals/retrievals/.control/sif_owned_externally"
        "RRS_XCO2/inversion/instrument/representative_stokes_coefficients.nc"
        "RRS_XCO2/inversion/instrument/representative_snr_coefficients.nc"
    )

    require_count 32 "$full_truth" 'hiressim_*.nc'
    require_count 32 "$full_truth/aerosol_chunked" 'hiressim_*.nc'
    require_count 40 "$bottom_truth" 'hiressim_*.nc'
    require_count 40 "$bottom_truth/aerosol_chunked" 'hiressim_*.nc'
    require_count 80 "$measurement_root" 'OCO2sims_*.nc'
    require_count 80 "$noise_root" 'OCO2noise_*.nc'

    while IFS= read -r -d '' path; do
        REPO_FILES+=("${path#"$repo_root/"}")
    done < <(find "$full_truth" "$full_truth/aerosol_chunked" \
        -maxdepth 1 -type f -name 'hiressim_*.nc' -print0 | sort -z)
    while IFS= read -r -d '' path; do
        REPO_FILES+=("${path#"$repo_root/"}")
    done < <(find "$bottom_truth" "$bottom_truth/aerosol_chunked" \
        -maxdepth 1 -type f -name 'hiressim_*.nc' -print0 | sort -z)
    while IFS= read -r -d '' path; do
        REPO_FILES+=("${path#"$repo_root/"}")
    done < <(find "$measurement_root" -maxdepth 1 -type f \
        -name 'OCO2sims_*.nc' -print0 | sort -z)
    while IFS= read -r -d '' path; do
        REPO_FILES+=("${path#"$repo_root/"}")
    done < <(find "$noise_root" -maxdepth 1 -type f \
        -name 'OCO2noise_*.nc' -print0 | sort -z)

    mapfile -t REPO_FILES < <(printf '%s\n' "${REPO_FILES[@]}" | LC_ALL=C sort -u)
    for path in "${REPO_FILES[@]}"; do
        [[ "$path" != /* && "$path" != ../* && "$path" != */../* ]] ||
            die "unsafe repository payload path: $path"
        require_file "$repo_root/$path"
    done
}

package_bundle() {
    local output_directory="$1"
    local ownership="$repo_root/RRS_XCO2/bottom_layer_XCO2_retrievals/retrievals/.control/sif_owned_externally"
    local expected_ownership="${EXPECTED_SIF_OWNERSHIP_SHA256:-}"
    local full_truth="$repo_root/RRS_XCO2/truth_map"
    local bottom_truth="$repo_root/RRS_XCO2/bottom_layer_XCO2_retrievals/truth"
    local obsolete_root="$repo_root/RRS_XCO2/obsolete"
    local legacy_root="$obsolete_root/$archive_name"
    local source_sha packager_sha256 created_utc path

    [[ "$expected_ownership" =~ ^[0-9a-f]{64}$ ]] ||
        die "EXPECTED_SIF_OWNERSHIP_SHA256 must be a lowercase SHA-256"
    packager_sha256="$(require_packager_at_head)"
    source_sha="$(git -C "$repo_root" rev-parse HEAD)"
    [[ "$source_sha" =~ ^[0-9a-f]{40}$ ]] || die "cannot resolve source checkpoint"
    case "$output_directory/" in
        "$repo_root/"*) die "bundle output must be outside the Git checkout" ;;
    esac

    require_absent "$full_truth/sif_v2_release_complete.dat"
    require_absent "$full_truth/.sif_v2_publication_in_progress"
    require_absent "$bottom_truth/bottom_layer_sif_v2_release_complete.dat"
    require_absent "$bottom_truth/.bottom_layer_sif_v2_publication_in_progress"
    require_file "$ownership"
    [[ "$(sha_of "$ownership")" == "$expected_ownership" ]] ||
        die "external-SIF ownership marker differs from its approved checksum"
    awk '!/^#/ {if ($7 == "total_0p5") n++} END {exit n == 32 ? 0 : 1}' \
        "$full_truth/true_states.dat" ||
        die "full-column table does not contain exactly 32 legacy SIF-on rows"
    awk '!/^#/ {if ($7 == "total_0p5") n++} END {exit n == 40 ? 0 : 1}' \
        "$bottom_truth/true_states.dat" ||
        die "bottom-layer table does not contain exactly 40 legacy SIF-on rows"
    validate_legacy_full_column_archive
    collect_archive_files
    collect_repo_files

    [[ ! -e "$output_directory" ]] ||
        die "output directory already exists: $output_directory"
    mkdir -p "$(dirname "$output_directory")"
    mkdir "$output_directory"

    (
        cd "$repo_root"
        sha256sum "${REPO_FILES[@]}"
    ) > "$output_directory/$repo_manifest_name"
    (
        cd "$obsolete_root"
        sha256sum "${ARCHIVE_FILES[@]}"
    ) > "$output_directory/$archive_manifest_name"
    (
        cd "$repo_root"
        printf '%s\0' "${REPO_FILES[@]}" | \
            tar --null --no-recursion --files-from=- -cf \
                "$output_directory/$repo_tar_name"
    )
    (
        cd "$obsolete_root"
        printf '%s\0' "${ARCHIVE_FILES[@]}" | \
            tar --null --no-recursion --files-from=- -cf \
                "$output_directory/$archive_tar_name"
    )

    # Detect any source mutation between manifest creation and tar capture and
    # prove that each completed tar contains exactly the bytes authorized by
    # its manifest before describing the bundle as verified.
    verify_packaged_archives "$output_directory"

    created_utc="$(date -u +%FT%TZ)"
    {
        printf 'bundle_schema=2\n'
        printf 'source_checkpoint=%s\n' "$source_sha"
        printf 'packager_path=%s\n' "$packager_relative_path"
        printf 'packager_sha256=%s\n' "$packager_sha256"
        printf 'created_utc=%s\n' "$created_utc"
        printf 'repository_file_count=%s\n' "${#REPO_FILES[@]}"
        printf 'archive_file_count=%s\n' "${#ARCHIVE_FILES[@]}"
        printf 'ownership_marker_sha256=%s\n' "$expected_ownership"
        printf 'repository_manifest_sha256=%s\n' \
            "$(sha_of "$output_directory/$repo_manifest_name")"
        printf 'archive_manifest_sha256=%s\n' \
            "$(sha_of "$output_directory/$archive_manifest_name")"
    } > "$output_directory/$metadata_name"
    (
        cd "$output_directory"
        sha256sum "$repo_tar_name" "$archive_tar_name" \
            "$repo_manifest_name" "$archive_manifest_name" "$metadata_name" \
            > "$envelope_name"
        sha256sum --check --strict "$envelope_name"
    )
    printf 'Static release bundle created and verified:\n  %s\n' "$output_directory"
}

safe_tar_members() {
    local archive="$1"
    local member
    while IFS= read -r member; do
        [[ -n "$member" && "$member" != /* && "$member" != ../* &&
           "$member" != */../* && "$member" != */.. ]] ||
            die "unsafe tar member in $archive: $member"
    done < <(tar -tf "$archive")
}

manifest_path() {
    local raw="$1"
    raw="${raw#\*}"
    printf '%s\n' "$raw"
}

preflight_tree() {
    local source_root="$1"
    local destination_root="$2"
    local manifest="$3"
    local destination_root_real digest raw_relative relative source destination
    destination_root_real="$(realpath -m "$destination_root")"
    while read -r digest raw_relative; do
        [[ "$digest" =~ ^[0-9a-f]{64}$ ]] || die "invalid manifest digest: $digest"
        relative="$(manifest_path "$raw_relative")"
        [[ -n "$relative" && "$relative" != /* && "$relative" != ../* &&
           "$relative" != */../* && "$relative" != */.. ]] ||
            die "unsafe manifest path: $relative"
        source="$source_root/$relative"
        require_file "$source"
        [[ "$(sha_of "$source")" == "$digest" ]] ||
            die "staged payload checksum mismatch: $source"
        destination="$destination_root/$relative"
        case "$(realpath -m "$destination")" in
            "$destination_root_real"/*) ;;
            *) die "payload target escapes its approved root: $destination" ;;
        esac
        if [[ -e "$destination" ]]; then
            [[ -f "$destination" ]] || die "existing target is not a file: $destination"
            [[ "$(sha_of "$destination")" == "$digest" ]] ||
                die "existing target differs; nothing was installed: $destination"
        fi
    done < "$manifest"
}

install_tree() {
    local source_root="$1"
    local destination_root="$2"
    local manifest="$3"
    local digest raw_relative relative source destination temporary
    while read -r digest raw_relative; do
        relative="$(manifest_path "$raw_relative")"
        source="$source_root/$relative"
        destination="$destination_root/$relative"
        if [[ -f "$destination" ]]; then
            continue
        fi
        mkdir -p "$(dirname "$destination")"
        temporary="$(mktemp "$(dirname "$destination")/.incoming.$(basename "$destination").XXXXXX")"
        cp -- "$source" "$temporary"
        [[ "$(sha_of "$temporary")" == "$digest" ]] ||
            die "temporary copy checksum mismatch: $temporary"
        if ln -- "$temporary" "$destination" 2>/dev/null; then
            rm -- "$temporary"
        else
            [[ -f "$destination" && "$(sha_of "$destination")" == "$digest" ]] ||
                die "target appeared during installation with different bytes: $destination"
            rm -- "$temporary"
        fi
    done < "$manifest"
}

verify_tree() {
    local destination_root="$1"
    local manifest="$2"
    (
        cd "$destination_root"
        sha256sum --check --strict --quiet "$manifest"
    )
}

verify_packaged_archives() {
    local bundle_directory="$1"
    local validation_root repo_stage archive_stage
    safe_tar_members "$bundle_directory/$repo_tar_name"
    safe_tar_members "$bundle_directory/$archive_tar_name"
    validation_root="$(mktemp -d \
        "$(dirname "$bundle_directory")/.sif_bundle_verify.XXXXXX")"
    repo_stage="$validation_root/repo"
    archive_stage="$validation_root/archive"
    mkdir "$repo_stage" "$archive_stage"
    tar --no-same-owner --no-same-permissions \
        -xf "$bundle_directory/$repo_tar_name" -C "$repo_stage"
    tar --no-same-owner --no-same-permissions \
        -xf "$bundle_directory/$archive_tar_name" -C "$archive_stage"
    verify_tree "$repo_stage" "$bundle_directory/$repo_manifest_name"
    verify_tree "$archive_stage" "$bundle_directory/$archive_manifest_name"
    rm -rf -- "$validation_root"
}

install_bundle() {
    local bundle_directory="$(realpath "$1")"
    local private_root="${RRS_PRIVATE_ROOT:-}"
    local bundle_schema source_sha actual_sha bundle_digest stage_root
    local metadata_packager_path metadata_packager_sha256 local_packager_sha256
    local repo_stage archive_stage

    [[ -n "$private_root" ]] || die "RRS_PRIVATE_ROOT is required in install mode"
    mkdir -p "$private_root" "$private_root/transfer_staging"
    private_root="$(realpath "$private_root")"
    [[ "$private_root" != "$repo_root" ]] || die "private root cannot equal the Git checkout"
    case "$private_root/" in
        "$repo_root/"*) die "private root cannot be inside the Git checkout" ;;
    esac

    for path in "$repo_tar_name" "$archive_tar_name" "$repo_manifest_name" \
                "$archive_manifest_name" "$metadata_name" "$envelope_name"; do
        require_file "$bundle_directory/$path"
    done
    (
        cd "$bundle_directory"
        sha256sum --check --strict "$envelope_name"
    )
    bundle_schema="$(awk -F= '$1 == "bundle_schema" {print $2}' \
        "$bundle_directory/$metadata_name")"
    [[ "$bundle_schema" == 2 ]] ||
        die "unsupported static-release bundle schema: $bundle_schema"
    metadata_packager_path="$(awk -F= '$1 == "packager_path" {print $2}' \
        "$bundle_directory/$metadata_name")"
    [[ "$metadata_packager_path" == "$packager_relative_path" ]] ||
        die "bundle names an unexpected packager path: $metadata_packager_path"
    metadata_packager_sha256="$(awk -F= \
        '$1 == "packager_sha256" {print $2}' \
        "$bundle_directory/$metadata_name")"
    [[ "$metadata_packager_sha256" =~ ^[0-9a-f]{64}$ ]] ||
        die "bundle lacks a valid packager checksum"
    local_packager_sha256="$(require_packager_at_head)"
    [[ "$local_packager_sha256" == "$metadata_packager_sha256" ]] ||
        die "bundle and checkout use different packager bytes"
    source_sha="$(awk -F= '$1 == "source_checkpoint" {print $2}' \
        "$bundle_directory/$metadata_name")"
    actual_sha="$(git -C "$repo_root" rev-parse HEAD)"
    [[ "$source_sha" =~ ^[0-9a-f]{40}$ && "$actual_sha" == "$source_sha" ]] ||
        die "bundle source checkpoint $source_sha differs from checkout $actual_sha"
    git -C "$repo_root" diff --quiet -- || die "tracked checkout has unstaged changes"
    git -C "$repo_root" diff --cached --quiet -- || die "tracked checkout has staged changes"

    safe_tar_members "$bundle_directory/$repo_tar_name"
    safe_tar_members "$bundle_directory/$archive_tar_name"
    bundle_digest="$(sha_of "$bundle_directory/$envelope_name")"
    stage_root="$(mktemp -d "$private_root/transfer_staging/sif_release_${bundle_digest:0:16}.XXXXXX")"
    repo_stage="$stage_root/repo"
    archive_stage="$stage_root/archive"
    mkdir "$repo_stage" "$archive_stage"
    tar --no-same-owner --no-same-permissions \
        -xf "$bundle_directory/$repo_tar_name" -C "$repo_stage"
    tar --no-same-owner --no-same-permissions \
        -xf "$bundle_directory/$archive_tar_name" -C "$archive_stage"
    verify_tree "$repo_stage" "$bundle_directory/$repo_manifest_name"
    verify_tree "$archive_stage" "$bundle_directory/$archive_manifest_name"

    # Preflight both roots before the first installation write. Existing
    # identical files are accepted; a single differing byte blocks everything.
    preflight_tree "$repo_stage" "$repo_root" \
        "$bundle_directory/$repo_manifest_name"
    preflight_tree "$archive_stage" "$private_root/archive" \
        "$bundle_directory/$archive_manifest_name"
    install_tree "$repo_stage" "$repo_root" \
        "$bundle_directory/$repo_manifest_name"
    install_tree "$archive_stage" "$private_root/archive" \
        "$bundle_directory/$archive_manifest_name"
    verify_tree "$repo_root" "$bundle_directory/$repo_manifest_name"
    verify_tree "$private_root/archive" "$bundle_directory/$archive_manifest_name"
    printf 'Static release dependencies installed and verified.\n'
    printf 'Audit staging retained at:\n  %s\n' "$stage_root"
}

[[ $# -eq 2 ]] || {
    usage >&2
    exit 2
}
mode="$1"
target="$2"
case "$mode" in
    package) package_bundle "$(realpath -m "$target")" ;;
    install) install_bundle "$target" ;;
    *) usage >&2; exit 2 ;;
esac
