#!/usr/bin/env bash
# Package new Round-6 workflow code and small inputs; never copy or modify truth spectra.
set -euo pipefail
umask 077
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
rrs_root="$(dirname "$script_dir")"
bottom="$rrs_root/bottom_layer_XCO2_retrievals"
base="${ROUND6_SOURCE_ROOT:-$HOME/code/github/uni_vSmartMOM_round5_jacobian}"
destination="${ROUND6_RELEASE_PARENT:-$bottom/round6_fixed_sif/deployment}"
release=round6_fixed_sif_standard_utls_v1
bundle="$destination/$release"
[[ ! -e "$bundle" && ! -e "$destination/$release.tar.gz" ]] || {
    echo "STOP: release already exists; do not change a deployed release" >&2; exit 1;
}
mkdir -p "$bundle/inversion/retrieval_setup" "$bundle/inputs"
for name in Round6FixedSIF.jl Round6RetrievalCampaign.jl run_round6_fixed_sif_retrievals.jl; do
    cp -p "$script_dir/$name" "$bundle/inversion/$name"
done
cp -p "$script_dir/retrieval_setup/build_round6_fixed_sif_apriori.jl" "$bundle/inversion/retrieval_setup/"
for name in round6_environment.sh run_round6_no_sif_partition.sh round6_sif_on_gattaca.sbatch submit_round6_sif_on_gattaca.sh; do
    cp -p "$script_dir/$name" "$bundle/$name"
done
cp -p "$base/Manifest.toml" "$bundle/inputs/Manifest.toml"
cp -p "$bottom/retrieval_setup/apriori_states_acos_mapped_tapered_vertical_correlation.nc" "$bundle/inputs/source_prior.nc"
for mode in off on; do
    prior="apriori_states_round6_fixed_sif_${mode}_standard_utls_acos_mapped_tapered_vertical_correlation"
    cp -p "$bottom/round6_fixed_sif/retrieval_setup/$prior.nc" "$bundle/inputs/"
    cp -p "$bottom/round6_fixed_sif/retrieval_setup/$prior.dat" "$bundle/inputs/"
done
cp -p "$HOME/RRS_XCO2_private/results/bottom_layer_sif_acos_mapped_tapered_vertical_correlation_v1/retrieval_setup/true_states_corrected_sif_v2.dat" "$bundle/inputs/"
cp -p "$bottom/truth/scene_components.dat" "$bundle/inputs/"
cp -p "$script_dir/instrument/representative_stokes_coefficients.nc" "$bundle/inputs/"
cp -p "$rrs_root/../src/SIF_emission/sif-spectra.csv" "$bundle/inputs/"
chmod 755 "$bundle/round6_environment.sh" "$bundle/run_round6_no_sif_partition.sh" \
    "$bundle/round6_sif_on_gattaca.sbatch" "$bundle/submit_round6_sif_on_gattaca.sh"
(
    cd "$bundle"
    find . -type f ! -name SHA256SUMS -print0 | LC_ALL=C sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum --check --strict SHA256SUMS
)
tar -czf "$destination/$release.tar.gz" -C "$destination" "$release"
(cd "$destination" && sha256sum "$release.tar.gz" > "$release.tar.gz.sha256")
echo "Released $destination/$release.tar.gz"
sha256sum "$destination/$release.tar.gz"
