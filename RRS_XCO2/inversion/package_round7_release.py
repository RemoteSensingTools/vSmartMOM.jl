#!/usr/bin/env python3
"""Package or independently verify Round 7; never build observations or launch Julia.

After approval: python3 package_round7_release.py pack
Read-only:     python3 <bundle>/package_round7_release.py verify --bundle <bundle>

The observation builder must first publish observations/SHA256SUMS containing
exactly 80 NC files and provenance.json. Packaging seals that manifest, freezes
the actual Round-5 small inputs, and refuses to replace an existing release.
Compatible with the hosts' Python 3.6; verification needs only the stdlib.
"""

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile


MODEL = "round7_fixed_sif_imperfect_correction"
RELEASE = MODEL + "_v1"
CHECKPOINT = "7acab57000dae259207a6760faae156cfa1734f6"
# SHA256 of the exact bytes of `git ls-tree -r --full-tree HEAD`.
BASE_TREE = "1b2c4061bd9d7b3131523d2c76269cf103f0cc4f014904d69b4e9e161cc940f3"
OBSERVATION_SHA = "a281bedd941037177a1958fa1455ce3e5f53674ccc1854a20796dfd2aafe8f47"
SOLAR_SHA = "9e44e424788ce9f91c654398a789cd9da80205b9d6dc0c00df5e20f500ee8644"
CONFIG_SHA = "30fac998ea635ed0bd994bab7b89d173042b24c715ba0e8e8b22e96e38647c81"
CONFIG_RELATIVE = "sandbox/workflows/RRS_XCO2/config/oco_grass_3aerosol.yaml"
PRIOR_NAMES = {mode: "apriori_states_round5_fixed_sif_%s_tight_utls_acos_mapped_tapered_vertical_correlation.nc" % mode
               for mode in ("off", "on")}
INPUT_HASHES = {
    "Manifest.toml": "a2981cc28cba51473ce5526131e7800285f76e45003ddfae7a283172080df1ae",
    PRIOR_NAMES["off"]: "37725edc467a9d5056489a4cb05a5823dfa8386b24556c7874d6ddf279defa76",
    PRIOR_NAMES["on"]: "578ab22d61badb1c8dbfce5f9e4d205025217b0a55b2bf644515437b16daab05",
    "source_prior.nc": "34c0e81b7a853af157b68bb879db0771342a9579460835675f92df8aab7f9375",
    "true_states_corrected_sif_v2.dat": "282592c1e9f22794216d16d331017b71cb130fb11bc4d2548a5d75c68e541913",
    "representative_stokes_coefficients.nc": "cd1ac380508f3e8b53bb759af5e35960af822af6ffe9cced3e6ac2e6a2b3dda0",
    "scene_components.dat": "65eb6bec7d3a0be6888d421e334d7ba1abdd7bb773c5740b1ba8668a4b53dff9",
    "sif-spectra.csv": "bfea28dc130f38c286a4a7c6091662b5df032685780e0291b5313fd022d60d77",
}
SCRIPTS = ("package_round7_release.py", "round7_environment.sh", "run_round7_partition.sh")
WORKFLOW = ("Round7RetrievalCampaign.jl", "run_round7_imperfect_correction_retrievals.jl",
            "build_round7_observations.py", "validate_round7_instrument.jl")
OBSERVATIONS = {"OCO2round7_%03d.nc" % index for index in range(1, 81)} | {"provenance.json"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def digest_lines(values):
    return hashlib.sha256(("\n".join(values) + "\n").encode("utf-8")).hexdigest()


def check_hash(path, expected):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), "missing or symlink file: " + str(path))
    require(sha256(path) == expected, "checksum mismatch: " + str(path))


def parse_manifest(text):
    entries = {}
    for line in text.splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  ([A-Za-z0-9_./-]+)", line)
        require(match is not None, "invalid SHA256SUMS entry: " + line)
        digest, name = match.groups()
        path = PurePosixPath(name)
        require(not path.is_absolute() and all(p not in ("", ".", "..") for p in name.split("/")),
                "unsafe manifest path: " + name)
        require(name not in entries and name != "SHA256SUMS", "duplicate/self manifest entry: " + name)
        entries[name] = digest
    require(entries, "empty SHA256SUMS")
    return entries


def verify_manifest(root, expected_names=None):
    root = Path(root)
    require(not root.is_symlink(), "symlink manifest directory: " + str(root))
    manifest = root / "SHA256SUMS"
    require(manifest.is_file() and not manifest.is_symlink(), "missing/symlink manifest: " + str(manifest))
    entries = parse_manifest(manifest.read_text(encoding="utf-8"))
    if expected_names is not None:
        require(set(entries) == set(expected_names), "manifest inventory differs: " + str(root))
    actual = set()
    for path in root.rglob("*"):
        require(not path.is_symlink(), "symlink in immutable tree: " + str(path))
        if path.is_file() and path != manifest:
            actual.add(path.relative_to(root).as_posix())
        else:
            require(path == manifest or path.is_dir(), "nonregular release entry: " + str(path))
    require(actual == set(entries), "unlisted or missing files: " + str(root))
    for name, digest in entries.items():
        check_hash(root / name, digest)
    return sha256(manifest)


def verify_base(root):
    root = Path(root).resolve()
    def git(*arguments):
        return subprocess.check_output(["git", "--no-optional-locks", "-C", str(root)] + list(arguments))
    require(git("rev-parse", "HEAD").decode().strip() == CHECKPOINT, "wrong pinned checkpoint")
    require(not git("status", "--porcelain", "--untracked-files=all"), "pinned checkout is dirty")
    require(hashlib.sha256(git("ls-tree", "-r", "--full-tree", "HEAD")).hexdigest() == BASE_TREE,
            "wrong pinned tree")
    require(not (root / "LocalPreferences.toml").exists(), "unsealed base LocalPreferences.toml")
    check_hash(root / "Manifest.toml", INPUT_HASHES["Manifest.toml"])
    check_hash(root / CONFIG_RELATIVE, CONFIG_SHA)
    return root


def identity_code_files(campaign):
    """Read the Julia literal contract without executing any Julia code."""
    text = Path(campaign).read_text(encoding="utf-8")
    match = re.search(r"const ROUND7_CODE_FILES\s*=\s*\(([^)]*)\)", text)
    require(match is not None, "missing literal ROUND7_CODE_FILES contract")
    names = re.findall(r'"([A-Za-z0-9_.-]+)"', match.group(1))
    remainder = re.sub(r'"[A-Za-z0-9_.-]+"|[\s,]', "", match.group(1))
    require(names and not remainder and len(names) == len(set(names)), "unsupported ROUND7_CODE_FILES contract")
    require(set(WORKFLOW[:2]).issubset(names), "identity must cover both runtime Julia files")
    return names


def verify_observations(directory, generator):
    manifest_sha = verify_manifest(directory, OBSERVATIONS)
    require(manifest_sha == OBSERVATION_SHA, "observations differ from the approved published manifest")
    record = json.loads((directory / "provenance.json").read_text(encoding="utf-8"))
    definition = record["definition"]
    expected = {
        "round7_definition_version": 1,
        "round7_correction_definition": "Cabannes+RRS-Rayleigh",
        "round7_correction_operation": "subtract",
        "round7_pressure_geometry_source": "independent_scene_metadata",
        "round7_noise_policy": "reuse_round5_uncorrected_noise_and_covariance",
        "round7_forward_code_checkpoint_sha": CHECKPOINT,
        "source_truth_table_sha256": INPUT_HASHES["true_states_corrected_sif_v2.dat"],
        "source_analyzer_sha256": INPUT_HASHES["representative_stokes_coefficients.nc"],
        "source_generator_sha256": sha256(generator),
    }
    for key, value in expected.items():
        require(definition.get(key) == value, "observation provenance mismatch: " + key)
    require(record["prior_hashes"] == {mode: INPUT_HASHES[name] for mode, name in PRIOR_NAMES.items()},
            "observations did not use the actual Round-5 priors")
    require(sorted(scene["state"] for scene in record["scenes"]) == list(range(1, 81)),
            "observation provenance must enumerate each of the 80 scenes once")
    return manifest_sha


def campaign_root(path):
    path = Path(path).resolve()
    require(path.name == MODEL, "Round-7 root must end in " + MODEL)
    return path


def verify_release(bundle, base, root, solar):
    bundle = Path(bundle).resolve()
    require(bundle.name == RELEASE, "wrong release directory name")
    release_sha = verify_manifest(bundle)
    metadata = json.loads((bundle / "release.json").read_text(encoding="utf-8"))
    require(metadata["release"] == RELEASE and metadata["checkpoint"] == CHECKPOINT and
            metadata["base_tree_sha256"] == BASE_TREE, "wrong release metadata")
    code_files = identity_code_files(bundle / "inversion/Round7RetrievalCampaign.jl")
    require(metadata["identity_code_files"] == code_files, "Julia identity contract changed")
    for name in SCRIPTS:
        require((bundle / name).is_file(), "missing deployment script: " + name)
    for name in set(WORKFLOW) | set(code_files):
        require((bundle / "inversion" / name).is_file(), "missing workflow file: " + name)
    for name, digest in INPUT_HASHES.items():
        check_hash(bundle / "inputs" / name, digest)
    base = verify_base(base)
    check_hash(solar, SOLAR_SHA)
    observations = campaign_root(root) / "observations"
    observation_sha = verify_observations(observations, bundle / "inversion/build_round7_observations.py")
    require((bundle / "inputs/round7_observation_manifest.sha256").read_text().strip() == observation_sha,
            "observation manifest differs from the sealed release")
    check_hash(bundle / "inputs/OBSERVATION_SHA256SUMS", observation_sha)
    check_hash(bundle / "inputs/round7_observation_provenance.json", sha256(observations / "provenance.json"))
    require(metadata["observation_manifest_sha256"] == observation_sha, "release observation seal differs")
    return base, observations, code_files, release_sha


def identity(bundle, base, observations, solar, mode, code_files):
    inputs = bundle / "inputs"
    ordered = [inputs / "true_states_corrected_sif_v2.dat", inputs / PRIOR_NAMES[mode],
               inputs / "source_prior.nc", inputs / "representative_stokes_coefficients.nc",
               inputs / "scene_components.dat", inputs / "sif-spectra.csv", solar,
               base / CONFIG_RELATIVE, observations / "SHA256SUMS"]
    codeset = digest_lines([BASE_TREE] + [sha256(bundle / "inversion" / name) for name in code_files])
    input_set = digest_lines([sha256(path) for path in ordered])
    observation_sha = sha256(observations / "SHA256SUMS")
    return {
        "ROUND7_CODE_CHECKPOINT_SHA": CHECKPOINT,
        "ROUND7_CODESET_SHA256": codeset,
        "ROUND7_INPUT_SET_SHA256": input_set,
        "ROUND7_CAMPAIGN_IDENTITY_SHA256": digest_lines([MODEL, CHECKPOINT, codeset, input_set, observation_sha]),
        "ROUND7_OBSERVATION_MANIFEST_SHA256": observation_sha,
    }


def write_manifest(root):
    paths = sorted(path for path in root.rglob("*") if path.is_file())
    with (root / "SHA256SUMS").open("x", encoding="utf-8") as stream:
        for path in paths:
            stream.write("%s  %s\n" % (sha256(path), path.relative_to(root).as_posix()))


def pack(args):
    workflow = Path(__file__).resolve().parent
    base = verify_base(args.source_root)
    root = campaign_root(args.round7_root)
    data = args.data_checkout.resolve() / "RRS_XCO2"
    check_hash(args.solar_out, SOLAR_SHA)
    code_files = identity_code_files(workflow / "Round7RetrievalCampaign.jl")
    tests = sorted(set(workflow.glob("test_round7*.py")) | set(workflow.glob("test_round7*.jl")))
    require(any(p.suffix == ".py" for p in tests) and any(p.suffix == ".jl" for p in tests),
            "both Python and Julia Round-7 tests must exist before packaging")
    sources = {name: workflow / name for name in SCRIPTS}
    sources.update({"inversion/" + name: workflow / name for name in set(WORKFLOW) | set(code_files)})
    sources.update({"inversion/" + p.name: p for p in tests})
    sources["ROUND7_IMPERFECT_CORRECTION.md"] = workflow / "ROUND7_IMPERFECT_CORRECTION.md"
    small_inputs = {
        "Manifest.toml": base / "Manifest.toml",
        "source_prior.nc": args.source_prior,
        "true_states_corrected_sif_v2.dat": args.truth_table,
        "representative_stokes_coefficients.nc": data / "inversion/instrument/representative_stokes_coefficients.nc",
        "scene_components.dat": data / "bottom_layer_XCO2_retrievals/truth/scene_components.dat",
        "sif-spectra.csv": args.data_checkout.resolve() / "src/SIF_emission/sif-spectra.csv",
    }
    small_inputs.update({name: data / "bottom_layer_XCO2_retrievals/round5_fixed_sif/retrieval_setup" / name
                         for name in PRIOR_NAMES.values()})
    for name, path in small_inputs.items():
        check_hash(path, INPUT_HASHES[name])
    sources.update({"inputs/" + name: path for name, path in small_inputs.items()})
    # Read-only control paths and their hashes are preserved in this record;
    # packaging never creates writable links back into the Round-5 results.
    sources["inputs/round7_observation_provenance.json"] = root / "observations/provenance.json"
    for path in sources.values():
        require(path.is_file() and not path.is_symlink(), "missing/symlink package source: " + str(path))
    # Snapshot digests before copying: concurrent agent edits cannot silently
    # change the builder whose digest the observations certify.
    source_hashes = {name: sha256(path) for name, path in sources.items()}
    observation_sha = verify_observations(root / "observations", workflow / "build_round7_observations.py")
    parent = root / "deployment"
    parent.mkdir(parents=True, exist_ok=True)
    with (parent / ".package.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        bundle = parent / RELEASE
        archive = parent / (RELEASE + ".tar.gz")
        checksum = parent / (RELEASE + ".tar.gz.sha256")
        require(not any(os.path.lexists(str(p)) for p in (bundle, archive, checksum)),
                "release already exists; never overwrite a deployed release")
        staging_parent = Path(tempfile.mkdtemp(prefix=".round7-package-", dir=str(parent)))
        staging = staging_parent / RELEASE
        staging.mkdir()
        try:
            for name, source in sources.items():
                target = staging / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(str(source), str(target))
                check_hash(target, source_hashes[name])
            shutil.copyfile(str(root / "observations/SHA256SUMS"), str(staging / "inputs/OBSERVATION_SHA256SUMS"))
            check_hash(staging / "inputs/OBSERVATION_SHA256SUMS", observation_sha)
            (staging / "inputs/round7_observation_manifest.sha256").write_text(observation_sha + "\n")
            metadata = {"release": RELEASE, "checkpoint": CHECKPOINT, "base_tree_sha256": BASE_TREE,
                        "identity_code_files": code_files, "observation_manifest_sha256": observation_sha}
            (staging / "release.json").write_text(json.dumps(metadata, sort_keys=True, indent=2) + "\n")
            write_manifest(staging)
            verify_release(staging, base, root, args.solar_out)
            temporary_archive = staging_parent / archive.name
            with tarfile.open(str(temporary_archive), "w:gz") as tar:
                tar.add(str(staging), arcname=RELEASE)
            archive_sha = sha256(temporary_archive)
            # Publication is under the packaging lock; files use non-overwriting
            # hardlinks. On failure retain evidence, never remove user data.
            staging.rename(bundle)
            os.link(str(temporary_archive), str(archive))
            with checksum.open("x") as stream:
                stream.write(archive_sha + "  " + archive.name + "\n")
            # Apply the immutable release permissions only after publication.
            # Some filesystems reject renaming a directory whose write bits were
            # removed while it is still inside the staging parent.
            for path in bundle.rglob("*"):
                path.chmod(0o500 if path.is_dir() or path.name in SCRIPTS else 0o400)
            bundle.chmod(0o500)
            archive.chmod(0o400)
            checksum.chmod(0o400)
            print("Published %s\nArchive SHA256: %s\nStaging evidence: %s" % (bundle, archive_sha, staging_parent))
        except Exception:
            print("Packaging stopped; retained staging: " + str(staging_parent), file=sys.stderr)
            raise


def main():
    home_dir = Path.home()
    checkout = home_dir / "code/github/uni_vSmartMOM"
    bottom = checkout / "RRS_XCO2/bottom_layer_XCO2_retrievals"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("pack", "verify"))
    parser.add_argument("--source-root", type=Path, default=home_dir / "code/github/uni_vSmartMOM_round5_jacobian")
    parser.add_argument("--data-checkout", type=Path, default=checkout)
    parser.add_argument("--round7-root", type=Path)
    parser.add_argument("--bundle", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--solar-out", type=Path, default=home_dir / "Raman_misc/workflows/worktrees/uni_vSmartMOM_sanghavi_2025-03-18/src/SolarModel/solar.out")
    parser.add_argument("--truth-table", type=Path, default=home_dir / "RRS_XCO2_private/results/bottom_layer_sif_acos_mapped_tapered_vertical_correlation_v1/retrieval_setup/true_states_corrected_sif_v2.dat")
    parser.add_argument("--source-prior", type=Path, default=bottom / "retrieval_setup/apriori_states_acos_mapped_tapered_vertical_correlation.nc")
    parser.add_argument("--mode", choices=("off", "on"), help="verify: print validated runner identity assignments")
    args = parser.parse_args()
    if args.round7_root is None:
        args.round7_root = args.data_checkout / "RRS_XCO2/bottom_layer_XCO2_retrievals" / MODEL
    try:
        if args.action == "pack":
            require(args.mode is None, "--mode is only for read-only verification")
            pack(args)
        else:
            bundle = args.bundle.resolve()
            base, observations, code_files, release_sha = verify_release(bundle, args.source_root, args.round7_root, args.solar_out)
            print("Verified release SHA256SUMS=" + release_sha, file=sys.stderr)
            if args.mode:
                for key, value in sorted(identity(bundle, base, observations, args.solar_out, args.mode, code_files).items()):
                    print(key + "=" + value)
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as exception:
        parser.exit(1, "STOP: %s\n" % exception)


if __name__ == "__main__":
    main()
