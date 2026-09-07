# Release checks

Run from the repository root unless a command explicitly changes directory.
Use Julia 1.10–1.12 and Python 3.11 or later.

```sh
bash tools/install_taplo.sh /path/to/local/bin
/path/to/local/bin/taplo lint
/path/to/local/bin/taplo format --check
python -m pip install -r tools/requirements-checks.txt
python tools/check_schema.py
julia --project=test tools/audit_api.jl api-inventory.json
(cd test && julia --project=. ../tools/check_examples.jl)
(cd test && julia --project=. runtests.jl)
(cd test && julia --project=. local/gpu/runtests.jl)
(cd test && julia --project=. local/workflows/runtests.jl)
npm ci --prefix docs
npm audit --prefix docs
(cd docs && CI=false julia --project=. make.jl)
```

Taplo 0.10.0 is pinned by SHA-256 for Linux x86_64. On another platform,
install the same Taplo version with that platform's package manager. The
configuration formats maintained TOML and applies the local RT schema only to
scene locations. The Python check validates all discovered YAML/TOML RT scenes
and negative cases; Julia parser tests cover runtime semantics.

The API audit checks exported bindings and parses registered docstrings in
package modules and loaded CUDA/Metal extensions. It records unloaded extensions
explicitly. This is mechanical coverage; numerical contract tests and review
are still required to establish semantic accuracy.

The workflow harness creates temporary synthetic measurement/truth products.
Without external scientific inputs it runs portable OE/campaign contracts and
reports data-dependent tests as skipped. To run the extended suite, set
`CO2_COVARIANCE_FILE` to the study covariance file and `VSMARTMOM_SIF_DATA_DIR`
to the directory containing the original `sif-spectra.csv`. Inputs are read-only;
test outputs and the synthetic producer manifest stay in temporary directories.
`test_forward_state_mapping.jl` is a separate full-OCO integration test requiring
real spectroscopy, solar, and prior data; it is not part of this harness.

The npm lockfile overrides Vite to 6.4.3 and xmldom to 0.9.12 to address the
September 2026 dependency audit. Vite lies outside VitePress 1.6.4's declared
Vite range; keep the strict build and development/preview smoke checks when
updating either package. Do not remove the override without repeating the audit.
