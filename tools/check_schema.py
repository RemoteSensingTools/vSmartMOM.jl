#!/usr/bin/env python3
"""Validate scene files and exercise negative cases with the shipped schema.

Run with Python 3.11+ and tools/requirements-checks.txt installed. This checks
static configuration structure; Julia tests separately check parser semantics.
"""
from copy import deepcopy
import json
from pathlib import Path
import tomllib

from jsonschema import Draft4Validator
import yaml

ROOT = Path(__file__).resolve().parents[1]
schema = json.loads((ROOT / "schemas/vsmartmom-parameters.schema.json").read_text())
Draft4Validator.check_schema(schema)
validator = Draft4Validator(schema)
patterns = ("config/**/*", "configs/**/*", "examples/**/*", "sandbox/**/*",
            "test/**/*", "src/CoreRT/DefaultParameters.yaml")
paths = sorted({path for pattern in patterns for path in ROOT.glob(pattern)
                if path.is_file() and path.suffix in (".yaml", ".yml", ".toml")
                and path.name not in ("Project.toml", "Manifest.toml")})
failures = []
scene_count = 0
for path in paths:
    data = (tomllib.loads(path.read_text()) if path.suffix == ".toml"
            else yaml.safe_load(path.read_text()))
    # Standalone aerosol input metadata and benchmark harness TOML are not RT scenes.
    if not isinstance(data, dict) or "radiative_transfer" not in data:
        continue
    scene_count += 1
    for error in validator.iter_errors(data):
        failures.append(f"{path.relative_to(ROOT)}:{'.'.join(map(str,error.path))}: {error.message}")

base = tomllib.loads((ROOT / "config/quickstart.toml").read_text())
cases = (
    ({"fourier_convergence": "intensity", "fourier_tolerance": 1e-5,
      "fourier_min_m": 3, "fourier_n_consecutive": 2}, True),
    ({"fourier_convergence": "stokes", "ss_correction": "tms"}, True),
    ({"fourier_convergence": "iqu", "dtau_max_threshold": 0.001}, True),
    ({"fourier_convergence": "all", "fourier_tolerance": 1e-5}, False),
    ({"fourier_tolerance": 1e-5}, False),
    ({"fourier_convergence": "intensity", "fourier_tolerance": 0}, False),
    ({"fourier_convergence": "stokes", "fourier_min_m": 2}, False),
    ({"fourier_convergence": "stokes", "fourier_n_consecutive": 0}, False),
    ({"fourier_convergence": "unknown"}, False),
    ({"ss_correction": "unknown"}, False),
    ({"fourier_tolernace": 1e-5}, False),
)
for numerics, valid in cases:
    data = deepcopy(base)
    data["radiative_transfer"]["numerics"] = numerics
    if validator.is_valid(data) != valid:
        failures.append(f"Numerics contract case failed: {numerics}, expected valid={valid}")
if failures:
    raise SystemExit("\n".join(failures))
print(f"Schema valid; {scene_count} scene files and {len(cases)} contract cases passed.")
