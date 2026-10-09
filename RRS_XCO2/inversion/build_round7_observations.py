#!/usr/bin/env python3
"""Freeze Round-7 observations without modifying Round-5 products.

Delta = LUT(Cabannes + RRS - Rayleigh); y7 = y5_uncorrected - H(Delta).
The LUT is aerosol-free/SIF-off for EVERY scene. Pressure and geometry come
from the independent scene metadata, not the retrieved state. Only the three
O2 albedo coefficients use the mean of noisy uncorrected Round-5 retrievals.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
from datetime import datetime, timezone

import netCDF4
import numpy as np

REPO = Path(__file__).resolve().parents[2]
BOTTOM = REPO / "RRS_XCO2/bottom_layer_XCO2_retrievals"
PRIVATE = Path("/home/sanghavi/RRS_XCO2_private/results")
BASE_SHA = "7acab57000dae259207a6760faae156cfa1734f6"
CODESET_SHA = "1b2c4061bd9d7b3131523d2c76269cf103f0cc4f014904d69b4e9e161cc940f3"
SIF_ON = "angular_integral760_0p5"
NOISE_POLICY = "reuse_round5_uncorrected_noise_and_covariance"


def sha256(path):
    digest = hashlib.sha256()
    with open(str(path), "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite(values, label):
    require(not np.any(np.ma.getmaskarray(values)), label + " contains missing values")
    result = np.asarray(values, dtype=np.float64)
    require(np.all(np.isfinite(result)), label + " contains nonfinite values")
    require(np.all(np.abs(result) < 1e15), label + " contains fill/sentinel values")
    return result


def bracket(axis, query):
    """Linear weights on an increasing axis; exact nodes need no neighbours."""
    axis = np.asarray(axis, dtype=np.float64)
    require(np.all(np.isfinite(axis)) and np.all(np.diff(axis) > 0), "invalid LUT axis")
    require(np.isfinite(query) and axis[0] <= query <= axis[-1],
            "query %r is outside LUT [%r, %r]" % (query, axis[0], axis[-1]))
    high = int(np.searchsorted(axis, query))
    if high < len(axis) and axis[high] == query:
        return [(high, 1.0)]
    low = high - 1
    weight = (query - axis[low]) / (axis[high] - axis[low])
    return [(low, 1.0 - weight), (high, weight)]


def spectral_weights(axis, target):
    require(np.all(np.diff(axis) > 0), "spectral LUT axis must increase")
    require(np.min(target) >= axis[0] and np.max(target) <= axis[-1],
            "spectral extrapolation is forbidden")
    high = np.clip(np.searchsorted(axis, target, side="right"), 1, len(axis)-1)
    low = high - 1
    weight = (target - axis[low]) / (axis[high] - axis[low])
    return low, high, weight


def albedo_spectrum(coefficients, nu, basis_limits):
    """Use the original complete-band wavenumber coordinate, never LUT limits."""
    lo, hi = basis_limits
    require(hi > lo and np.shape(coefficients) == (3,), "invalid surface polynomial")
    x = 2.0 * (nu - lo) / (hi - lo) - 1.0
    return coefficients[0] + coefficients[1]*x + coefficients[2]*(3*x*x-1)/2


class RamanLUT:
    """Read only required pressure/SZA nodes from the trusted nadir LUTs."""
    def __init__(self, directory):
        self.directory = Path(directory)
        self.pressures = np.array([500.0, 750.0, 1000.0])
        self.cache = {}
        self.provenance = {}

    def _read(self, pressure, sza):
        path = self.directory / ("o2a_raman_lut_psurf%d.nc" % pressure)
        key = (float(pressure), float(sza))
        if key in self.cache:
            return self.cache[key]
        with netCDF4.Dataset(str(path)) as ds:
            require(np.array_equal(ds["psurf"][:], [pressure]), "LUT pressure mismatch")
            require(np.array_equal(ds["sif_on"][:], [0]), "Round7 requires the SIF-off LUT")
            require(float(ds.vza_deg) == 0 and float(ds.vaz_deg) == 0,
                    "unexpected nadir geometry convention")
            require(np.all(finite(ds["tau_aer"][:], "LUT aerosol") == 0),
                    "Round7 requires aerosol-free LUT")
            require(np.all(finite(ds["SIF0"][:], "LUT SIF") == 0), "nonzero LUT SIF")
            require(np.array_equal(ds["stokes_index"][:], [1, 2, 3]), "wrong Stokes order")
            nu = finite(ds["wn"][:], "LUT wavenumber")
            albedo = finite(ds["albedo"][:], "LUT albedo")
            # Equal-mu0 SZA grid: interpolate in cos(SZA), in Float64.
            mu = np.cos(np.deg2rad(finite(ds["sza"][:], "LUT SZA")))
            order = np.argsort(mu)
            nodes = [(int(order[i]), w) for i, w in bracket(mu[order], np.cos(np.deg2rad(sza)))]
            components = {}
            for name in ("rayleigh", "cabannes", "rrs"):
                variable = ds["stokes_" + name]
                require(variable.dimensions == ("wn", "stokes", "sza", "albedo", "sif", "psurf"),
                        "unexpected LUT dimension ordering")
                values = np.zeros((len(nu), 3, len(albedo)), dtype=np.float64)
                for i, w in nodes:
                    values += w * finite(variable[:, :, i, :, 0, 0], "LUT " + name)
                components[name] = values
            self.provenance[str(path.resolve())] = {
                "sha256": sha256(path), "pressure_hpa": pressure,
                "sza_nodes_deg": [float(ds["sza"][i]) for i, _ in nodes],
                "mu0_weights": [w for _, w in nodes],
                "profile_reduction": int(ds.profile_reduction),
                "git_commit": str(ds.git_commit),
                "profile_source_note": str(ds.profile_source_note),
            }
        self.cache[key] = (nu, albedo, components)
        return self.cache[key]

    def evaluate(self, nu, albedo, pressure, sza, vza, raz):
        require(vza == 0 and raz == 0, "this release supports the existing nadir/RAZ0 scenes only")
        require(np.shape(nu) == np.shape(albedo), "albedo/grid shape mismatch")
        result = {name: np.zeros((3, len(nu))) for name in ("rayleigh", "cabannes", "rrs")}
        for index, pressure_weight in bracket(self.pressures, pressure):
            grid, alb_grid, components = self._read(self.pressures[index], sza)
            lo, hi, w = spectral_weights(grid, nu)
            alo, ahi, aw = spectral_weights(alb_grid, albedo)
            for name, cube in components.items():
                # Evaluate every albedo-node spectrum at the target wavenumber,
                # then use the albedo appropriate to that SAME wavenumber.
                at_nu = (1-w[:, None, None])*cube[lo] + w[:, None, None]*cube[hi]
                row = np.arange(len(nu))
                selected = ((1-aw[:, None])*at_nu[row, :, alo] + aw[:, None]*at_nu[row, :, ahi])
                result[name] += pressure_weight * selected.T
        return result


def instrument_operator(nu, stokes, coefficients, targets, fwhm=0.04):
    """Exact port of pinned SyntheticOCO2.process_stokes_spectrum.

    Cross-language comparison to that unchanged Julia function is mandatory
    before retrieval launch (validate_round7_instrument.jl).
    """
    require(stokes.shape == (3, len(nu)), "invalid Stokes shape")
    wavelength = 1e7 / np.asarray(nu, dtype=np.float64)
    scalar = (coefficients[0]*stokes[0] - coefficients[1]*stokes[1] + coefficients[2]*stokes[2])
    spectrum = scalar * (1e7 / wavelength**2)
    order = np.argsort(wavelength)
    wavelength, spectrum = wavelength[order], spectrum[order]
    require(np.all(np.diff(wavelength) > 0), "duplicate spectral coordinates")
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    radius = 6*sigma
    require(wavelength[0] <= np.min(targets)-radius and wavelength[-1] >= np.max(targets)+radius,
            "missing six-sigma ILS shoulders")
    weights = np.empty(len(wavelength))
    weights[0] = (wavelength[1]-wavelength[0])/2
    weights[-1] = (wavelength[-1]-wavelength[-2])/2
    weights[1:-1] = (wavelength[2:]-wavelength[:-2])/2
    result = np.empty(len(targets))
    for i, center in enumerate(targets):
        lo, hi = np.searchsorted(wavelength, center-radius), np.searchsorted(wavelength, center+radius, side="right")
        require(hi > lo, "empty ILS kernel")
        kernel = np.exp(-0.5*((wavelength[lo:hi]-center)/sigma)**2)*weights[lo:hi]
        result[i] = np.sum(kernel*spectrum[lo:hi])/np.sum(kernel)
    return finite(result, "processed correction")


def read_truth_table(path):
    header = None
    result = []
    with open(str(path)) as stream:
        for line in stream:
            if line.startswith("# index "):
                header = line[2:].split()
            elif line.strip() and not line.startswith("#"):
                require(header is not None, "missing truth-table header")
                fields = line.split()
                require(len(fields) == len(header), "malformed truth row")
                result.append(dict(zip(header, fields)))
    require([int(row["index"]) for row in result] == list(range(1, 81)), "expected all 80 scenes")
    require(set(row["sif_case"] for row in result) == {"off", SIF_ON}, "stale SIF truth table")
    return result


def read_round5_scene(root, row, expected_prior_sha):
    index = int(row["index"])
    records, sources = [], []
    for perturbation in range(1, 12):
        path = root / "uncorrected" / ("retrieval_state%03d_perturbation%02d.nc" % (index, perturbation))
        with netCDF4.Dataset(str(path)) as ds:
            attrs = dict(ds.__dict__)
            expected = {"truth_state_index": index, "perturbation_index": perturbation,
                        "measurement_class": "uncorrected", "surface": row["surface"],
                        "aerosol_case": row["aerosol_case"], "sif_case": row["sif_case"],
                        "retrieval_complete": 1, "retrieval_state_model": "round5_fixed_sif",
                        "round5_code_checkpoint_sha": BASE_SHA, "round5_codeset_sha256": CODESET_SHA,
                        "round5_prior_sha256": expected_prior_sha, "state_dimension": 28,
                        "round5_stratospheric_aerosol_sigma_scale": 0.1}
            for key, value in expected.items():
                require(attrs.get(key) == value, "%s: unexpected %s" % (path, key))
            if perturbation <= 10:
                require(attrs.get("converged") == 1 and attrs.get("fit_quality_ok") == 1,
                        "every scene must have ten valid noisy uncorrected members: " + str(path))
            names = attrs["parameter_names"].split()
            require(names[19:22] == ["o2a_surface_P0", "o2a_surface_P1", "o2a_surface_P2"],
                    "unexpected O2 surface coordinates")
            variables = {name: finite(ds[name][:], str(path) + ":" + name) for name in (
                "final_state", "measurement_noiseless", "measurement_perturbed",
                "normalized_noise_draw", "injected_measurement_noise", "noise_standard_deviation",
                "Se_diagonal", "wavelength", "band_start_index", "band_end_index")}
            variables["attrs"] = attrs
            require(np.array_equal(variables["measurement_perturbed"],
                                   variables["measurement_noiseless"] + variables["injected_measurement_noise"]),
                    "stored Round5 measurement/noise identity failed")
            require(np.array_equal(variables["injected_measurement_noise"],
                                   variables["noise_standard_deviation"] * variables["normalized_noise_draw"]),
                    "stored Round5 noise-draw identity failed")
            require(np.array_equal(variables["Se_diagonal"], variables["noise_standard_deviation"]**2),
                    "stored Round5 covariance mismatch")
            require(np.all(variables["Se_diagonal"] > 0), "invalid covariance")
        records.append(variables)
        sources.append({"path": str(path.resolve()), "sha256": sha256(path), "perturbation": perturbation})
    ref = records[0]
    for record in records[1:]:
        for key in ("measurement_noiseless", "noise_standard_deviation", "Se_diagonal", "wavelength",
                    "band_start_index", "band_end_index"):
            require(np.array_equal(ref[key], record[key]), "inconsistent members: " + key)
    require(np.all(records[-1]["injected_measurement_noise"] == 0), "p11 must be noiseless")
    ideal_path = root / "corrected" / ("retrieval_state%03d_perturbation11.nc" % index)
    with netCDF4.Dataset(str(ideal_path)) as ds:
        require(ds.truth_state_index == index and ds.measurement_class == "corrected"
                and ds.sif_case == row["sif_case"] and ds.retrieval_complete == 1,
                "wrong ideal-corrected diagnostic source")
        ideal = finite(ds["measurement_noiseless"][:], "ideal corrected reference")
        require(np.array_equal(ds["wavelength"][:], ref["wavelength"]), "ideal-grid mismatch")
    sources.append({"path": str(ideal_path.resolve()), "sha256": sha256(ideal_path), "role": "ideal_reference_only"})
    return records, sources, ideal


def write_scene(path, row, records, ideal, correction, components, nu, albedo, coeffs,
                analyzer, shared_attributes, source_records, lut_provenance):
    ref = records[0]
    y = ref["measurement_noiseless"]
    y7 = y - correction
    starts, stops = ref["band_start_index"].astype(int), ref["band_end_index"].astype(int)
    require(np.array_equal(starts, [1, 935, 1742]) and np.array_equal(stops, [934, 1741, 2742]),
            "unexpected observation band grid")
    require(np.all(correction[934:] == 0) and np.array_equal(y7[934:], y[934:]), "CO2 bands changed")
    require(np.array_equal(ideal[934:], y[934:]), "Round5 CO2 observation classes differ")
    noise = np.stack([r["injected_measurement_noise"] for r in records])
    old_perturbed = np.stack([r["measurement_perturbed"] for r in records])
    # Use precisely the arithmetic order required by the existing result writer.
    new_perturbed = y7[None, :] + noise
    require(np.allclose(new_perturbed, old_perturbed-correction, rtol=0, atol=2e-12),
            "paired subtraction identity failed")
    with netCDF4.Dataset(str(path), "w", format="NETCDF4") as ds:
        for name, count in (("measurement", len(y)), ("perturbation", 11), ("band", 3),
                            ("coefficient", 3), ("noisy_member", 10), ("high_resolution", len(nu)),
                            ("stokes", 3), ("analyzer_coefficient", 4)):
            ds.createDimension(name, count)
        def variable(name, values, dimensions, units="1", dtype="f8"):
            v = ds.createVariable(name, dtype, dimensions, zlib=True, complevel=4)
            v.units = units
            v[:] = values
        radiance_units = "mW m-2 sr-1 nm-1"
        for name, values in (("measurement_uncorrected", y), ("measurement_ideal_corrected", ideal),
                             ("measurement_imperfectly_corrected", y7), ("correction_observation", correction),
                             ("noise_standard_deviation", ref["noise_standard_deviation"])):
            variable(name, values, ("measurement",), radiance_units)
        variable("wavelength", ref["wavelength"], ("measurement",), "nm")
        variable("Se_diagonal", ref["Se_diagonal"], ("measurement",), "(mW m-2 sr-1 nm-1)^2")
        for name in ("band_start_index", "band_end_index"):
            variable(name, ref[name].astype(int), ("band",), dtype="i4")
        variable("perturbation_index", np.arange(1, 12), ("perturbation",), dtype="i4")
        variable("random_seed_uint64", np.array([int(r["attrs"]["random_seed_uint64"]) for r in records], dtype=np.uint64),
                 ("perturbation",), dtype="u8")
        for name, values in (("uncorrected_perturbed", old_perturbed),
                             ("imperfectly_corrected_perturbed", new_perturbed),
                             ("injected_measurement_noise", noise),
                             ("normalized_noise_draw", np.stack([r["normalized_noise_draw"] for r in records]))):
            variable(name, values, ("perturbation", "measurement"),
                     "1" if name == "normalized_noise_draw" else radiance_units)
        variable("mean_o2a_surface_coefficients", coeffs, ("coefficient",))
        variable("round5_uncorrected_o2a_surface_coefficients", np.stack([r["final_state"][19:22] for r in records[:10]]),
                 ("noisy_member", "coefficient"))
        variable("correction_wavenumber", nu, ("high_resolution",), "cm-1")
        variable("correction_surface_albedo", albedo, ("high_resolution",))
        variable("o2a_analyzer_coefficients", analyzer, ("analyzer_coefficient",))
        for name, values in components.items():
            variable("lut_" + name, values, ("stokes", "high_resolution"), "mW m-2 sr-1 (cm-1)-1")
        delta = components["cabannes"] + components["rrs"] - components["rayleigh"]
        variable("correction_stokes", delta, ("stokes", "high_resolution"), "mW m-2 sr-1 (cm-1)-1")
        attrs = dict(shared_attributes)
        attrs.update({k: v for k, v in ref["attrs"].items() if k.startswith("sif_")})
        attrs.update({
            "state_index": int(row["index"]), "surface": row["surface"], "aerosol_case": row["aerosol_case"],
            "sif_case": row["sif_case"], "psurf_hpa": float(row["psurf_hpa"]),
            "sza_deg": float(row["sza_deg"]), "vza_deg": float(row["vza_deg"]),
            "relative_azimuth_deg": float(row["relative_azimuth_deg"]),
            "round5_prior_sha256": ref["attrs"]["round5_prior_sha256"],
            "round5_campaign_identity_sha256": ref["attrs"]["round5_campaign_identity_sha256"],
            "round7_round5_sources_json": json.dumps(source_records, sort_keys=True),
            "round7_lut_sources_json": json.dumps(lut_provenance, sort_keys=True),
            "round7_observation_complete": 1,
        })
        ds.setncatts(attrs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=BOTTOM / "round7_fixed_sif_imperfect_correction/observations")
    parser.add_argument("--truth-table", type=Path, default=PRIVATE / "bottom_layer_sif_acos_mapped_tapered_vertical_correlation_v1/retrieval_setup/true_states_corrected_sif_v2.dat")
    parser.add_argument("--nosif-root", type=Path, default=BOTTOM / "round5_fixed_sif/retrievals_nosif")
    parser.add_argument("--sif-root", type=Path, default=PRIVATE / "bottom_layer_round5_fixed_sif_on_tight_utls_acos_mapped_tapered_vertical_correlation_v1/retrievals")
    parser.add_argument("--lut-directory", type=Path, default=Path("/home/sanghavi/data/RamanSIFgrid/O2ABand"))
    parser.add_argument("--grid", type=Path, default=BOTTOM / "truth/sim_wavelength.nc")
    parser.add_argument("--analyzer", type=Path, default=REPO / "RRS_XCO2/inversion/instrument/representative_stokes_coefficients.nc")
    args = parser.parse_args()
    require(not args.output.exists(), "output exists; immutable observations must never be overwritten")
    rows = read_truth_table(args.truth_table)
    with netCDF4.Dataset(str(args.grid)) as ds:
        basis_nu = finite(ds["o2a_wavenumber"][:], "canonical O2 grid")
    require(np.all(np.diff(basis_nu) > 0) and len(basis_nu) == 2735, "unexpected canonical O2 grid")
    # Retain native LUT spectral coordinates: its 0.1 cm-1 grid has an offset
    # from the truth grid. Unnecessary spectral interpolation would smooth
    # narrow lines before the real instrument convolution.
    with netCDF4.Dataset(str(args.lut_directory / "o2a_raman_lut_psurf1000.nc")) as ds:
        native_nu = finite(ds["wn"][:], "native LUT wavenumber")
    nu = native_nu[(native_nu >= basis_nu[0]) & (native_nu <= basis_nu[-1])]
    with netCDF4.Dataset(str(args.analyzer)) as ds:
        require(ds.projection_convention == "M11*I - M12*Q + M13*U", "wrong analyzer convention")
        analyzer = finite(ds["representative_stokes_coefficients"][0, :], "O2 analyzer")
    prior_paths = {mode: BOTTOM / ("round5_fixed_sif/retrieval_setup/apriori_states_round5_fixed_sif_%s_tight_utls_acos_mapped_tapered_vertical_correlation.nc" % mode) for mode in ("off", "on")}
    prior_hashes = {mode: sha256(path) for mode, path in prior_paths.items()}
    lut = RamanLUT(args.lut_directory)
    shared = {
        "round7_definition_version": 1,
        "round7_correction_definition": "Cabannes+RRS-Rayleigh",
        "round7_correction_operation": "subtract",
        "round7_pressure_geometry_source": "independent_scene_metadata",
        "round7_noise_policy": NOISE_POLICY,
        "round7_albedo_estimator": "arithmetic mean of ten valid Round5 uncorrected noisy members 01:10; excludes p11",
        "round7_lut_policy": "same aerosol-free SIF-off LUT for every aerosol/SIF combination",
        "round7_interpolation": "linear in pressure_hPa, cos(SZA), spectral albedo; native LUT wavenumbers; no extrapolation; nadir RAZ0",
        "round7_spectral_grid": "native LUT nodes within canonical O2 band; lambda=1e7/Float64(wn); no pre-convolution spectral interpolation",
        "round7_lut_limitations": "aerosol-free, SIF-off, 12-layer AFGL atmosphere; scalar-albedo LUT evaluated with spectral albedo; historical mixed-commit promoted LUT, identified by whole-file checksum",
        "round7_lut_radiance_units": "mW m-2 sr-1 (cm-1)-1; inferred from physical solar-source recipe; no extra pi or M11 normalization",
        "round7_forward_code_checkpoint_sha": BASE_SHA,
        "source_truth_table": str(args.truth_table.resolve()), "source_truth_table_sha256": sha256(args.truth_table),
        "source_analyzer_sha256": sha256(args.analyzer), "source_grid_sha256": sha256(args.grid),
        "source_generator_sha256": sha256(Path(__file__)),
        "surface_basis_min_wavenumber": float(basis_nu[0]), "surface_basis_max_wavenumber": float(basis_nu[-1]),
        "instrument_operator": "pinned SyntheticOCO2: M11*I-M12*Q+M13*U; per-cm-1 to per-nm; Gaussian 0.04nm FWHM; +/-6sigma; direct sampling",
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".round7-observations-", dir=str(args.output.parent)))
    print("staging", staging, flush=True)
    audit = []
    for row in rows:
        mode = "off" if row["sif_case"] == "off" else "on"
        root = args.nosif_root if mode == "off" else args.sif_root
        records, sources, ideal = read_round5_scene(root, row, prior_hashes[mode])
        mean_coeff = np.mean(np.stack([r["final_state"][19:22] for r in records[:10]]), axis=0)
        albedo = albedo_spectrum(mean_coeff, nu, (basis_nu[0], basis_nu[-1]))
        components = lut.evaluate(nu, albedo, float(row["psurf_hpa"]), float(row["sza_deg"]),
                                  float(row["vza_deg"]), float(row["relative_azimuth_deg"]))
        delta = components["cabannes"] + components["rrs"] - components["rayleigh"]
        correction = np.zeros_like(records[0]["measurement_noiseless"])
        correction[:934] = instrument_operator(nu, delta, analyzer, records[0]["wavelength"][:934])
        path = staging / ("OCO2round7_%03d.nc" % int(row["index"]))
        write_scene(path, row, records, ideal, correction, components, nu, albedo, mean_coeff,
                    analyzer, shared, sources, lut.provenance)
        noise_std = records[0]["noise_standard_deviation"][:934]
        before = records[0]["measurement_noiseless"][:934] - ideal[:934]
        after = before - correction[:934]
        record = {"state": int(row["index"]), "surface": row["surface"], "aerosol": row["aerosol_case"],
                  "sif": row["sif_case"], "mean_coefficients": mean_coeff.tolist(),
                  "albedo_range": [float(np.min(albedo)), float(np.max(albedo))],
                  "correction_rms_sigma": float(np.sqrt(np.mean((correction[:934]/noise_std)**2))),
                  "before_rms_sigma": float(np.sqrt(np.mean((before/noise_std)**2))),
                  "after_rms_sigma": float(np.sqrt(np.mean((after/noise_std)**2))),
                  "source_retrievals": sources}
        audit.append(record)
        print("state=%03d SIF=%s aerosol=%s albedo=[%.5f,%.5f] residual RMS/noise %.4f -> %.4f" %
              (record["state"], mode, record["aerosol"], *record["albedo_range"], record["before_rms_sigma"], record["after_rms_sigma"]), flush=True)
    manifest = {"definition": shared, "lut_sources": lut.provenance, "prior_hashes": prior_hashes, "scenes": audit}
    with open(str(staging / "provenance.json"), "w") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    with open(str(staging / "SHA256SUMS"), "w") as stream:
        for path in sorted(staging.iterdir()):
            if path.name != "SHA256SUMS":
                stream.write(sha256(path) + "  " + path.name + "\n")
    os.rename(str(staging), str(args.output))
    print("Published", args.output, "manifest_sha256=" + sha256(args.output / "SHA256SUMS"), flush=True)


if __name__ == "__main__":
    main()
