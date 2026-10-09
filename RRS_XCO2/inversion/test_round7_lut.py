"""CPU-only RamanLUT regression tests; all NetCDF inputs are in-memory mocks.

Run from test/ with:
    python3 -B -m unittest discover -s ../RRS_XCO2/inversion -p test_round7_lut.py

Only file I/O and hashing are mocked: the production reader, validation,
interpolation, cache, and provenance paths run unchanged. The oracle is a
multi-affine function with cross terms, not another interpolation routine.
"""

from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

import build_round7_observations as builder


COMPONENTS = ("rayleigh", "cabannes", "rrs")
DIMS = ("wn", "stokes", "sza", "albedo", "sif", "psurf")


def analytic_radiance(component, stokes, nu, albedo, pressure, mu):
    """Exactly multilinear, including spectral-albedo and four-axis coupling.

    Signed, nonzero Q/U distinguish Stokes-axis mistakes. This is a numerical
    oracle, deliberately not a physical principal-plane radiance spectrum.
    """
    x = (np.asarray(nu) - 13000.0) / 10.0
    p = (pressure - 500.0) / 500.0
    scale = (1.0, -0.37, 0.19)[stokes] * (component + 1)
    return scale * (2.0 + x + 3.0*albedo + 5.0*mu + 7.0*p
                    + 11.0*x*albedo + 13.0*albedo*mu
                    + 17.0*x*p + 19.0*x*albedo*mu*p)


class MockVariable:
    def __init__(self, values, dimensions):
        self.values = np.ma.array(values, copy=True)
        self.dimensions = dimensions

    def __getitem__(self, key):
        return self.values[key]


class MockDataset:
    def __init__(self, pressure):
        self.vza_deg = 0.0
        self.vaz_deg = 0.0
        self.profile_reduction = 12
        self.git_commit = "synthetic-parent-not-a-uniform-solver-claim"
        self.profile_source_note = "synthetic multi-affine fixture"
        self.variables = {}
        self.add("psurf", [pressure], ("psurf",))
        self.add("sif_on", [0], ("sif",))
        self.add("stokes_index", [1, 2, 3], ("stokes",))
        # Nonuniform axes expose accidental index-based interpolation.
        self.add("wn", [13000., 13004., 13011., 13019.], ("wn",))
        self.add("albedo", [0., .2, .65, 1.], ("albedo",))
        self.add("sza", [0., 20., 60., 70.], ("sza",))
        mu = np.cos(np.deg2rad(self["sza"][:]))
        self.add("mu0", mu, ("sza",))
        self.add("tau_aer", np.zeros((4, 1)), ("wn", "psurf"))
        self.add("SIF0", np.zeros((4, 1)), ("wn", "sif"))
        for component, name in enumerate(COMPONENTS):
            values = np.empty((4, 3, 4, 4, 1, 1))
            for stokes in range(3):
                values[:, stokes, :, :, 0, 0] = analytic_radiance(
                    component, stokes, self["wn"][:][:, None, None],
                    self["albedo"][:][None, None, :], pressure,
                    mu[None, :, None])
            self.add("stokes_" + name, values, DIMS)

    def add(self, name, values, dimensions):
        self.variables[name] = MockVariable(values, dimensions)

    def __getitem__(self, key):
        return self.variables[key]

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False


class RamanLUTTests(unittest.TestCase):
    def setUp(self):
        self.datasets = {p: MockDataset(p) for p in (500, 750, 1000)}

        def open_dataset(path):
            pressure = int(Path(path).stem.rsplit("psurf", 1)[1])
            if pressure not in self.datasets:
                raise FileNotFoundError(str(path))
            return self.datasets[pressure]

        opener = patch.object(builder.netCDF4, "Dataset", side_effect=open_dataset)
        self.open_dataset = opener.start()
        self.addCleanup(opener.stop)
        hasher = patch.object(builder, "sha256", return_value="a" * 64)
        self.hash_file = hasher.start()
        self.addCleanup(hasher.stop)
        self.lut = builder.RamanLUT("/synthetic-round7-lut")

    def evaluate(self, nu=None, albedo=None, pressure=1000., sza=30., vza=0., raz=0.):
        nu = np.array([13000., 13002., 13011., 13018., 13019.]) if nu is None else np.asarray(nu)
        albedo = np.full(nu.shape, .3) if albedo is None else np.asarray(albedo)
        return self.lut.evaluate(nu, albedo, pressure, sza, vza, raz)

    def assert_oracle(self, nu, albedo, pressure=1000., sza=30.):
        actual = self.evaluate(nu, albedo, pressure, sza)
        self.assertEqual(set(actual), set(COMPONENTS))
        for component, name in enumerate(COMPONENTS):
            expected = np.array([analytic_radiance(
                component, s, nu, albedo, pressure, np.cos(np.deg2rad(sza)))
                for s in range(3)])
            self.assertEqual(actual[name].shape, (3, len(nu)))
            self.assertEqual(actual[name].dtype, np.dtype("float64"))
            np.testing.assert_allclose(actual[name], expected, rtol=3e-14, atol=3e-13)

    def test_every_exact_pressure_sza_albedo_and_spectral_node(self):
        for pressure, ds in self.datasets.items():
            for iz, sza in enumerate(ds["sza"][:]):
                for ia, albedo in enumerate(ds["albedo"][:]):
                    with self.subTest(pressure=pressure, sza=sza, albedo=albedo):
                        actual = self.evaluate(ds["wn"][:], np.full(4, albedo), pressure, sza)
                        for name in COMPONENTS:
                            expected = ds["stokes_" + name][:, :, iz, ia, 0, 0].T
                            np.testing.assert_array_equal(actual[name], expected)

    def test_spectral_interpolation_on_nonuniform_grid(self):
        self.assert_oracle(np.array([13001., 13007., 13017.]), np.full(3, .2), sza=20.)

    def test_different_albedo_for_each_spectral_sample(self):
        # Distinguishes pointwise spectral albedo from a scalar or outer product.
        self.assert_oracle(np.array([13001., 13007., 13017.]), np.array([.8, .1, .4]), sza=20.)

    def test_mu_interpolation_is_not_linear_in_sza_degrees(self):
        nu, albedo = np.array([13004.]), np.array([.2])
        self.assert_oracle(nu, albedo, sza=40.)
        actual = self.evaluate(nu, albedo, sza=40.)["rrs"]
        wrong = .5 * (self.evaluate(nu, albedo, sza=20.)["rrs"]
                      + self.evaluate(nu, albedo, sza=60.)["rrs"])
        self.assertGreater(np.max(np.abs(actual - wrong)), .01)

    def test_pressure_interpolation_in_both_intervals(self):
        for pressure in (625., 875.):
            with self.subTest(pressure=pressure):
                self.assert_oracle(np.array([13004.]), np.array([.65]), pressure, 20.)

    def test_joint_spectral_albedo_mu_and_pressure_interpolation(self):
        self.assert_oracle(np.array([13001., 13008., 13016.]), np.array([.11, .48, .91]), 812.5, 43.)

    def test_query_order_is_preserved(self):
        self.assert_oracle(np.array([13019., 13001., 13011., 13001.]), np.array([1., .3, .6, .8]))

    def test_native_node_subset_is_not_smoothed(self):
        ds = self.datasets[1000]
        # A spike must survive unchanged when selecting native coordinates.
        for name in COMPONENTS:
            ds["stokes_" + name].values[1, :, 1, 1, 0, 0] += 1000.
        actual = self.evaluate([13004., 13011.], [.2, .2], sza=20.)
        for name in COMPONENTS:
            np.testing.assert_array_equal(actual[name], ds["stokes_" + name][1:3, :, 1, 1, 0, 0].T)

    def test_exact_pressure_needs_no_other_pressure_files(self):
        self.datasets.pop(500)
        self.datasets.pop(750)
        self.assert_oracle(np.array([13004.]), np.array([.2]), sza=20.)
        self.assertEqual(self.open_dataset.call_count, 1)

    def test_exact_sza_does_not_read_unneeded_slabs(self):
        ds = self.datasets[1000]
        for name in COMPONENTS:
            for iz in (0, 2, 3):
                ds["stokes_" + name].values[:, :, iz, :, :, :] = np.nan
        self.assert_oracle(np.array([13004.]), np.array([.2]), sza=20.)

    def test_cache_does_not_freeze_albedo_or_query_grid(self):
        self.assert_oracle(np.array([13001.]), np.array([.15]))
        self.assert_oracle(np.array([13002., 13017.]), np.array([.5, .95]))
        self.assertEqual(self.open_dataset.call_count, 1)
        self.assertEqual(self.hash_file.call_count, 1)
        provenance = next(iter(self.lut.provenance.values()))
        self.assertEqual(provenance["sha256"], "a" * 64)
        self.assertEqual(provenance["pressure_hpa"], 1000.)
        self.assertAlmostEqual(sum(provenance["mu0_weights"]), 1.)

    def test_missing_required_pressure_file(self):
        self.datasets.pop(750)
        with self.assertRaises(FileNotFoundError):
            self.evaluate(pressure=875.)

    def test_missing_required_component(self):
        del self.datasets[1000].variables["stokes_rrs"]
        with self.assertRaises(KeyError):
            self.evaluate()

    def test_masked_required_interpolation_node(self):
        # SZA 30 brackets indices 1/2; albedo .3 brackets indices 1/2.
        self.datasets[1000]["stokes_rrs"].values[1, 0, 2, 2, 0, 0] = np.ma.masked
        with self.assertRaisesRegex(ValueError, "missing"):
            self.evaluate()

    def test_nonfinite_and_finite_fill_values_are_rejected(self):
        for name in COMPONENTS:
            variable = self.datasets[1000]["stokes_" + name]
            original = variable.values.copy()
            for value in (np.nan, np.inf, -np.inf, 9.96921e36, -1e20):
                with self.subTest(component=name, value=value):
                    variable.values[1, 0, 1, 1, 0, 0] = value
                    self.lut.cache.clear()
                    with self.assertRaises(ValueError):
                        self.evaluate()
                    variable.values = original.copy()

    def test_out_of_range_and_nonfinite_queries(self):
        cases = ({"pressure": p} for p in (499., 1001., np.nan, np.inf))
        cases = list(cases) + [{"sza": s} for s in (71., 90., np.nan)]
        cases += [{"nu": [n], "albedo": [.3]} for n in (12999., 13020., np.nan, np.inf)]
        cases += [{"nu": [13004.], "albedo": [a]} for a in (-.01, 1.01, np.nan, np.inf)]
        for kwargs in cases:
            with self.subTest(query=kwargs), self.assertRaises(ValueError):
                self.evaluate(**kwargs)

    def test_shape_mismatch_and_unsupported_geometry(self):
        for kwargs in ({"nu": [13004., 13011.], "albedo": [.3]},
                       {"vza": 1.}, {"raz": 180.}, {"vza": np.nan}, {"raz": np.nan}):
            with self.subTest(query=kwargs), self.assertRaises(ValueError):
                self.evaluate(**kwargs)

    def test_duplicate_or_reversed_interpolation_axes(self):
        ds = self.datasets[1000]
        for name in ("wn", "albedo", "sza"):
            original = ds[name].values.copy()
            for duplicate in (True, False):
                # SZA is intentionally sorted by the reader; descending or
                # ascending storage are valid, but duplicate cosines are not.
                if name == "sza" and not duplicate:
                    continue
                with self.subTest(axis=name, duplicate=duplicate):
                    ds[name].values = original.copy()
                    if duplicate:
                        ds[name].values[1] = ds[name].values[0]
                    else:
                        ds[name].values = original[::-1].copy()
                    self.lut.cache.clear()
                    with self.assertRaises(ValueError):
                        self.evaluate()
            ds[name].values = original

    def test_reader_rejects_wrong_physics_and_dimension_metadata(self):
        mutations = (
            ("psurf", [999.]), ("sif_on", [1]), ("stokes_index", [1, 3, 2]),
            ("tau_aer", np.ones((4, 1)) * .01), ("SIF0", np.ones((4, 1))),
        )
        ds = self.datasets[1000]
        for name, values in mutations:
            original = ds[name].values
            with self.subTest(variable=name):
                ds[name].values = np.ma.array(values)
                with self.assertRaises(ValueError):
                    self.evaluate()
                ds[name].values = original
        for attribute in ("vza_deg", "vaz_deg"):
            with self.subTest(attribute=attribute):
                setattr(ds, attribute, 1.)
                with self.assertRaises(ValueError):
                    self.evaluate()
                setattr(ds, attribute, 0.)
        ds["stokes_cabannes"].dimensions = ("stokes", "wn", "sza", "albedo", "sif", "psurf")
        with self.assertRaisesRegex(ValueError, "dimension ordering"):
            self.evaluate()

    def test_evaluation_does_not_modify_query_arrays(self):
        nu, albedo = np.array([13002., 13015.]), np.array([.3, .8])
        original_nu, original_albedo = nu.copy(), albedo.copy()
        nu.flags.writeable = albedo.flags.writeable = False
        self.evaluate(nu, albedo)
        np.testing.assert_array_equal(nu, original_nu)
        np.testing.assert_array_equal(albedo, original_albedo)


if __name__ == "__main__":
    unittest.main()
