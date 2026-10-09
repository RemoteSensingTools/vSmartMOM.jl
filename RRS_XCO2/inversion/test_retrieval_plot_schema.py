#!/usr/bin/env python3
"""Regression tests for legacy/round-4/round-5 plotting-state reconstruction."""

import unittest

import numpy as np

from retrieval_plot_schema import (
    CORE_SIF_REFERENCE_WAVELENGTH_NM,
    MSIF_NAME,
    ROUND4_BASE_NAMES,
    ROUND4_KNOWN_WAVELENGTH_NM,
    ROUND4_KNOWN_LNU759,
    ROUND4_SIF_ON_CASE,
    ROUND4_TRUTH_ANGULAR_INTEGRAL760,
    ROUND4_TRUTH_SIF760,
    ROUND4_TRUTH_MSIF,
    SIF760_NAME,
    RetrievalPlotSchema,
)


class FakeDataset:
    def __init__(self, attributes, active_core):
        self.attributes = dict(attributes)
        self.variables = {
            "active_core_parameter_index": np.asarray(active_core, dtype=int),
        }

    def ncattrs(self):
        return tuple(self.attributes)

    def getncattr(self, name):
        return self.attributes[name]

    def __getitem__(self, name):
        return self.variables[name]

    def filepath(self):
        return "synthetic-retrieval.nc"


def round4_dataset(mode):
    delta_nu = (
        1.0e7 / CORE_SIF_REFERENCE_WAVELENGTH_NM -
        1.0e7 / ROUND4_KNOWN_WAVELENGTH_NM
    )
    active_names = list(ROUND4_BASE_NAMES)
    active_core = list(range(1, 29))
    if mode == "on":
        active_names.append(MSIF_NAME)
        active_core.append(30)
    attributes = {
        "parameter_names": " ".join(active_names),
        "state_dimension": len(active_names),
        "retrieval_state_model": "round4_known_sif759",
        "round4_known_sif_enabled": 1,
        "round4_sif_case": mode,
        "sif_case": "off" if mode == "off" else ROUND4_SIF_ON_CASE,
        "round4_known_sif_wavelength_nm": ROUND4_KNOWN_WAVELENGTH_NM,
        "round4_core_sif_reference_wavelength_nm": (
            CORE_SIF_REFERENCE_WAVELENGTH_NM
        ),
        "round4_known_sif_Lnu_mW_m-2_sr-1_per_cm-1": (
            0.0 if mode == "off" else ROUND4_KNOWN_LNU759
        ),
        "round4_delta_nu_760_minus_759_cm-1": delta_nu,
    }
    return FakeDataset(attributes, active_core)


def round5_dataset(mode):
    anchor = ROUND4_KNOWN_LNU759 if mode == "on" else 0.0
    slope = ROUND4_TRUTH_MSIF if mode == "on" else 0.0
    delta_nu = 1.0e7 / 760.0 - 1.0e7 / 759.0
    return FakeDataset({
        "parameter_names": " ".join(ROUND4_BASE_NAMES),
        "state_dimension": 28,
        "retrieval_state_model": "round5_fixed_sif",
        "round5_sif_case": mode,
        "sif_case": "off" if mode == "off" else ROUND4_SIF_ON_CASE,
        "round5_known_sif_wavelength_nm": 759.0,
        "round5_known_sif_Lnu_mW_m-2_sr-1_per_cm-1": anchor,
        "round5_fixed_mSIF_mW_m-2_sr-1_per_cm-2": slope,
        "round5_core_SIF760_mW_m-2_sr-1_per_cm-1": anchor + delta_nu * slope,
        "round5_active_state_dimension": 28,
        "round5_core_state_dimension": 30,
        "round5_active_to_core_parameter_index": " ".join(map(str, range(1, 29))),
        "round5_active_to_full_parameter_index": " ".join(map(str, [1] + list(range(6, 33)))),
    }, range(1, 29))


class RetrievalPlotSchemaTests(unittest.TestCase):
    def test_sif_off_expands_fixed_zeros(self):
        schema = RetrievalPlotSchema.from_dataset(round4_dataset("off"))
        active = np.arange(28, dtype=float)
        expanded = schema.expand_absolute_state(active)
        self.assertEqual(expanded.shape, (30,))
        np.testing.assert_array_equal(expanded[:28], active)
        np.testing.assert_array_equal(expanded[28:], 0.0)
        np.testing.assert_array_equal(schema.expand_tangent(active)[28:], 0.0)

    def test_sif_on_absolute_and_tangent_maps_differ_by_intercept(self):
        schema = RetrievalPlotSchema.from_dataset(round4_dataset("on"))
        active = np.arange(29, dtype=float)
        active[-1] = -7.0e-7
        absolute = schema.expand_absolute_state(active)
        tangent = schema.expand_tangent(active)
        self.assertAlmostEqual(absolute[-1], active[-1])
        self.assertAlmostEqual(tangent[-1], active[-1])
        self.assertAlmostEqual(
            absolute[-2], schema.known_lnu + schema.delta_nu * active[-1]
        )
        self.assertAlmostEqual(tangent[-2], schema.delta_nu * active[-1])
        self.assertAlmostEqual(absolute[-2] - tangent[-2], schema.known_lnu)

    def test_round4_truth_uses_the_affine_retrieval_state(self):
        schema = RetrievalPlotSchema.from_dataset(round4_dataset("on"))
        row = {
            "sif_case": ROUND4_SIF_ON_CASE,
            # The nonlinear template value is deliberately distinct from the
            # affine retrieval-space value reconstructed below.
            "sif_angular_integral760": ROUND4_TRUTH_ANGULAR_INTEGRAL760,
            SIF760_NAME: ROUND4_TRUTH_SIF760,
            MSIF_NAME: ROUND4_TRUTH_MSIF,
        }
        truth = schema.truth_sif_coordinates(row)
        self.assertEqual(truth["SIF759"], ROUND4_KNOWN_LNU759)
        self.assertEqual(truth[MSIF_NAME], ROUND4_TRUTH_MSIF)
        self.assertAlmostEqual(
            truth[SIF760_NAME],
            ROUND4_KNOWN_LNU759 + schema.delta_nu * ROUND4_TRUTH_MSIF,
        )
        self.assertNotAlmostEqual(truth[SIF760_NAME], row[SIF760_NAME])

    def test_round4_accepts_ascii_rounded_truth_msif(self):
        schema = RetrievalPlotSchema.from_dataset(round4_dataset("on"))
        row = {
            "sif_case": ROUND4_SIF_ON_CASE,
            "sif_angular_integral760": ROUND4_TRUTH_ANGULAR_INTEGRAL760,
            SIF760_NAME: 4.596394756494e-3,
            MSIF_NAME: 1.229123068146e-5,
        }
        truth = schema.truth_sif_coordinates(row, "ASCII truth table")
        self.assertEqual(truth[MSIF_NAME], row[MSIF_NAME])

    def test_legacy_coordinates_are_unchanged(self):
        names = ROUND4_BASE_NAMES + (SIF760_NAME, MSIF_NAME)
        dataset = FakeDataset({
            "parameter_names": " ".join(names),
            "state_dimension": len(names),
        }, [])
        dataset.variables.pop("active_core_parameter_index")
        schema = RetrievalPlotSchema.from_dataset(dataset)
        state = np.linspace(-1.0, 1.0, len(names))
        np.testing.assert_array_equal(schema.expand_absolute_state(state), state)
        np.testing.assert_array_equal(schema.expand_tangent(state), state)

    def test_round4_rejects_legacy_sif_case_label(self):
        dataset = round4_dataset("on")
        dataset.attributes["sif_case"] = "total_0p5"
        with self.assertRaisesRegex(ValueError, "expected"):
            RetrievalPlotSchema.from_dataset(dataset)

    def test_round4_rejects_drifted_known_sif_anchor(self):
        dataset = round4_dataset("on")
        dataset.attributes[
            "round4_known_sif_Lnu_mW_m-2_sr-1_per_cm-1"
        ] *= 1.001
        with self.assertRaisesRegex(ValueError, "known_sif_Lnu"):
            RetrievalPlotSchema.from_dataset(dataset)


class Round5PlotSchemaTests(unittest.TestCase):
    def test_absolute_fixed_values_and_zero_tangents(self):
        for mode in ("off", "on"):
            with self.subTest(mode=mode):
                schema = RetrievalPlotSchema.from_dataset(round5_dataset(mode))
                self.assertTrue(schema.is_round5)
                self.assertFalse(schema.is_round4)
                active = np.arange(56, dtype=float).reshape(2, 28)
                absolute = schema.expand_absolute_state(active)
                tangent = schema.expand_tangent(active)
                self.assertEqual(absolute.shape, (2, 30))
                np.testing.assert_array_equal(absolute[:, :28], active)
                np.testing.assert_array_equal(tangent[:, :28], active)
                np.testing.assert_array_equal(tangent[:, 28:], 0.0)
                np.testing.assert_allclose(
                    absolute[:, -2], schema.known_lnu + schema.delta_nu * schema.fixed_msif,
                    rtol=0, atol=0,
                )
                np.testing.assert_array_equal(absolute[:, -1], schema.fixed_msif)
                native = schema.native_dict(active[0])
                self.assertEqual(native["SIF759"], schema.known_lnu)
                self.assertEqual(native[MSIF_NAME], schema.fixed_msif)
                native_tangent = schema.native_dict(active[0], tangent=True)
                for name in ("SIF759", SIF760_NAME, MSIF_NAME):
                    self.assertEqual(native_tangent[name], 0.0)
                self.assertIn("fixed", schema.msif_status)
                self.assertIn("round 5", schema.description())

    def test_truth_uses_fixed_slope_despite_ascii_rounding(self):
        schema = RetrievalPlotSchema.from_dataset(round5_dataset("on"))
        row = {
            "sif_case": ROUND4_SIF_ON_CASE,
            "sif_angular_integral760": 0.5,
            SIF760_NAME: 4.596394756494e-3,
            MSIF_NAME: 1.229123068146e-5,
        }
        truth = schema.truth_sif_coordinates(row)
        native = schema.native_dict(np.zeros(28))
        for name in ("SIF759", SIF760_NAME, MSIF_NAME):
            self.assertEqual(truth[name], native[name])

    def test_off_truth_is_zero(self):
        schema = RetrievalPlotSchema.from_dataset(round5_dataset("off"))
        row = {"sif_case": "off", SIF760_NAME: 0.0, MSIF_NAME: 0.0}
        self.assertEqual(schema.truth_sif_coordinates(row),
                         {"SIF759": 0.0, SIF760_NAME: 0.0, MSIF_NAME: 0.0})

    def test_rejects_wrong_metadata(self):
        changes = {
            "round5_sif_case": "legacy",
            "sif_case": "total_0p5",
            "state_dimension": 29,
            "parameter_names": " ".join(reversed(ROUND4_BASE_NAMES)),
            "round5_active_state_dimension": 29,
            "round5_core_state_dimension": 29,
            "round5_known_sif_wavelength_nm": 760.0,
            "round5_known_sif_Lnu_mW_m-2_sr-1_per_cm-1": 0.0,
            "round5_fixed_mSIF_mW_m-2_sr-1_per_cm-2": 0.0,
            "round5_core_SIF760_mW_m-2_sr-1_per_cm-1": ROUND4_TRUTH_SIF760,
            "round5_active_to_core_parameter_index": " ".join(map(str, range(2, 30))),
            "round5_active_to_full_parameter_index": " ".join(map(str, range(1, 29))),
        }
        for key, value in changes.items():
            with self.subTest(attribute=key):
                dataset = round5_dataset("on")
                dataset.attributes[key] = value
                with self.assertRaises(ValueError):
                    RetrievalPlotSchema.from_dataset(dataset)

    def test_rejects_missing_slope_and_incorrect_core_map(self):
        dataset = round5_dataset("on")
        del dataset.attributes["round5_fixed_mSIF_mW_m-2_sr-1_per_cm-2"]
        with self.assertRaisesRegex(ValueError, "missing"):
            RetrievalPlotSchema.from_dataset(dataset)
        dataset = round5_dataset("on")
        dataset.variables["active_core_parameter_index"][0] = 2
        with self.assertRaisesRegex(ValueError, "active_core_parameter_index"):
            RetrievalPlotSchema.from_dataset(dataset)

    def test_rejects_nonzero_off_coefficients_and_wrong_truth(self):
        for name in ("round5_known_sif_Lnu_mW_m-2_sr-1_per_cm-1",
                     "round5_fixed_mSIF_mW_m-2_sr-1_per_cm-2"):
            dataset = round5_dataset("off")
            dataset.attributes[name] = 1.0e-30
            with self.assertRaises(ValueError):
                RetrievalPlotSchema.from_dataset(dataset)
        schema = RetrievalPlotSchema.from_dataset(round5_dataset("on"))
        with self.assertRaises(ValueError):
            schema.validate_truth_row({"sif_case": "off"})


if __name__ == "__main__":
    unittest.main()
