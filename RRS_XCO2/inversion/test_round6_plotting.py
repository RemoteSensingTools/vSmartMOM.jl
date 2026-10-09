"""Round-6 fixed-SIF reconstruction, provenance, routing, and plot labels."""

from pathlib import Path
import runpy
import unittest

import numpy as np

import plot_physical_retrieval_ensembles as plots
import test_physical_ensemble_round5 as physical_fixtures
from retrieval_plot_schema import RetrievalPlotSchema
from test_retrieval_plot_schema import round5_dataset


def round6_dataset(mode):
    dataset = round5_dataset(mode)
    dataset.attributes = {
        key.replace("round5_", "round6_"): value
        for key, value in dataset.attributes.items()
    }
    dataset.attributes["retrieval_state_model"] = "round6_fixed_sif"
    dataset.attributes["round6_stratospheric_aerosol_sigma_scale"] = 1.0
    return dataset


class Round6PlottingTests(unittest.TestCase):
    def test_absolute_and_tangent_reconstruction(self):
        for mode in ("off", "on"):
            with self.subTest(mode=mode):
                schema = RetrievalPlotSchema.from_dataset(round6_dataset(mode))
                self.assertTrue(schema.is_round6)
                self.assertFalse(schema.is_round5)
                self.assertTrue(schema.is_fixed_sif_model)
                self.assertIn("round 6:", schema.description())
                active = np.arange(56).reshape(2, 28)
                absolute = schema.expand_absolute_state(active)
                tangent = schema.expand_tangent(active)
                np.testing.assert_array_equal(absolute[:, :28], active)
                np.testing.assert_array_equal(tangent[:, :28], active)
                np.testing.assert_array_equal(tangent[:, 28:], 0.0)
                np.testing.assert_array_equal(absolute[:, -1], schema.fixed_msif)
                np.testing.assert_array_equal(
                    absolute[:, -2], schema.known_lnu + schema.delta_nu * schema.fixed_msif)
                round5 = RetrievalPlotSchema.from_dataset(round5_dataset(mode))
                self.assertNotEqual(schema.identity, round5.identity)
                np.testing.assert_array_equal(absolute, round5.expand_absolute_state(active))

    def test_rejects_wrong_or_missing_round6_metadata(self):
        for key, wrong in (
                ("round6_stratospheric_aerosol_sigma_scale", 0.1),
                ("round6_fixed_mSIF_mW_m-2_sr-1_per_cm-2", 0.0),
                ("round6_known_sif_wavelength_nm", 760.0),
                ("round6_active_state_dimension", 29),
                ("round6_core_state_dimension", 28),
                ("round6_active_to_full_parameter_index", "1 2"),
                ("round6_active_to_core_parameter_index", "1 2"),
                ("sif_case", "total_0p5")):
            for missing in (False, True):
                with self.subTest(key=key, missing=missing):
                    dataset = round6_dataset("on")
                    if missing:
                        del dataset.attributes[key]
                    else:
                        dataset.attributes[key] = wrong
                    with self.assertRaises(ValueError):
                        RetrievalPlotSchema.from_dataset(dataset)

    def test_physical_truth_and_fixed_sif_panel(self):
        fixture = physical_fixtures.Round5PhysicalEnsembleTests()
        for mode in ("off", "on"):
            with self.subTest(mode=mode):
                schema, record, truth = fixture.records(mode, round6_dataset)
                self.assertEqual(truth["sif"]["state_model"], "round6_fixed_sif")
                np.testing.assert_array_equal(plots.sif_curve(record), plots.sif_curve(truth))
                self.assertEqual(plots.sif_known_anchor_radiance(record),
                                 plots.sif_known_anchor_radiance(truth))
                ensemble = {
                    "plot_schema": schema, "indices": [1, 2, 3],
                    "corrected": [record] * 3, "uncorrected": [record] * 3,
                    "noiseless_paired": False,
                }
                figure, (axis, annotation_axis) = plots.plt.subplots(2, 1)
                try:
                    plots.draw_sif_card(
                        axis, annotation_axis,
                        {"ensemble": ensemble, "truth": truth, "state_index": 1},
                        (-0.1, 0.1), False, False)
                    self.assertIn("fixed", axis.get_title())
                    self.assertNotIn("retrieved slope", axis.get_title())
                    if mode == "on":
                        self.assertIn("no SIF retrieval spread",
                                      annotation_axis.texts[0].get_text())
                finally:
                    plots.plt.close(figure)
        handles, _ = plots.co2_figure_legend_groups(round4=True, round6=True)
        self.assertIn("round-6 fixed state", handles[0].get_label())

    def test_campaign_routing_preserves_round_identity(self):
        scripts = Path(__file__).resolve().parents[1] / "visualization"
        frontend = runpy.run_path(str(scripts / "plot_physical_retrieval_ensembles.py"))
        comparison = runpy.run_path(str(scripts / "plot_round3_vs_round4_xco2_psurf.py"))
        for sif in (0, 1):
            regime = frontend["resolve_regime"](6, sif)
            self.assertEqual(regime.state_model, "round6_fixed_sif")
            self.assertIn("standard_utls", regime.source_prior_basename)
            self.assertIn("original UTLS", regime.label)
            campaigns = comparison["COMPARISON_CAMPAIGNS"][6][sif]
            self.assertEqual(campaigns["Round 6"][0], regime.retrieval_root)
            self.assertEqual(campaigns["Round 6"][1], regime.truth_table)
            self.assertNotEqual(regime.output_root,
                                frontend["resolve_regime"](5, sif).output_root)


if __name__ == "__main__":
    unittest.main()
