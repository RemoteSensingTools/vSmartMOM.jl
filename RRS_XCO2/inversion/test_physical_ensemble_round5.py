"""Round-5 ensemble routing, physical SIF reconstruction, and label regressions."""

from pathlib import Path
import runpy
import unittest

import numpy as np

import plot_physical_retrieval_ensembles as plots
from retrieval_plot_schema import RetrievalPlotSchema
from test_retrieval_plot_schema import round5_dataset


class Round5PhysicalEnsembleTests(unittest.TestCase):
    def records(self, mode, dataset_factory=round5_dataset):
        schema = RetrievalPlotSchema.from_dataset(dataset_factory(mode))
        state = np.zeros(28)
        retrieved = plots.native_to_physical(schema, state, 400.0)
        row = {
            "index": "1", "aerosol_case": "none", "xco2_ppm": "400",
            "psurf_hpa": "1000", "sif_angular_integral760": 0.0,
            "sif_case": "off", "SIF760": 0.0, "mSIF": 0.0,
        }
        if mode == "on":
            row.update({
                "sif_case": plots.SIF_CASE_ON,
                "sif_angular_integral760": 0.5,
                "SIF760": 4.596394756494e-3,
                "mSIF": 1.229123068146e-5,
            })
        for _, truth_band, _, _, _ in plots.BANDS:
            for order in range(3):
                row["%s_P%d" % (truth_band, order)] = "0"
        aerosol_cases = {"none": {
            "sulfate_AOD760": 0.0, "organic_AOD760": 0.0,
            "utls_sulfate_AOD760": 0.0,
        }}
        vertical = {species: {"z0": 1.0}
                    for species, _, _, _, _ in plots.SPECIES}
        truth = plots.truth_record(row, aerosol_cases, vertical,
                                   plot_schema=schema)
        return schema, retrieved, truth

    def test_fixed_sif_truth_and_retrieval_curves_agree(self):
        for mode in ("off", "on"):
            with self.subTest(mode=mode):
                schema, retrieved, truth = self.records(mode)
                for record in (retrieved, truth):
                    self.assertEqual(record["sif"]["state_model"],
                                     "round5_fixed_sif")
                    self.assertIn("fixed", record["sif"]["msif_status"])
                    self.assertEqual(record["sif"]["known_wavelength_nm"], 759.0)
                    self.assertEqual(record["sif"]["known_Lnu"], schema.known_lnu)
                np.testing.assert_array_equal(plots.sif_curve(retrieved),
                                              plots.sif_curve(truth))
                self.assertEqual(plots.sif_known_anchor_radiance(retrieved),
                                 plots.sif_known_anchor_radiance(truth))
                if mode == "off":
                    np.testing.assert_array_equal(plots.sif_curve(truth), 0.0)

    def test_sif_panel_labels_both_parameters_fixed(self):
        for mode in ("off", "on"):
            with self.subTest(mode=mode):
                schema, retrieved, truth = self.records(mode)
                ensemble = {
                    "plot_schema": schema, "indices": [1, 2, 3],
                    "corrected": [retrieved] * 3,
                    "uncorrected": [retrieved] * 3,
                    "noiseless_paired": False,
                }
                figure, (axis, annotation_axis) = plots.plt.subplots(2, 1)
                try:
                    plots.draw_sif_card(
                        axis, annotation_axis,
                        {"ensemble": ensemble, "truth": truth, "state_index": 1},
                        (-0.1, 0.1), False, False,
                    )
                    self.assertIn("fixed", axis.get_title())
                    self.assertNotIn("retrieved slope", axis.get_title())
                    annotation = " ".join(text.get_text()
                                          for text in annotation_axis.texts)
                    if mode == "on":
                        self.assertIn("no SIF retrieval spread", annotation)
                    else:
                        self.assertIn("not retrieved", annotation)
                finally:
                    plots.plt.close(figure)

    def test_round5_routes_use_separate_output_and_prior_names(self):
        frontend = runpy.run_path(str(
            Path(__file__).resolve().parents[1] / "visualization" /
            "plot_physical_retrieval_ensembles.py"
        ))
        off = frontend["resolve_regime"](5, 0)
        on = frontend["resolve_regime"](5, 1)
        for regime, mode in ((off, "off"), (on, "on")):
            self.assertEqual(regime.state_model, "round5_fixed_sif")
            self.assertEqual(regime.co2_coordinate, "bottom_co2")
            self.assertIn("round5_fixed_sif_%s_tight_utls" % mode,
                          regime.source_prior_basename)
        self.assertEqual(off.sif_case, "off")
        self.assertEqual(on.sif_case, plots.SIF_CASE_ON)
        self.assertNotEqual(off.output_root, on.output_root)
        self.assertEqual(on.truth_table, frontend["CORRECTED_SIF_TRUTH_TABLE"])

    def test_round5_legend_does_not_claim_round4(self):
        handles, _ = plots.co2_figure_legend_groups(round4=True, round5=True)
        self.assertIn("round-5 fixed state", handles[0].get_label())
        self.assertNotIn("round-4", handles[0].get_label())
        handles, _ = plots.co2_figure_legend_groups(round4=True)
        self.assertIn("round-4 state-space truth", handles[0].get_label())


if __name__ == "__main__":
    unittest.main()
