"use strict";

// Paths are relative to RRS_XCO2/bottom_layer_XCO2_retrievals.
const ROUNDS = [3, 4, 5, 6];
const CO2_LEVELS = [360, 380, 400, 420, 440];
const PRODUCTS = [
  { key: "profiles_surface", label: "Aerosol + surface" },
  { key: "co2_profiles", label: "CO₂ + SIF" },
];

function plotRoot(round, sif) {
  if (!ROUNDS.includes(round)) throw new RangeError("Unknown retrieval round");
  if (round === 3) {
    return "retrievals_acos_mapped_tapered_vertical_correlation_" +
      (sif ? "sif" : "nosif") + "/plots/physical_ensembles";
  }
  const campaign = round === 4 ? "round4_known_sif759" : "round" + round + "_fixed_sif";
  return campaign + "/plots_" + (sif ? "sif" : "nosif") + "/physical_ensembles";
}

function buildGroups() {
  const groups = [];
  // Truth-table order: [001,021,041,061], ..., [020,040,060,080].
  for (const aerosol of ["no_aerosol", "with_aerosol"]) {
    for (const sif of [false, true]) {
      for (const ppm of CO2_LEVELS) {
        const number = groups.length + 1;
        const states = [0, 20, 40, 60].map((offset) =>
          String(number + offset).padStart(3, "0"));
        const label = ppm + " ppm · " +
          (aerosol === "no_aerosol" ? "no aerosol" : "with aerosol") +
          " · " + (sif ? "with SIF" : "no SIF");
        const slug = sif ? "sif_angular_integral760_0p5" : "nosif";
        const panes = PRODUCTS.flatMap((product) => ROUNDS.map((round) => ({
          round,
          product: product.key,
          label: "Round " + round + " — " + product.label,
          path: plotRoot(round, sif) + "/bottom_co2_" + ppm + "_" + slug +
            "_" + aerosol + "_" + product.key + ".png",
        })));
        groups.push({ number, states, ppm, aerosol, sif, label, panes });
      }
    }
  }
  return groups;
}

function wrapIndex(index, count) {
  if (!Number.isInteger(index) || !Number.isInteger(count) || count <= 0) return 0;
  return ((index % count) + count) % count;
}

if (typeof module !== "undefined") module.exports = { buildGroups, plotRoot, wrapIndex };
if (typeof window !== "undefined") window.rrsModel = { wrapIndex };
