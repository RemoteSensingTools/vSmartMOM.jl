# Phase 6 script port list — triage (2026-04-22)

Source: [vSmartMOM.jl/test/benchmarks/](/home/sanghavi/code/github/vSmartMOM.jl/test/benchmarks/) (sanghavi @ 9ee9a75, READ-ONLY worktree).

Target: [test/benchmarks/](../test/benchmarks/) in sanghavi-unified.

Phase 6 goal: bring sanghavi's research-workflow scripts into unified where they exercise RT paths we want to keep regression-testable. Not all scripts are worth porting — many are one-off research notebooks with deep data-file / external-package dependencies (Insolation, ClimaParams, per-campaign NetCDF blobs). The rule of thumb is:

- **port** if it exercises a supported RT path (noRS / RRS / VS) with a public config and adds scenario coverage beyond what's already in `test/benchmarks/`;
- **drop** if it's a plotting script, research one-off, or depends on non-public data;
- **defer** (port later, gated) if large config data is needed and the scenario is valuable but not in the critical path for Phase 6.

**User sign-off required** before any porting starts (per the pickup-plan heuristic: research-code counts drift from 8 → 23 quickly). The budget should reflect the *port* column below, not the full 48-script count.

## Already ported (ignore)

| Script | Lines | Status |
|---|---|---|
| 6SV1_R_trues.jl | 109 | Already in [test/benchmarks/](../test/benchmarks/). |
| benchmark_6SV1.jl | 87 | Already in [test/benchmarks/](../test/benchmarks/). |
| benchmark_natraj.jl | 85 | Natraj.yaml + natraj_trues.jl already in unified. |
| natraj_trues.jl | 52 | Already in unified. |
| prototype_EMIT_aer_ht.jl | 863 | Ported Phase 3c (`6f54002`). |
| creategrid_O2Aband_RamanSIF.jl | 927 | Ported Phase 3d (`6b01912`). |
| raman_batched_ops_benchmark.jl | 365 | Already exists as `batched_*_benchmark.jl` in unified. |
| raman_optimization_baseline.jl / _compare.jl | 68 / 182 | Superseded by harness (scenarios.toml + run_benchmarks.jl). |

## Port (exercise supported RT paths, small-medium)

These are the minimum-viable additions — each demonstrates a scenario not already covered by the harness and has self-contained inputs.

| Script | Lines | Why | Risk |
|---|---|---|---|
| prototype_inelastic.jl | 92 | Smallest RRS demo against [O2Parameters.yaml](../test/test_parameters/O2Parameters.yaml). Useful as a tutorial-style forward-only example. | Low. |
| prototype_O2ABand_RRS.jl | 586 | Full O₂-A band RRS workflow. Gives wider spectral + physics coverage than Phase1b's toy. | Low-medium — check Plots/JLD2 deps. |
| prototype_O2ABand_RRS_SIF.jl | 808 | RRS + Lambertian SIF; parallel to RamanSIF acceptance scripts but without the grid sweep. | Low-medium. |
| prototype_VS_O2Aband.jl | 442 | Vibrational Scattering (VS) smoke — currently no VS script in unified test/benchmarks/. Unlocks VS regression coverage. | Medium — VS paths less battle-tested after Phase 1b. |
| prototype_O2BBand_RRS_SIF.jl | 833 | O₂-B band mirror of O2A + SIF. Same shape as already-ported scripts. | Low (once O2A_SIF is up). |

**Sub-total: 5 scripts, ~2800 lines.**

## Consider (larger grids / broader scenarios; port if Phase 6 has budget)

| Script | Lines | Why | Risk |
|---|---|---|---|
| creategrid_O2Bband_RamanSIF.jl | 907 | O₂-B band grid sweep — mirror of already-ported O2A. | Low. |
| creategrid_betwnAB_RamanSIF.jl | 908 | Between-band sweep for SIF retrievals. | Low. |
| creategrid_CO2Wband_OCORayl.jl | 943 | OCO CO₂-W band Rayleigh grid. Stretches coverage to CO₂ spectra. | Medium — needs OCO LUT data. |
| creategrid_O2Aband_OCORayl.jl | 927 | OCO O₂-A Rayleigh grid (noRS). | Medium — same data dependency. |
| compare_rt_EMIT.jl | 350 | Comparison driver for EMIT scenarios. Not yet in unified. | Low — check NetCDF input paths. |
| emit_modtran_noRS_scenarios.jl | 803 | EMIT+MODTRAN noRS acceptance script. | Medium. |
| plot_ramanXS.jl | 205 | Plots Raman cross-sections; good for docs/tutorial. Stripped of Plots deps it becomes a data-dump utility. | Low. |
| sif_raman.jl | 163 | Small SIF + Raman demo. | Low. |

**Sub-total: 8 scripts, ~5200 lines.**

## Drop (research notebooks / external deps / superseded)

These are intentionally not ported. Pull individually into Phase 7 docs or a `contrib/` directory on user request.

| Script | Lines | Reason |
|---|---|---|
| O3Huggins_polRaman.jl | 928 | O₃ Huggins polarized-Raman research. Specific to one campaign. |
| RT_vs_hydrolight.jl | 833 | Comparison against Hydrolight (external radiative-transfer code). |
| RamanSIFPolyCoeffMaps.jl | 1022 | Insolation / ClimaParams deps + campaign data. |
| RamanSIF_oco_retr.jl | 1648 | OCO retrieval driver. Deep data deps. |
| RamanSIFmaps.jl | 275 | Plotting script with Insolation. |
| RamanSIFspectra.jl | 1924 | Research one-off, Insolation deps. |
| RamanSIFworldmaps.jl | 1641 | Research one-off, Insolation deps. |
| TestSZA.jl | 134 | Specific SZA sanity check, narrow scope. |
| compLelli.jl / compLelli2.jl | 250 / 495 | Comparison vs Lelli et al. Campaign-specific. |
| plot_MODvsMOM.jl | 33 | MOD vs MOM plotting wrapper. |
| prototype_CarbonI.jl | 825 | CarbonI mission research. |
| prototype_fraunhofer_VS_spectrum.jl | 174 | Fraunhofer VS demo; overlap with `prototype_VS_O2Aband.jl`. |
| prototype_inelastic_OCO2.jl | 378 | OCO-2 research prototype. |
| prototype_inelastic_VS2.jl | 754 | VS variant prototype, overlap with `prototype_VS_O2Aband.jl`. |
| prototype_inelastic_ms.jl | 474 | Multi-sensor prototype; requires `rt_run_multisensor.jl`. |
| prototype_inelastic_sunspot_VS.jl | 777 | Sunspot VS research. |
| testAer_O2Aband_Raman.jl | 465 | Research eval against specific aerosol parameterizations. |
| testAerosol.jl | 482 | Ditto. |
| testCPU.jl | 891 | Research CPU dump with mixed scenarios. |
| testRayCabRaman.jl | 106 | Tiny Rayleigh/Cabannes/Raman demo; superseded by Phase 1b test. |
| evalAer_O2ABbands_Raman.jl | 1084 | Aerosol evaluation research. |
| evalAer_O2Aband_Raman.jl | 927 | Ditto. |
| test_creategrid_O2Aband_RamanSIF.jl | 928 | Test wrapper for already-ported grid script. |
| test_creategrid_noabs_RamanSIF.jl | 958 | No-absorption variant of the same. |
| create_HITRAN_LUTs.jl | 158 | One-time LUT generation. Move to `scripts/` if wanted, not benchmarks. |

**Sub-total: 25 scripts, ~16 000 lines — explicitly NOT in Phase 6 scope.**

## Recommended Phase 6 shape

- **Scope**: the 5 "Port" entries above (~2800 lines, 5 acceptance scripts).
- **Budget**: 2-3 sessions (~1 week). Each script needs an input-YAML review, a dependency strip (Plots / Revise / `CUDA.device!(1)`), and a smoke-test assertion added to `test/runtests.jl`.
- **Exit criteria**:
  1. Each ported script runs end-to-end against current `sanghavi-unified` tip.
  2. A `prototype_*_smoke.jl` test exists that runs the first ~30 lines of setup through `rt_run` and checks finite output (similar to `test_forward_noRS.jl`).
  3. Output JLD2 under `test/reference/` for any script whose result is non-trivial to re-derive.
- **Deferred to Phase 6b / 7** (user call):
  - The 8 "Consider" entries if user wants broader scenario coverage.
  - Any of the 25 "Drop" entries explicitly requested.

User sign-off needed on this shape before the port work begins. If the shape changes ("also port creategrid_O2Bband_RamanSIF and compare_rt_EMIT"), update the Port / Consider split here first.
