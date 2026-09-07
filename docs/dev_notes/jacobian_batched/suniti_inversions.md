# Suniti's RRS/XCO₂ inversion workflow and optimization transfer

Inspected on 2026-09-07: `/home/sanghavi/code/github/uni_vSmartMOM`, branch
`suniti_multi_sensor`, HEAD `cb90b83b`. The workflow has local changes and
untracked campaign code; this is not a claim about the committed tree alone.
The checkout, campaign products, and running retrievals were left untouched.
The optimization work is in `perf/jacobian-batched-propagation`, starting at
`075b3412` for this comparison. SIF support is committed as `00e670e1`; the
optical-precision and absorption-phase repairs are committed as `4284e13a`.

## What the inversions compute

The study measures the XCO₂ retrieval bias caused by rotational Raman scattering.
Both the corrected and uncorrected inversions use the **elastic analytic
Jacobian solver** (`noRS`). Corrected measurements come from Rayleigh truth;
uncorrected measurements come from Cabannes + RRS truth, fitted with the same
elastic forward model. Thus elastic Jacobian optimizations transfer to both
retrieval classes. Accelerating Raman truth generation is a separate problem.

The live implementation is in `RRS_XCO2/inversion/VSmartMOMForward.jl`,
`OptimalEstimation.jl`, `Round4KnownSIF.jl`, and the corresponding campaign
runners. Scene setup is in `scripts/common.jl` and
`config/oco_grass_3aerosol.yaml`.

- Float32 CUDA, IQU, nine weighted streams and a nadir viewing node: 30 × 30
  diffuse operators with external direct solar illumination.
- Sixteen atmospheric layers and three aerosol species: sulfate, organic
  carbon, and upper-atmosphere sulfate. Aerosol size distributions and indices
  are fixed; loading and vertical location are retrieved.
- Three complete spectral batches: 2,735 O₂, 1,281 weak-CO₂, and 995 strong-CO₂
  points, including instrument shoulders. Instrument projection produces
  2,742 measurements through Stokes analysis, spectral-density conversion,
  Gaussian ILS convolution, and resampling.
- Twelve retrieved CO₂ layers, surface pressure, three log AODs, three log
  aerosol heights, nine Legendre surface coefficients, and two SIF parameters:
  30 active coordinates. The upper four CO₂ layers are fixed.
- `OCO_RRS_synth` already selects band-local Jacobian columns: 12 for O₂ and
  22 for each CO₂ band. It already skips Mie microphysics and H₂O derivatives.
- O₂ carries `SolarBeam + SurfaceSIF`; the CO₂ bands carry `SolarBeam` alone.
  The explicit physical solar spectrum retains Fraunhofer structure. Fourier
  convergence is enabled; inspected logs finish at `m_used=3` rather than
  evaluating all orders through the configured ceiling of 15.

The optimal-estimation loop uses a Gaussian prior and diagonal measurement-noise
covariance with damped Gauss–Newton/Levenberg–Marquardt steps. Each trial,
including a rejected trial, computes a complete forward model and Jacobian.
The accepted-state evaluation is already retained when retrying a rejected
step. A converged retrieval performs an additional terminal evaluation.

The bottom-layer study spans four surfaces, two aerosol conditions, two SIF
conditions, and five bottom-layer CO₂ values: 80 scenes. Each has corrected and
uncorrected retrievals with ten noise realizations plus one noiseless case,
giving 1,760 retrievals in a complete campaign. Known-SIF Round 4 wraps the
30-coordinate evaluator and then removes fixed SIF coordinates/columns; it
does not yet remove their work inside the RT/instrument calculations.

## Where campaign time goes

The audit sampled 64 completed NetCDF retrievals from
`retrievals_acos_mapped_tapered_vertical_correlation_nosif`: 32 aerosol and
32 clear cases, covering both correction classes and noise/noiseless results.
These are pooled campaign observations across workers, not a controlled GPU A/B.

| Median recorded quantity | Aerosol | Clear |
|---|---:|---:|
| Trial forward + Jacobian evaluation | 108.17 s | 62.71 s |
| RT/Jacobian solve | 105.16 s | 61.10 s |
| Model construction timer | 1.30 s | 0.74 s |
| Instrument calculation | 0.052 s | 0.030 s |
| OE matrix algebra | 0.00049 s | 0.00040 s |
| Sum of evaluation times per retrieval, including converged terminal evaluation | 691.87 s | 298.86 s |

RT accounts for about **96.5%** of aggregate evaluation time. The 177 sampled
aerosol trials include 34 rejections; the 112 clear trials include none.
Instrument convolution and the small OE matrix solve are secondary targets.

`evaluate_oco_forward` starts with `deepcopy(evaluator.base_parameters)` **outside
the model-construction timer**. A direct isolated measurement copied
3,767,189,504 host bytes in 2.98 s, including the loaded absorption interpolators.
The evaluator loads its LUTs once, but deep-copies their coefficient arrays on
every trial. Sharing read-only spectroscopy while copying mutable atmospheric,
aerosol-profile and surface state is a concrete remaining opportunity.

The loaded tables are legacy ABSCO O₂/CO₂/H₂O LUTs plus an O₂-band HITRAN H₂O
LUT. Our parsed-HITRAN cache does not remove this deepcopy. Likewise the native
Mie derivative speedup does not explain gains here: this study already disables
those derivatives. Forward Mie optics are still reconstructed each trial even
though their microphysical inputs are fixed.

## Replay methodology

An isolated Julia driver loads the live workflow adapter against our package
and replays the archived final state in corrected
`retrieval_state035_perturbation10.nc`. It preserves her state mapping, LUTs,
solar spectrum, grids, geometry, Fourier convergence, and instrument operator.
It never includes a campaign runner. All generated products go to a separate
scratch directory. Timings exclude one complete warm-up per mode and use three
complete evaluations, with CUDA synchronization and garbage collection between
samples. Hardware: curry GPU 0, A100 PCIe 40 GB; Julia 1.12.6, four Julia threads,
one OpenBLAS thread.

Before enabling source adding for SIF, the optimized matrix path took a median
**17.040 s** end to end, including **14.009 s** RT. Using source adding on the two
solar-only CO₂ bands reduced this to **14.675 s**, including **11.669 s** RT.
The archived terminal evaluation for this same state recorded **117.289 s**.
That historical comparison suggests about an eightfold reduction, but combines
code/environment and run-time differences; it is not a freshly controlled
benchmark of the old branch.

The optimized matrix replay differed from the archived measurement by at most
0.000890 detector standard deviations. The largest relative L2 difference in
any Jacobian column was 3.33e-5. These are single-state numerical comparisons,
not validation of a complete migrated retrieval campaign or its convergence.

Enabling source adding for O₂ initially exposed a separate optical-assembly
issue. Local matrix adding reproduced the same discrepancy from physical
columns (about 0.13% in the worst raw O₂ Jacobian column), while source adding
agreed with local matrix adding within 1.24e-5 per-column relative L2. All three
stopped at Fourier order 3. The SIF boundary implementation was not responsible.
Stronger optical-level tests found a literal `1.0` in the physical Rayleigh
albedo, promoting Float32 mixtures to Float64. The repair uses `one(FT)` and
retains the same successive forward mixing order in the local path. Its compact
derivative basis is unchanged. This intentionally corrects the precision of the
old physical path; archived Float32 campaign results should therefore be
compared with a scientific tolerance, not bitwise equality.
Pure absorption now also preserves the incoming scattering-phase derivatives
and appends exact zero gas-phase columns, eliminating a cancelling quotient in
the physical path. The forward optical arrays are checked for exact equality
between physical and local coordinates. Phase-derivative comparisons retain
relative tolerances and allow 10 Float32 epsilons per successive scatterer
addition near cancelled entries; the measured three-aerosol differences there
are approximately 1.2–2.2e-6 for derivatives around 1e-3. Exact gas zeros and
the stricter Float64/finite-difference checks are retained.

The combined final regression runs pass **4,100 CPU checks** (3,724 optical
basis, 284 SIF, 92 solar-only source adding) and **1,397 CUDA checks** (1,021
optical basis, 284 SIF, 92 solar-only source adding). A strict
Documenter/Vitepress build also passes. Metal was not hardware-tested.

## Final replay after the precision repairs

Same archived state, all 5,011 wavelengths and the complete instrument operator;
medians of three warmed evaluations on the A100:

| Current solver configuration | Complete evaluation | RT/Jacobian | Model construction | Instrument |
|---|---:|---:|---:|---:|
| Matrix adding, automatic basis selection | 16.038 s | 12.523 s | 0.451 s | 0.0485 s |
| Local basis, source adding in all three bands including SIF | 14.766 s | 11.333 s | 0.479 s | 0.0473 s |

Component medians need not sum to the total median. The total includes parameter
deepcopy and other adapter work outside the component timers. Host allocations
are 4.750 GB and 4.723 GB per evaluation, respectively. The all-band source path
is **1.086× faster** than the current matrix path for this selected layout.
SIF support is now available, but its addition alone does not materially improve
on the earlier mixed O₂-matrix/CO₂-source configuration; pruning the unused
microphysics directions remains the next main RT opportunity.

Current source/matrix measurement differences are at most **5.0e-5 detector
standard deviations** (relative L2 2.00e-8). The largest relative L2 difference
in any instrument-level Jacobian column is **1.123e-5**, or 0.00113%.
The replay checks measurement differences below 0.001 noise units and every
Jacobian column below 3e-4 relative L2.

Relative to the archived state, both implementations after the precision repair differ
by at most 0.0138 detector standard deviations and about 0.146% in the most
affected Jacobian column. That shared shift follows the intentional precision
repair, not a remaining source/matrix disagreement. The historical terminal
timing, 117.289 s, is **7.94×** the new 14.766 s evaluation. This supports a large
retrieval speedup, while retaining the historical-environment caveat above.
Complete inversion convergence and campaign products have not been revalidated.

Raw timings, per-column comparisons, source hashes and reproduction commands are
in [the evidence directory](evidence/suniti_replay/README.md).

## Integration priorities

1. **Completed in `00e670e1`:** enable surface SIF in equivalent-source adding,
   with 284 SIF checks plus 92 solar-only checks passing on each of CPU and CUDA.
   The boundary
   derivative must distinguish reflected solar attenuation from locally emitted
   fluorescence. See [the source-adding derivation](source_adding.md).
2. Prune fixed microphysics directions from the local basis. With three aerosols
   the current basis has 17 directions, although this study needs only the two
   scalar optical directions and three phase-mixture directions. O₂ has only
   seven selected atmospheric columns already, so forcing 17 local directions
   is not automatically faster. Pruning could reduce its doubling basis to five.
3. Share immutable loaded LUTs across trial-state copies. Cache fixed forward
   Mie/truncated Greek/phase-node data across iterations with explicit keys for
   all microphysical, spectral and truncation settings.
4. Cache repeated initial-state evaluations where the entire forward state is
   identical. The inspected 64 files contain only four distinct initial vectors,
   each repeated 16 times. Include fixed upper CO₂ and all configuration/source
   dependencies in the cache key; the state vector alone is insufficient.
5. Remove known-SIF columns at the early selection boundary in Round 4. Consider
   forward-only rejected-trial screening only with a measured cost/reuse model;
   computing forward and then recomputing it for accepted Jacobians can erase
   the saving.

Port through an isolated adapter/worktree and compare complete representative
retrievals before changing an active campaign. Preserve the priors, covariance,
convergence tests, instrument shoulders and physical source spectrum: changing
these would confound the RRS-bias experiment.
