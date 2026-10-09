# Analytical-Jacobian speedup integration audit

Audit date: 2026-09-09
Current study branch: `suniti_multi_sensor`
Inspected optimized branch: `origin/integration/release-candidate` at
`d441061ba29cec333ec5c57228c5e0f408c86906`

## Outcome

Christian's analytical-Jacobian optimizations are applicable to the
RRS-XCO2 retrievals and were explicitly benchmarked against this workflow.
They should be incorporated through the integrated release-candidate lineage,
not copied as isolated lines into the current branch.

The visible retrieval-adapter change is small:

```julia
params = copy_parameters(evaluator.base_parameters; share_luts=true)

result = rt_run_lin(
    model, lin_model; i_band=iband, sources,
    jacobian_basis=:local,
    jacobian_adding=:source)
```

However, the current `suniti_multi_sensor` source does not contain the APIs or
implementation files required by those calls, including
`copy_parameters`, local optical Jacobian bases, equivalent-source adding,
batched/blocked Jacobian products, and the reusable spectroscopy cache. The
adapter patch itself passes `git apply --check`, but applying it alone would
leave a runtime-incompatible program.

## What the optimized solver changes

The optimization keeps the public physical Jacobian layout and analytical
MOM equations, while changing how tangent work is represented and executed:

1. Products over wavelength and parameter axes are batched instead of creating
   per-column GPU views and launches.
2. The forward matrix inverse is reused for every derivative direction.
3. Lambertian source tangents are batched over wavelength and parameter.
4. A local optical basis avoids propagating unused physical directions.
5. Fixed aerosol microphysics is pruned from the local basis while retaining
   the requested AOD and height directions.
6. Equivalent-source adding propagates forcing vectors without constructing
   every full source-matrix derivative.
7. Trial parameter copies share read-only ABSCO/LUT storage rather than
   deep-copying several gigabytes per OE evaluation.
8. Parsed HITRAN spectroscopy is cached for workflows that do not use
   caller-owned LUTs.

The `OCO_RRS_synth` selection remains supported. Christian's workflow audit
records the same 12 band-local O2 columns and 22 band-local columns in each CO2
band, scattered into the common 30-column core state. Round 5 subsequently
removes both fixed SIF columns at its retrieval boundary.

## Reported RRS-XCO2 performance

On an A100 PCIe 40 GB, the archived-workflow comparison reports:

| Case | Historical evaluation total | Optimized retrieval | Ratio |
|---|---:|---:|---:|
| Clear state 001, uncorrected | 247.42 s | 20.69 s | 11.96x |
| Aerosol state 035, corrected | 698.15 s | 36.65 s | 19.05x |
| Aerosol state 035, uncorrected | 702.04 s | 35.80 s | 19.61x |

For the matched archived aerosol state, one complete forward-plus-Jacobian
evaluation fell from 117.289 s to a warmed median of 5.36763 s (21.85x).
These are shared-machine, cross-version workflow comparisons, not guaranteed
throughput multipliers.

The controlled 10,000-wavelength synthetic benchmark reports 38.50x for scalar
and 16.33x for IQU relative to the first recorded reference. A realistic
20-layer, 69-column, 10,000-point direct-HITRAN construction-plus-solve case
reports 3.174 s forward and 6.473 s forward plus Jacobians after spectroscopy
reuse. The ratio is scene- and timing-boundary-dependent; it must not be
advertised as universally below 2x.

## Independent checks performed here

The release candidate was exported into an isolated temporary checkout and
tested with Julia 1.12.5 on CPU. Compiled modules were disabled after the
shared `/tmp` filesystem caused a cache/linker failure. The following all
passed:

| Regression | Passed |
|---|---:|
| Batched Jacobian propagation | 416 |
| Batched Lambertian source Jacobians | 1,536 |
| Local optical basis | 5,164 |
| Independent parameter copies/read-only LUTs | 53 |
| Selective and `OCO_RRS_synth` Jacobians | 51 |
| Equivalent-source adding with SIF | 284 |
| **Total** | **7,504** |

The selective tests include two-sided finite differences for aerosol loading,
aerosol height, layer absorption, surface pressure, all three surface
coefficients, `SIF760`, and `mSIF`. The local-basis and source-adding tests also
compare against the retained analytical reference implementation. GPU results
were not independently rerun in this audit; the optimized branch retains its
documented CUDA evidence.

## Numerical caveat

The optimized lineage also contains intentional Float32 optical-assembly
repairs. Therefore an old archived retrieval and a new optimized retrieval are
not expected to be bit-identical. Christian's same-state replay found the new
matrix/source implementations mutually consistent within about `5e-5`
detector-noise units and a worst-column relative L2 difference of
`1.123e-5`. Relative to the archived calculation, the intentional precision
repair produced up to `0.0138` detector-noise units and about `0.146%` in the
most affected Jacobian column.

A broader six-pair retrieval comparison preserved convergence decisions and
kept the maximum XCO2 shift to `0.000449 ppm`, but three aerosol pairs reached
`0.02076` noise units and therefore failed a predeclared `0.01`-noise spectral
parity gate. This is small, diagnosed as precision-sensitive, and does not
invalidate the speedup; it does mean that campaign migration requires a fresh
round-5 closure test rather than assuming bitwise equivalence.

## Safe incorporation path

The two branches diverge before the surface-split integration. The optimized
release candidate already contains integrated/cherry-picked versions of the
multi-sensor and round-4 work, while the current branch lacks its required
solver foundation. Consequently:

1. preserve and validate the round-5 work as a small, reviewable patch;
2. base the optimized campaign on the current release-candidate lineage;
3. port round 5 into its relocated `sandbox/workflows/RRS_XCO2` tree;
4. apply the two-line retrieval adapter activation shown above;
5. run one clear and one aerosol, SIF-off and SIF-on noiseless closure case;
6. compare measurements in detector-noise units, every retained Jacobian
   column, OE decisions, final XCO2, and final state in prior-scaled units;
7. only then use the optimized path for production round-5 retrievals.

Directly cherry-picking only the late performance commits onto
`suniti_multi_sensor` is not recommended: the optimization series depends on
the integrated surface/source/cache types and spans many source files. A full
integration should be handled as a branch migration with the tests above, not
as a local hot-path edit.
