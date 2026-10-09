# OCO-2 XCO2 retrieval convergence criteria

## Key distinction

An OCO-2 full-physics retrieval does **not** converge merely because every
spectral residual lies within the instrumental noise. The algorithm treats
the following as separate questions:

1. **Has the nonlinear optimal-estimation iteration converged?** This is
   decided from the size of the state-vector update relative to the posterior
   state uncertainty.
2. **Is the converged spectral fit consistent with the assumed measurement
   uncertainties?** This is assessed afterward from a band-by-band normalized
   chi-square statistic.

Consequently, a retrieval can converge mathematically while retaining a poor
spectral fit. This distinction is explicit in the OCO-2 outcome codes.

## Optimal-estimation problem

For state vector `x`, measurement vector `y`, forward model `F(x)`, prior
state `x_a`, prior covariance `S_a`, and measurement-error covariance `S_e`,
the retrieval minimizes

```text
J(x) = (y - F(x))' S_e^-1 (y - F(x))
     + (x - x_a)' S_a^-1 (x - x_a).
```

OCO-2 uses a Levenberg-Marquardt modification of the Gauss-Newton update. The
Levenberg-Marquardt parameter starts at 10 in the documented standard
configuration and is adjusted according to the ratio between the actual and
linear-model-predicted reduction in the cost function.

## State-space convergence test

After a successful, non-divergent iteration, the solver evaluates an error
variance derivative of the form

```text
d_sigma_sq = Delta_x' * S_hat^-1 * Delta_x,
```

where `S_hat` is the posterior state covariance and `Delta_x` is the proposed
state update used for the convergence test. Thus `d_sigma_sq` is the squared
state update measured in posterior-standard-deviation units.

The OCO-2 Level-2 ATBD writes the stopping rule as

```text
d_sigma_sq < f * N_active,
```

where `N_active` is the number of active state-vector elements and `f` is a
factor of order unity. The public NASA implementation explicitly constructs

```text
d_sigma_sq_scaled = d_sigma_sq / N_active
```

and declares convergence when this scaled value is below the configured
threshold. The public standard OCO configuration sets

```text
d_sigma_sq_scaled < 2.0
maximum iterations       = 7
maximum divergent steps  = 2
initial LM gamma          = 10.0.
```

A trial step is classified as divergent when the ratio of actual to predicted
cost reduction is less than or equal to `1e-4`. Divergent steps are rejected
and the damping is increased.

These criteria measure whether further iterations would materially change the
retrieved state. They do not impose a channel-by-channel radiance-residual
limit.

## Spectral-fit quality test

Once state-space convergence has been reached, OCO-2 separately calculates a
normalized or reduced chi-square for each spectral band:

```text
chi2_red[band] = residual[band]' * S_e[band]^-1 * residual[band]
                 / N_band.
```

The public standard configuration uses

```text
chi2_red[band] < 1.4
```

for every fitted band to classify the converged fit as provisionally
acceptable. With a diagonal covariance and correctly specified independent
noise, this is equivalent to a noise-normalized weighted RMS residual below

```text
sqrt(1.4) = 1.183 noise standard deviations.
```

It does **not** require every channel to lie within `+/-1 sigma`. For Gaussian
noise, about 31.7% of individual residuals are expected to lie outside
`+/-1 sigma` even for a statistically correct model. A channel-by-channel
`1 sigma` stopping rule would therefore reject almost every sufficiently long
spectrum.

The Level-2 ATBD also notes that spectroscopy and instrumental systematics
make residuals inconsistent with detector-noise estimates, especially at high
signal-to-noise ratio. The `max_chi2` value is therefore empirical and is a
fit-quality classification rather than the optimizer's stopping condition.

## Retrieval outcomes

The documented OCO-2 outcome values expose the separation directly:

| Outcome | Meaning |
|---:|---|
| 1 | State-space convergence reached and all bands pass the spectral-fit quality test. |
| 2 | State-space convergence reached, but at least one band has a poor fit relative to the assumed uncertainties. |
| 3 | Convergence was not reached within the maximum number of iterations. |
| 4 | The maximum number of divergent steps was exceeded. |

Outcome 2 demonstrates that convergence does not imply a noise-consistent
spectral fit.

## Implication for the RRS-XCO2 retrieval experiment

The corrected and uncorrected synthetic retrievals should likewise keep
iteration convergence and spectral-fit diagnostics separate:

1. Stop the nonlinear solve using a posterior-normalized state-step test.
2. Record the total cost, number of iterations, number of rejected/divergent
   steps, and `d_sigma_sq_scaled`.
3. Record `chi2_red` and noise-normalized RMS separately for the O2 A, weak
   CO2, and strong CO2 bands.
4. Assign separate convergence and fit-quality flags rather than merging them
   into one boolean.

Whether an RRS-minus-no-RS spectral difference lies inside a plotted noise
cloud is not by itself a convergence test. The `noRS` forward model can adjust
CO2, aerosols, surface properties, SIF, or surface pressure to absorb a
spectral discrepancy. It can therefore converge to a biased state even when
the final residual is comparable to instrumental noise. Conversely, a
noise-free corrected synthetic measurement can legitimately produce
`chi2_red << 1` because no random noise realization was added.

For an OCO-like implementation, the documented public defaults provide a
reasonable initial specification:

```text
converged:              d_sigma_sq_scaled < 2.0
fit provisionally OK:  chi2_red[band] < 1.4 for all bands
maximum iterations:    7
maximum divergences:   2
divergent-step ratio:   actual/predicted cost reduction <= 1e-4.
```

These should remain configurable rather than hard-coded, and retrieval output
should preserve the underlying diagnostics so that alternative thresholds can
be evaluated without rerunning the forward model.

## Documentation scope

NASA currently identifies OCO-2 v11.2r/v11.3r as the reference data record,
but the latest publicly posted Level-2 Full Physics ATBD is Version 3.0,
Revision 1 (January 2021), describing the v10 algorithm. NASA continues to
point users to this ATBD and the public RT Retrieval Framework. No public
v11 document located for this review specifies a replacement for the
state-step convergence logic described above. The exact numerical thresholds
used by any production version should nevertheless be recorded with that
version's processing configuration.

## Primary sources

- [OCO-2 and OCO-3 Level 2 Full Physics Retrieval ATBD, Version 3.0 Rev. 1](https://docserver.gesdisc.eosdis.nasa.gov/public/project/OCO/OCO_L2_ATBD.pdf): Section 3.5 describes the inverse method, state-step convergence rule, outcome values, and bandwise goodness-of-fit test.
- [NASA RT Retrieval Framework: `connor_convergence.cc`](https://github.com/nasa/RtRetrievalFramework/blob/master/lib/Implementation/connor_convergence.cc): implements divergent-step handling, the `d_sigma_sq_scaled` stopping test, iteration limits, and the separate bandwise chi-square quality classification.
- [NASA RT Retrieval Framework: `connor_solver.cc`](https://github.com/nasa/RtRetrievalFramework/blob/master/lib/Implementation/connor_solver.cc): defines the optimal-estimation update and computes `d_sigma_sq`, `d_sigma_sq_scaled`, and measurement chi-square.
- [NASA standard OCO configuration: `oco_base_config.lua`](https://github.com/nasa/RtRetrievalFramework/blob/master/input/oco/config/oco_base_config.lua): supplies `threshold=2.0`, `max_iteration=7`, `max_divergence=2`, `max_chisq=1.4`, and `gamma_initial=10.0`.
- [OCO-2 Data Center](https://ocov2.jpl.nasa.gov/science/oco-2-data-center/): identifies the current reference products and links NASA's public Level-2 Full Physics code and documentation.
- [NASA GES DISC OCO-2/3 document index](https://docserver.gesdisc.eosdis.nasa.gov/public/project/OCO/OCO-2-3_Documents.pdf): identifies the latest public Level-2 Full Physics ATBD and current product documentation.
