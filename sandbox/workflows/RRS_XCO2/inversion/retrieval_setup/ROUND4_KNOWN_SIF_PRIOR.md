# Round-4 known-SIF prior

Round 4 treats the SIF spectral radiance at 759 nm as exactly known. The
retrieval source remains referenced internally at 760 nm, so the one free SIF
coordinate in a SIF-on retrieval is `mSIF = dLnu/dnu` and

```text
SIF760 = Lnu759 + mSIF * (nu760 - nu759).
```

For corrected-v2 SIF-on truth, the full truth template gives

```text
Lnu759     = 0.004818031987713776 mW m-2 sr-1 (cm-1)-1
Llambda759 = 0.0836346275560863   mW m-2 sr-1 nm-1
```

The slope remains constrained exactly as in round 3 in its physical
wavelength-space definition at 760 nm:

```text
dLlambda/dlambda = 0 +/- 0.002625 mW m-2 sr-1 nm-2 (one sigma).
```

With the 759-nm amplitude fixed, this becomes

```text
mSIF mean  = -7.34275712494799e-7
mSIF sigma =  8.780708772523119e-6
```

in native `mW m-2 sr-1 (cm-1)-2` units. The corresponding *derived* 760-nm
coordinate has mean `0.004830761266413151` and sigma
`0.00015222087186261526` in `mW m-2 sr-1 (cm-1)-1`.

For SIF-off scenes, `SIF760=0` and `mSIF=0` are retained in the canonical
34-element record with exactly zero variance. Neither is included in the
numerical solve.

The active full-state mappings are therefore:

```text
SIF-on:  [1; 6:32; 34]  (29 parameters)
SIF-off: [1; 6:32]      (28 parameters)
```

## Production rule

Use only `build_round4_known_sif_apriori.jl` to construct round-4 campaign
priors. It copies `xa[1:32]` and `Sa[1:32,1:32]` exactly from the approved
tapered round-3 prior and records that file's SHA-256.

Do not instead Gaussian-condition the old two-coordinate SIF covariance. The
old covariance includes an uncertain SIF amplitude. Conditioning that joint
distribution on a known 759-nm amplitude changes the mean and variance of the
physical slope, whereas the round-4 requirement is to retain the existing
wavelength-slope constraint.

The generator emits separate SIF-on and SIF-off NetCDF files plus matching
human-readable `.dat` audits under
`bottom_layer_XCO2_retrievals/round4_known_sif759/retrieval_setup/` by
default. It refuses to replace an existing file unless
`ROUND4_PRIOR_OVERWRITE=1` is explicitly set.
