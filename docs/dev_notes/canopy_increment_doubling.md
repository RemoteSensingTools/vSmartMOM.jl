# Canopy transmission increments

The directional canopy elemental operator approaches the identity as its
thickness tends to zero. Rounding `T = I + O(δ)` in Float32 and then repeatedly
squaring can amplify lost diagonal increments into apparent absorption or gain.
Promoting only the doubling stage cannot recover information lost in the seed.
Coarser seeds avoid that cancellation but increase finite single-scattering
initialization error. The canopy path therefore carries the increment directly.

## Algebra

All operators below use the existing D-transformed basis during doubling.
Let `Δ = T-I`, `S = (I-RR)⁻¹`, `H = S-I = S RR` and
`C = T S-I = Δ+H+ΔH`. For identical slabs, standard adding-doubling gives

```
Q = R + RΔ
R_new = R + Q + C Q
Δ_new = C + Δ + C Δ
```

No commutativity is assumed. `H` is obtained by solving `(I-RR)H=RR`.
For the elemental diagonal, form `expm1(-x) + exp(-x)*scattering`
directly, where `x = δτ G/μ`; never subtract one from an already rounded T.
Off-diagonal transmission is already an increment.

The solar attenuation increment is `β=expm1(-δτ G₀/μ₀)`. Its doubling
update is `β_new=β(2+β)`. If `j_d` and `j_u` are the current slab sources,
the sources of the lower copy are `j_d1=j_d+β j_d`, `j_u1=j_u+β j_u`.
Using old R and old sources throughout each right-hand side:

```
v_u = R j_d + j_u1
j_u_new = j_u + v_u + C v_u
v_d = R j_u1 + j_d
j_d_new = j_d1 + v_d + C v_d
```

Only after all doublings is `T=I+Δ` reconstructed. Existing D-matrix
restoration then supplies the reverse operators and upward solar source.
There is no phase, flux, absorption or operator normalization correction.

## Scope and tests

`test_canopy_increment.jl` compares this algebra with ordinary doubling at
matched Float64 seed thickness and compares Float32 with Float64 at thin
seeds. It spans Fourier orders 0–2, Stokes I/IQU/IQUV, diffuse and combined
diffuse/specular optics, one and three spectral points, and multiple LAIs.
The optional CUDA run requires `VSMART_CANOPY_INCREMENT_GPU=true` and disables
scalar indexing. `test_canopy_precision_surface.jl` checks the actual canopy
wrapper against a conservative photon budget, including multiple physical
layers over a perfectly reflecting soil. Settings and tolerances are in TOML.

Local validation on Julia 1.12.6 and an NVIDIA L40S: 3,696 operator checks on
CPU and 3,456 on CUDA passed, including zero-thickness layers. The production
surface tests passed 48 CPU and 8 CUDA checks; the existing canopy suite passed
135 checks, including its soil-albedo derivative regression. For the dense
conservative fixture, direct-budget error at 64 streams was approximately
1.2e-6 in Float32 and 3.2e-6 in Float64. The 16-stream angular discretization
retained approximately 1.8e-4 error in both precisions; the strict 1e-4 budget
target is checked at 32 and 64 streams. No result is rescaled to hit that target.

The combined specular Float32 cases used the separate CanopyOptics fix in
[PR #64](https://github.com/RemoteSensingTools/CanopyOptics.jl/pull/64).
These local checks do not replace the complete vSmartMOM suite or upstream CI.

These tests establish numerical agreement over their stated cases, not
physical correctness of every specular parameterization or infinite-resolution
conservation. Angular resolution, leaf-angle quadrature and elemental thickness
must be converged separately. Metal and arbitrary dual-number differentiation
have not been validated for this path. The existing canopy soil-albedo
linearization uses finite differences of the forward wrapper; this change does
not add leaf or canopy analytic Jacobians.

This implementation carries the solar source only and explicitly rejects
additional source slots. A distributed fluorescence source cannot use the
solar attenuation recurrence: its two sublayers emit according to their own
source distribution. It needs a separate source law and elemental integral.
The reflection/transmission algebra is reusable when those sources are added.

The current `batch_solve!` backend and temporary allocation strategy prioritize
validation. Throughput and allocation measurements are required before making
performance claims or selecting this as an ESM production backend.
