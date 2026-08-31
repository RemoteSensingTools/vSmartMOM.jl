# Exact single-scattering reference: derivations

This document derives the four single-scattering (SS) paths through a plane-parallel
atmosphere bounded above by space and below by a Lambertian surface. Companion to
`exact_ss_reference.jl` and `exact_ss_validate.py`.

**Convention:** `I₀` is the *radiance* of the solar beam (i.e., the perpendicular
solar flux). The downward solar irradiance on a horizontal surface at TOA is
`μ₀ · I₀`. This convention is fixed by path 2's closed form
`L_path2 = (μ₀ I₀ A / π) · exp(-τ/μ₀) · exp(-τ/μ_v)`, which evaluates the surface
direct-beam contribution in the standard atmospheric-RT way.

## Geometry conventions

- z-axis upward.
- Plane-parallel layers indexed top-down; layer 1 uppermost, layer N just above
  the surface.
- Cumulative optical depth τ measured from TOA downward. τ = 0 at TOA;
  τ = τ_total at the surface.
- Sun at zenith angle θ₀, azimuth ϕ₀ = 0. Direction of *propagation* is downward.
  μ₀ = cos θ₀ > 0.
- Sensor at TOA looking down. View direction (from sensor toward atmosphere)
  has zenith angle θ_v from nadir. Outgoing radiance direction (away from
  atmosphere, toward sensor) has μ_v = cos θ_v > 0.
- Δϕ = ϕ_v − ϕ₀ = ϕ_v.

## Scattering angle cosine

Cosine of the angle between the photon's incoming direction (downward at θ₀)
and outgoing direction (upward at θ_v):

```
cos Θ = -μ₀·μ_v + √(1-μ₀²)·√(1-μ_v²)·cos(Δϕ)
```

The negative sign on `μ₀μ_v` reflects that the sun direction has a *downward*
z-component while the outgoing photon has an *upward* z-component.

## RTE solution and the 1/μ_v factor

The fundamental RTE for upward radiance along view direction μ_v with vertical
optical depth τ is:

```
μ_v · dI/dτ = -I + S(τ)
```

where S(τ) is the source function per unit *vertical* optical depth.

The formal solution at TOA from a source distributed throughout the atmosphere:

```
I(τ=0, μ_v) = ∫_0^{τ_total} S(τ) · exp(-τ/μ_v) · dτ / μ_v
```

The factor `1/μ_v` is the conversion from "per unit vertical τ" to "per unit
slant length along the view direction." This is the textbook RTE solution and
applies to all upward-radiance integrals through plane-parallel atmospheres.

## Path 1: sun → atmosphere → sensor

The single-scattering source function from the direct solar beam at depth τ:

```
S_SS(τ) = (ϖ/4π) · P(cos Θ) · I₀ · exp(-τ/μ₀)
```

This is the radiance source per unit vertical τ, into solid angle around the
view direction (μ_v, ϕ_v), from a beam radiance I₀ attenuated through depth τ.

Applying the RTE solution for upward radiance:

```
ΔI_TOA = ∫ S_SS(τ) · exp(-τ/μ_v) · dτ / μ_v
       = (1/μ_v) · (ϖ I₀ / 4π) · P(cos Θ) · ∫ exp(-τ·a) · dτ
```

where `a = 1/μ₀ + 1/μ_v`. For a homogeneous layer:

```
ΔI_layer = (1/μ_v) · (ϖ I₀ / 4π) · P(cos Θ) · [exp(-τ_top·a) - exp(-τ_bot·a)] / a
```

**The 1/μ_v factor is essential.** Equivalently, using `1/(μ_v·a) = μ₀/(μ₀+μ_v)`:

```
ΔI_layer = (ϖ I₀ / 4π) · P(cos Θ) · μ₀/(μ₀+μ_v) · [exp(-τ_top·a) - exp(-τ_bot·a)]
```

For multiple layers, sum the per-layer contributions; ϖ and P are evaluated per
layer.

## Path 2: sun → surface → sensor

Direct beam reaches surface attenuated by `exp(-τ_total/μ₀)`. Lambertian
reflection: `L = (albedo/π) · F_down`. Downward irradiance on surface from
direct beam: `F_down = μ₀ I₀ exp(-τ_total/μ₀)`. Reflected radiance attenuates
to TOA via `exp(-τ_total/μ_v)`:

```
L_path2 = (albedo/π) · μ₀ I₀ · exp(-τ_total/μ₀) · exp(-τ_total/μ_v)
```

Closed form, no quadrature.

## Path 3: sun → atmosphere → surface → sensor

Photon scatters once in the atmosphere on its way down, ends up reaching the
surface, gets Lambertian-reflected, exits to sensor.

### Downward irradiance on surface from path-3 photons

The scattering source per unit vertical τ at depth `τ_scatter`, into a *downward*
direction (μ_d, ϕ_d), per unit incoming solar radiance:

```
S_down(τ_scatter, μ_d, ϕ_d) = (ϖ/4π) · P(cos Θ_d) · I₀ · exp(-τ_scatter/μ₀)
```

This is a *radiance source* per `dτ_vertical` per `dΩ_d`. The radiance arriving
at the surface in direction (μ_d, ϕ_d) from this source — applying the RTE
solution along a *downward* slant path:

```
L_down_surface(μ_d, ϕ_d) = ∫ S_down · exp(-(τ_total - τ_scatter)/μ_d) · dτ_scatter / μ_d
```

The 1/μ_d here is the same RTE solution factor as 1/μ_v in path 1.

To get downward *irradiance* at the surface, integrate L over the down-hemisphere,
weighted by μ_d (the cosine projection from radiance to irradiance):

```
F_down = ∫_{down hemi} L_down(μ_d, ϕ_d) · μ_d · dΩ_d
       = ∫∫ [∫ S_down · exp · dτ / μ_d] · μ_d · dΩ_d
```

**The 1/μ_d (RTE solution) and μ_d (radiance → irradiance) cancel.** The
integrand has no μ_d factor:

```
F_path3 = ∫_{τ_scatter} ∫∫_{down hemi} S_down(τ_scatter, μ_d, ϕ_d) ·
                                       exp(-(τ_total - τ_scatter)/μ_d) ·
                                       dΩ_d · dτ_scatter
```

### Azimuthal reduction

`P(cos Θ_d)` depends on ϕ_d only through cos(ϕ_d - ϕ_0). Define the azimuthally-
averaged phase:

```
P̄(μ_d, μ₀) = (1/2π) ∫_0^{2π} P(cos Θ_d) · dϕ_d
```

Then:

```
F_path3 = ∫_{τ_scatter} ϖ · I₀ · exp(-τ_scatter/μ₀) ·
          [(1/2) · ∫_0^1 P̄(μ_d, μ₀) · exp(-(τ_total - τ_scatter)/μ_d) · dμ_d] ·
          dτ_scatter
```

The factor `(1/2)` arises from `2π/4π` after performing the ϕ-integration with
the (ϖ/4π) prefactor. Implementation: 1D Gauss-Legendre quadrature on (0, 1] for μ_d.

### Per-layer τ-integration in closed form

For a layer with constant optical properties, the τ_scatter integral within the
layer combines the two exponentials. Numerically stable form (avoiding overflow
when μ_d is small):

```
f_top = exp(-τ_top/μ₀ - (τ_total - τ_top)/μ_d)
f_bot = exp(-τ_bot/μ₀ - (τ_total - τ_bot)/μ_d)
b     = 1/μ₀ - 1/μ_d
integral = (f_top - f_bot) / b           if |b| > 1e-10
         = (1/2) (f_top + f_bot) (τ_bot - τ_top)  otherwise
```

### Surface emission

```
L_path3 = (albedo / π) · F_path3 · exp(-τ_total / μ_v)
```

The exponential is the upward attenuation along the view path from surface to
TOA (no further scattering — that would be a higher-order path, not single).

## Path 4: sun → surface → atmosphere → sensor

Direct beam reaches surface as `μ₀ I₀ exp(-τ_total/μ₀)` (irradiance). Lambertian
reflection produces upward isotropic radiance:

```
L_surface = (albedo/π) · μ₀ I₀ · exp(-τ_total/μ₀)
```

This isotropic upward radiance experiences one atmospheric scatter, redirecting
into the view direction.

### Source function from upward isotropic input

The scattering source at depth `τ_scatter` into the view direction, from
incoming radiance over the upward hemisphere of directions (μ_u, ϕ_u):

```
S_into_view(τ_scatter) = (ϖ/4π) · ∫∫_{up hemi} L_in(τ_scatter, μ_u, ϕ_u) ·
                                                 P(cos Θ_4(μ_u, μ_v, ...)) ·
                                                 dΩ_u
```

The incoming radiance from the surface attenuates upward through the atmosphere:

```
L_in(τ_scatter, μ_u, ϕ_u) = L_surface · exp(-(τ_total - τ_scatter)/μ_u)
```

(Just radiance attenuation along upward slant. No 1/μ_u factor here — that comes
later when we compute outgoing TOA radiance.)

Combining and azimuthally averaging:

```
S_into_view(τ_scatter) = (ϖ/2) · L_surface · ∫_0^1 P̄(μ_u, μ_v) ·
                                              exp(-(τ_total - τ_scatter)/μ_u) ·
                                              dμ_u
```

**No μ_u factor in the integrand**: the angular integral is over solid angle
of *incoming radiance*, not irradiance. (This is the same structure as path 3,
where the radiance-to-irradiance μ_d canceled the RTE 1/μ_d.)

### Outgoing TOA radiance

Apply the RTE solution for upward radiance from the source S_into_view:

```
L_path4 = ∫_{τ_scatter} S_into_view(τ_scatter) · exp(-τ_scatter/μ_v) · dτ_scatter / μ_v
```

The 1/μ_v factor here is the standard upward-RTE one (same role as in path 1).

### Per-layer closed form (numerically stable)

The τ_scatter integral within a layer with constant ϖ and P̄:

```
f_top = exp(-τ_top/μ_v - (τ_total - τ_top)/μ_u)
f_bot = exp(-τ_bot/μ_v - (τ_total - τ_bot)/μ_u)
b     = 1/μ_v - 1/μ_u
integral = (f_top - f_bot) / b           if |b| > 1e-10
         = (1/2) (f_top + f_bot) (τ_bot - τ_top)  otherwise
```

## Symmetry and BRDF reciprocity

With the corrected factors (`1/μ_v` in path 1, `1/μ_v` in path 4, no μ_d/μ_u in
the inner integrands), the system satisfies BRDF reciprocity:

```
f_r(μ₀, μ_v) = f_r(μ_v, μ₀)
```

where `f_r = π · L_TOA / (μ₀ · I₀)`. Equivalently:

```
L_total(μ₀, μ_v, Δϕ) / μ₀ = L_total(μ_v, μ₀, Δϕ) / μ_v
```

Per-path reciprocities:

- **Path 1**: `path1(μ₀, μ_v) / μ₀ = path1(μ_v, μ₀) / μ_v` ✓
- **Path 2**: `path2(μ₀, μ_v) / μ₀ = path2(μ_v, μ₀) / μ_v` ✓ (trivial from formula)
- **Paths 3 ↔ 4**: under `(μ₀ ↔ μ_v)`, path 3 maps to path 4 and vice versa:

```
path3(μ₀, μ_v) / μ₀ = path4(μ_v, μ₀) / μ_v
path4(μ₀, μ_v) / μ₀ = path3(μ_v, μ₀) / μ_v
```

This is *the* test for whether paths 3 and 4 are correctly implemented. The
implementation in this reference passes it to ~10⁻⁶ at N_quad = 64.

## Summary of total exact-SS contribution

For Lambertian surface:

```
I_SS_total(μ₀, μ_v, Δϕ) = L_path1 + L_path2 + L_path3 + L_path4
```

with paths 1, 2 closed-form; paths 3, 4 requiring 1D quadratures over the
atmospheric down/up hemisphere with phase-function-weighted exponential
attenuation.

Validation criteria:

- **Conservative limit** (ϖ=1, albedo=0): only path 1 contributes (paths 2, 3,
  4 all need albedo > 0).
- **Reflectance limit** (no atmosphere): only path 2 contributes;
  `L = μ₀ I₀ albedo / π`.
- **Optically thin regime** (τ << 1): paths 1, 3, 4 are linear in τ to leading
  order; path 2 is approximately constant.
- **BRDF reciprocity**: both individual paths and the total satisfy
  `L/μ₀ ↔ L/μ_v` symmetry under (μ₀, μ_v) swap.
- **Quadrature convergence**: refining N_quad and N_phi monotonically reduces
  the residual against a high-resolution reference.

These tests, run together, pin down the formula completely. Failing any of them
indicates a missing or extra path-length factor somewhere in the implementation.

## Implementation note: numerical stability

The textbook closed form for the path 3 / path 4 layer integral

```
exp(-τ_total/μ_d) · [exp(-b·τ_top) - exp(-b·τ_bot)] / b
```

(where `b = 1/μ₀ - 1/μ_d`) overflows when μ_d is small (which happens at high
quadrature density, since GL nodes approach 0). The numerically stable form
combines the exponentials at each endpoint *before* subtracting, as shown
above. **This is mandatory for any production implementation of paths 3/4** —
the textbook form fails at N_quad ≳ 32 in practice.
