# 2 · Vector RTE & Discretization

> **For:** readers building the equation-level mental model — atmospheric scientists, theorists, retrieval developers who want to understand what the solver actually solves.
>
> **Prev:** [1 · The Problem & MOM Thesis](01_overview.md) · **Next:** [3 · Layer Optical Properties](03_layer_optics.md)

This page writes the equation that vSmartMOM solves, then discretizes it.
By the end, the diffuse per-layer problem reduces to square operators of shape
`(NquadN, NquadN, nSpec)`. An opt-in source-function path may evaluate the
phase matrix on one additional, external solar direction without adding that
direction to the diffuse operators.

## The vector RTE

For a plane-parallel atmosphere illuminated by direct solar flux ``\mathbf{S}_0`` along
``(\mu_0, \phi_0)``, the Stokes vector ``\mathbf{L}(\tau, \mu, \phi)`` for diffuse light obeys
(Sanghavi 2014, Eq. 2):

```math
\mu \frac{d\mathbf{L}(\tau,\mu,\phi;\mu_0,\phi_0)}{d\tau} =
- \mathbf{L}(\tau,\mu,\phi;\mu_0,\phi_0)
+ (1-\varpi_0)\mathbf{B}(T)
+ \frac{\varpi_0}{4\pi}\mathbf{Z}(\mu,\phi;\mu_0,\phi_0)\mathbf{S}_0\,e^{-\tau/\mu_0}
+ \frac{\varpi_0}{4\pi}\int_0^{2\pi}\!\!\int_{-1}^{1}\mathbf{Z}(\mu,\phi;\mu',\phi')\,\mathbf{L}(\tau,\mu',\phi';\mu_0,\phi_0)\,d\mu'\,d\phi'.
```

Reading the right-hand side term by term:
1. ``-\mathbf{L}`` — extinction along the path.
2. ``(1-\varpi_0)\mathbf{B}(T)`` — thermal emission (zero in solar bands).
3. The single-scatter direct-beam source, exponentially attenuated from TOA.
4. The diffuse multiple-scatter integral over all incoming directions ``(\mu',\phi')``.

**Sign convention:** ``\mu > 0`` is downward, ``\mu < 0`` is upward.
``\varpi_0 \in [0,1]`` is the single-scattering albedo.
``\mathbf{Z}`` is the 4×4 phase matrix (Sanghavi 2014, Eq. 8). The diffuse
solution ``\mathbf{L}`` adds to the direct delta-source ``\mathbf{S}_0\delta(\mu-\mu_0)\delta(\phi-\phi_0)e^{-\tau/\mu_0}``
to give the total radiation field; we solve only for the diffuse part to avoid
the singularity.

### Boundary conditions

At the top of atmosphere (Sanghavi 2014, Eq. 6):

```math
\mathbf{L}(0,\mu,\phi;\mu_0,\phi_0) = \mathbf{S}_0\,\delta(\mu-\mu_0)\,\delta(\phi-\phi_0)
\quad\text{for}\;\mu \geq 0.
```

At the surface, an emissivity ``\boldsymbol{\varepsilon}_g`` and a polarized
BRDF ``\mathbf{R}_p(\mu,\phi;\mu',\phi')`` close the system (Sanghavi 2014,
Eqs. 33–37). [Concepts/05](05_surfaces.md) covers the surface in detail.

## Fourier decomposition in azimuth

Because the problem is azimuthally symmetric about the solar direction, ``\mathbf{Z}``
can be expanded as a finite Fourier series in ``(\phi-\phi_0)`` (Sanghavi 2014, Eq. 8):

```math
\mathbf{Z}(\mu,\phi;\mu',\phi') = \tfrac12\,\mathbf{C}_0(\mu,\mu')
+ \sum_{m=1}^{M}\!\bigl[\bar{\mathbf{C}}_m(\mu,\mu')\cos m(\phi-\phi_0) + \bar{\mathbf{S}}_m(\mu,\mu')\sin m(\phi-\phi_0)\bigr].
```

``\bar{\mathbf{C}}_m`` is block-diagonal (mixing intensity ``I`` with linear
polarization ``Q``); ``\bar{\mathbf{S}}_m`` is block-off-diagonal (mixing ``U``
with the rest). Each Fourier moment ``m`` solves an *independent* RTE — the
``m``-loop in `rt_run.jl` is therefore embarrassingly parallel.

## Quadrature in μ

After Fourier decomposition, each per-moment integral over ``\mu' \in [-1, 1]``
is replaced by a finite sum. vSmartMOM offers two quadrature schemes
([`src/CoreRT/tools/rt_set_streams.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/tools/rt_set_streams.jl)):

- **`GaussLegQuad`** — half-space Gauss-Legendre on ``[0,1]``. By default
  (`external_solar=false`), the solar direction is appended to the diffuse
  operator as a zero-weight node alongside the viewing directions. With the
  opt-in `external_solar=true` (TOA-only), only viewing directions are
  appended; exact ``μ_0`` remains a scalar source direction and a
  phase-evaluation column.
- **`RadauQuad`** — block-Radau composite on ``[0, \mu_0] \cup [\mu_0, 1]``
  with ``\mu_0`` as a true quadrature node carrying non-zero weight (Sanghavi
  2014 App. B, Eqs. B.1–B.2). This is the historical DNI-oriented solution to
  off-quadrature solar geometry. It retains embedded ``μ_0`` and does not support
  external-solar mode.

::: tip Historical DNI design and the current external-solar SFI path
Sanghavi 2014 §2.2 deliberately excluded direct-solar single scattering from
``\mathbf{J}``. That historical formulation uses the shorter matrix-operator
recurrences (Eqs. 23–28) and Block-Radau (App. B) to place ``μ_0`` on the
quadrature. This explains the DNI/Radau lineage, but it is not a statement
about today's SFI kernel.

The current solver computes direct-beam source-function integration (SFI).
Its backward-compatible default keeps ``μ_0`` in the operator grid as a
zero-weight node because several solver paths still depend on that layout.
The opt-in external-solar form (`external_solar=true`, a TOA-only path
consumed through `rt_run_toa`) removes this representational relic: ``μ_0`` is
**not a diffuse stream**. The elemental solver receives it as a scalar
propagation direction. Exact rectangular ``Z_0(μ_i,μ_0)`` columns generate
finite-layer ``R_0/T_0`` operators, which are then contracted into
``\mathbf J``. Diffuse R/T operators remain square on the diffuse grid.

This path is restricted to Gauss quadrature, elastic forward/linearized
`noRS`, forward rotational Raman `RRS`, Lambertian surfaces, and endpoint TOA
output. Call `rt_run_toa(model)` or `rt_run_toa(rs,model)`; BOA, HDR, and BHR
arrays are not allocated. Raman Jacobians, VRS, single-scatter-only,
non-Lambertian, and interior-sensor calculations retain embedded ``μ_0``.
:::

## The supermatrix form (per Fourier moment)

For each Fourier moment ``m``, the RTE becomes (Sanghavi 2014, Eq. 12):

```math
\mu \frac{d\mathbf{I}_m(\tau,\mu,\mu_0)}{d\bar{\tau}}
= -\mathbf{I}_m(\tau,\mu,\mu_0)
+ \frac{\varpi_0}{2}\!\int_{-1}^{1}\!\bar{\mathbf{Z}}_m(\mu,\mu')\mathbf{I}_m(\tau,\mu',\mu_0)\,d\mu'
+ (1-\varpi_0)\mathbf{B}(T)\,\delta_{0m}.
```

After quadrature in ``\mu``, the upwelling and downwelling streams are stacked
into supervectors ``\mathbf{I}_m^+, \mathbf{I}_m^-`` of length ``4N_\mathrm{quad}``
(Sanghavi 2014, Eqs. 14–15):

```math
\frac{d}{d\bar{\tau}}\!\begin{bmatrix}\mathbf{I}_m^+\\\mathbf{I}_m^-\end{bmatrix}
= (1-\varpi_0)\mathbf{B}(T)\delta_{0m}\!\begin{bmatrix}\mathcal{M}^{-1}\\-\mathcal{M}^{-1}\end{bmatrix}
+ \frac{\varpi_0}{2}\,\mathcal{M}^{-1}\!\begin{bmatrix}\bar{\mathbf{Z}}_m^{(++)}\mathcal{C} - \mathbb{E} & \bar{\mathbf{Z}}_m^{(+-)}\mathcal{C}\\-\bar{\mathbf{Z}}_m^{(-+)}\mathcal{C} & -\bar{\mathbf{Z}}_m^{(--)}\mathcal{C}+\mathbb{E}\end{bmatrix}\!\begin{bmatrix}\mathbf{I}_m^+\\\mathbf{I}_m^-\end{bmatrix}.
```

``\mathcal{M} = \operatorname{diag}(\mu_i)`` and ``\mathcal{C} = \operatorname{diag}(c_i)``
are diagonal supermatrices of quadrature angles and weights. The four blocks
``\bar{\mathbf{Z}}_m^{(\pm\pm)}`` are the directional pairings of the phase matrix.

The ``+``/``-`` superscripts on ``\mathbf{I}, \mathbf{R}, \mathbf{T}, \mathbf{J}``
that you see throughout the source code (`r⁻⁺`, `t⁺⁺`, `R⁻⁺`, etc.) are exactly
this convention: ``+`` is downward, ``-`` is upward.

## The four arrays each layer reduces to

For each Fourier moment ``m``, every atmospheric layer is captured by two
spectral vectors and two phase-matrix arrays. Define
``\text{NquadN} = N_\mathrm{quad} \cdot n_\mathrm{stokes}`` for the diffuse
operator and ``\text{NphaseN} = N_\mathrm{phase} \cdot n_\mathrm{stokes}``
for phase evaluation:

| Symbol | Meaning | Shape |
|---|---|---|
| ``\tau`` | optical depth (per spectral point) | `(nSpec,)` |
| ``\varpi_0`` | single-scattering albedo | `(nSpec,)` |
| ``\mathbf{Z}^{++}`` | forward (downward → downward) Fourier-moment phase matrix | `(NphaseN, NphaseN, nSpec)` |
| ``\mathbf{Z}^{-+}`` | backscatter (downward → upward) Fourier-moment phase matrix | `(NphaseN, NphaseN, nSpec)` |

**The third dimension is spectral.** That is the whole batched-matmul story —
see [Concepts/07](07_architecture.md). Normally `NphaseN == NquadN`. With
external-solar SFI, phase evaluation may have one additional direction; the
diffuse elemental/doubling/adding operators remain
`(NquadN, NquadN, nSpec)`, and only the exact solar phase column is consumed
outside that block. Five weighted streams plus one distinct VZA in
`Stokes_IQU` therefore produces an ``18×18`` diffuse operator, not ``21×21``.

## Polarization is a type, not a runtime branch

`Stokes_I` (``n=1``, scalar intensity), `Stokes_IQ` (``n=2``), `Stokes_IQU` (``n=3``),
`Stokes_IQUV` (``n=4``) each subtype `AbstractPolarizationType` in
[`src/Scattering/types.jl:92–143`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/types.jl#L92-L143). Each carries the polarization symmetry vector
``\mathbf{D}`` (used in [Concepts/04](04_mom_solver.md#doubling)) and the unit
solar Stokes vector ``\mathbf{I}_0``.

The kernels specialize on the type at compile time. There is no
`if pol_type.n == 4` branch inside the inner loops — Julia's multiple dispatch
emits separate code paths for each polarization at JIT time. Adding `Stokes_IQ`
took adding one struct definition; the kernels picked it up automatically.

## Code anchors

| Concept | Source |
|---|---|
| Top-level RT loop | [`src/CoreRT/rt_run.jl:53–329`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/rt_run.jl#L53-L329) |
| Fourier-moment dispatch | [`src/CoreRT/rt_run.jl:208`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/rt_run.jl#L208) (`for m = 0:max_m-1`) |
| Quadrature construction | [`src/CoreRT/tools/rt_set_streams.jl:24–110`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/tools/rt_set_streams.jl#L24-L110) |
| `QuadPoints` struct | `src/CoreRT/types.jl::QuadPoints` |
| External direct-solar phase coupling | `src/CoreRT/CoreKernel/elemental.jl::get_elem_rt_SFI!` |
| Lean TOA-only entry point | `src/CoreRT/rt_run.jl::rt_run_toa` |
| Polarization types | [`src/Scattering/types.jl:92–143`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Scattering/types.jl#L92-L143) |
| Phase matrix Z from Greek | `src/Scattering/compute_Z_matrices.jl::compute_Z_moments` |
| Supermatrix layer types | `src/CoreRT/types.jl::CompositeLayer/AddedLayer` |

## Hands-on tutorials

Runnable examples with Plotly figures:

- [CoreRT walkthrough](../tutorials/Tutorial_CoreRT.md)
- [Quick Start](../tutorials/Tutorial_QuickStart.md)

## References

- Sanghavi et al. (2014), *vSmartMOM*, JQSRT **133**:412–433, [doi:10.1016/j.jqsrt.2013.09.004](https://doi.org/10.1016/j.jqsrt.2013.09.004). Eqs. (2), (4), (6), (8)–(10), (12), (14)–(15); App. B (Block-Radau). **Primary reference for everything on this page.**
- Chandrasekhar (1950), *Radiative Transfer*. (Original vector RTE in this form.)
- Hovenier (1971), *Symmetry relationships for scattering of polarized light…*, J. Atmos. Sci. **28**:488. (Block diagonal/anti-diagonal structure of ``\bar{\mathbf{Z}}_m``.)
- Full crib sheet: `docs/dev_notes/theory_references.md` §A, §E.
