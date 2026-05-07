# 2 · Vector RTE & Discretization

> **For:** readers building the equation-level mental model — atmospheric scientists, theorists, retrieval developers who want to understand what the solver actually solves.
>
> **Prev:** [1 · The Problem & MOM Thesis](01_overview.md) · **Next:** [3 · Layer Optical Properties](03_layer_optics.md)

This page writes the equation that vSmartMOM solves, then discretizes it.
By the end, the per-layer problem reduces to four arrays of shape
`(NquadN, NquadN, nSpec)` — the arrays the rest of the Concepts arc
manipulates.

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

## Quadrature in ``\mu``

After Fourier decomposition, each per-moment integral over ``\mu' \in [-1, 1]``
is replaced by a finite sum. vSmartMOM offers two quadrature schemes
([`src/CoreRT/tools/rt_set_streams.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/tools/rt_set_streams.jl)):

- **`GaussLegQuad`** — half-space Gauss-Legendre on ``[0,1]``; the solar
  zenith angle and viewing zenith angles are appended as zero-weight nodes.
  Cheap; interpolation between roots loses accuracy at off-quadrature angles.
- **`RadauQuad`** — block-Radau composite on ``[0, \mu_0] \cup [\mu_0, 1]``
  with ``\mu_0`` as a true quadrature node carrying non-zero weight (Sanghavi
  2014 App. B, Eqs. B.1–B.2). Required when the solar direction does not
  coincide with a Gauss root, because vSmartMOM excludes solar SFI from ``\mathbf{J}`` (see callout below) and therefore relies on direct ray-tracing
  through ``\mu_0``.

::: tip Design choice — no solar SFI in ``\mathbf{J}``
Sanghavi 2014 §2.2: vSmartMOM does *not* include single scattering of the
direct solar beam in the source term ``\mathbf{J}``. Two consequences:

1. The matrix-operator equations reduce from Sanghavi 2014 Eqs. (23)–(33) to
   Eqs. (23)–(28) for solar bands (no ``\mathbf{J}`` term). Doubling and
   adding are correspondingly cheaper.
2. For thermal-only (``m=0``) ``\mathbf{J}``, the source is isotropic, so its
   computation can be reused.

The price is needing Block-Radau quadrature for off-quadrature solar/viewing
angles. This is a deliberate departure from VLIDORT, SCIATRAN, and the
Plass/Fell/Hollstein matrix-operator formulations.
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

For each Fourier moment ``m``, every atmospheric layer is captured by four
arrays of shape `(NquadN, NquadN, nSpec)` where ``\text{NquadN} = N_\mathrm{quad} \cdot n_\mathrm{stokes}``:

| Symbol | Meaning | Shape |
|---|---|---|
| ``\tau`` | optical depth (per spectral point) | `(nSpec,)` |
| ``\varpi_0`` | single-scattering albedo | `(nSpec,)` |
| ``\mathbf{Z}^{++}`` | forward (downward → downward) Fourier-moment phase matrix | `(NquadN, NquadN, nSpec)` |
| ``\mathbf{Z}^{-+}`` | backscatter (downward → upward) Fourier-moment phase matrix | `(NquadN, NquadN, nSpec)` |

**The third dimension is spectral.** That is the whole batched-matmul story —
see [Concepts/07](07_architecture.md). Once these four arrays exist, the
elemental/doubling/adding kernels of [Concepts/04](04_mom_solver.md) are pure
linear algebra over the spectral axis.

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
