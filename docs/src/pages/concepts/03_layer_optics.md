# 3 · Layer Optical Properties — the bridge

> **For:** anyone reading the RT basics arc top-to-bottom. This is the page that connects "physics inputs" (gas absorption, Mie scattering) to "what the solver does." Centerpiece of the thread.
>
> **Prev:** [2 · Vector RTE & Discretization](02_rt_theory.md) · **Next:** [3a · Gas Absorption](03a_absorption.md)

[Concepts/02](02_rt_theory.md) ended on a four-array abstraction:
``(\tau, \varpi_0, \mathbf{Z}^{++}, \mathbf{Z}^{-+})``, one set per Fourier
moment per atmospheric layer. This page is the hinge: **everything before it
is "how to build those arrays" and everything after it is "what the solver
does with them."**

## Every atmospheric layer reduces to four arrays

For each Fourier moment ``m``, per layer:

| Symbol | Meaning | Shape |
|---|---|---|
| ``\tau`` | optical depth | `(nSpec,)` |
| ``\varpi_0`` | single-scattering albedo | `(nSpec,)` |
| ``\mathbf{Z}^{++}`` | forward-direction phase matrix | `(NquadN, NquadN, nSpec)` |
| ``\mathbf{Z}^{-+}`` | backscatter phase matrix | `(NquadN, NquadN, nSpec)` |

These are also the **three core RT variables** the linearized solver
differentiates against (``\mathbf{Z}^{++}`` and ``\mathbf{Z}^{-+}`` count as
one "core variable" because they share the same Greek expansion). Everything
upstream is plumbing to build them; everything downstream is linear algebra
over them. See [Concepts/06](06_linearization.md).

## How the four arrays are built (the operator chain)

The Julia type that carries them is `CoreScatteringOpticalProperties`, and
the constructor is `constructCoreOpticalProperties` in
[`src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:11–65`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl#L11-L65). The
build happens once per band, per Fourier moment, before the layer loop:

```@raw html
<pre class="mermaid">
flowchart LR
    G["Gas opacity τ_abs<br/>(HITRAN LBL)"]
    R["Rayleigh<br/>(greek + ϖ_Cab)"]
    A1["Aerosol_1<br/>(GreekCoefs, ϖ, k, fᵗ)"]
    A2["Aerosol_2<br/>(GreekCoefs, ϖ, k, fᵗ)"]
    C["<b>constructCoreOpticalProperties</b>"]
    L["Per-layer<br/>(τ, ϖ₀, Z⁺⁺, Z⁻⁺)"]
    M["MOM solver"]
    G --> C
    R --> C
    A1 --> C
    A2 --> C
    C --> L --> M
</pre>
<script type="module">
  if (!window.__mermaidLoader) {
    window.__mermaidLoader = import('https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs')
      .then(m => { window.__mermaid = m.default; m.default.initialize({ startOnLoad: false, theme: 'default' }); return m.default; });
  }
  window.__mermaidLoader.then(m => m.run({ querySelector: '.mermaid' }));
</script>
```

The build proceeds in five steps (the relevant lines from
`compEffectiveLayerProperties.jl:11–65`):

1. **Pick the Rayleigh greek source.** For pure-elastic runs (`noRS`), use the
   full Rayleigh greek — rotational Raman is rolled into the effective
   depolarization. For Raman-aware runs (`RRS`, `VS_*`), use the Cabannes
   greek and handle rotational Raman explicitly. Selector at lines 8–9. See
   [Concepts/08](08_inelastic.md).
2. **Build the Rayleigh layer's ``\mathbf{Z}^{++}, \mathbf{Z}^{-+}``** from
   the chosen greek and the quadrature points via
   `Scattering.compute_Z_moments(pol_type, μ, greek, m)`.
3. **For each aerosol**, run `compute_Z_moments` again (with that aerosol's
   `GreekCoefs`), wrap in a `CoreScatteringOpticalProperties` via
   `createAero` (which also applies δ-M truncation — see
   [Concepts/03c](03c_mixing.md)), and **mix it into the layer with `+`**.
4. **Add gas absorption** to the total optical depth:
   ``\tau_\lambda = \tau_\mathrm{abs} + \tau_\mathrm{scat}``,
   ``\varpi_\lambda = \tau_\mathrm{scat}\,\varpi_0 / \tau_\lambda``. Coded as
   `combo + CoreAbsorptionOpticalProperties(τ_abs)` (line 58).
5. **Stack layers vertically with `*`** (line 62):
   `prod([band_layer_props[i][iz] for i = 1:length(iBand)])`. The result is a
   `Vector{CoreScatteringOpticalProperties}`, one entry per atmospheric
   layer, ready for [Concepts/04](04_mom_solver.md).

The walkthrough in [Concepts/03c](03c_mixing.md) goes through the body of the
loop in detail, including the scattering-vs-total τ split that makes the
doubling step cheap.

## The two operators on `CoreScatteringOpticalProperties`

Mixing scatterers in a layer ([`src/CoreRT/types.jl:1063–1093`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/types.jl#L1063-L1093)):

```julia
function Base.:+(x::CoreScatteringOpticalProperties, y::CoreScatteringOpticalProperties)
    τ  = x.τ .+ y.τ                  # additive optical depth
    wx = x.τ .* x.ϖ                  # scattering-weight from each
    wy = y.τ .* y.ϖ
    w  = wx .+ wy
    ϖ  = w ./ τ                       # τϖ-weighted SSA
    Z⁺⁺ = (wx .* xZ⁺⁺ .+ wy .* yZ⁺⁺) ./ w   # τϖ-weighted phase matrix
    Z⁻⁺ = (wx .* xZ⁻⁺ .+ wy .* yZ⁻⁺) ./ w
    CoreScatteringOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺)
end
```

The weighting is ``w_x = \tau_x \cdot \varpi_x`` — proportional to the *scattering*
optical depth contribution from each species. That's what conserves the total
scattered intensity per stream pair when species are combined.

Stacking layers vertically (`src/CoreRT/types.jl:1096+`):

```julia
function Base.:*(x::CoreScatteringOpticalProperties, y::CoreScatteringOpticalProperties)
    arr_type = array_type(architecture(x.τ))
    x = expandOpticalProperties(x, arr_type)
    y = expandOpticalProperties(y, arr_type)
    CoreScatteringOpticalProperties(
        [x.τ; y.τ], [x.ϖ; y.ϖ],
        cat(x.Z⁺⁺, y.Z⁺⁺, dims=3),
        cat(x.Z⁻⁺, y.Z⁻⁺, dims=3))
end
```

This stacks the spectral axis (`cat(..., dims=3)`) for ``\mathbf{Z}`` and
concatenates the per-spectral ``\tau`` and ``\varpi``. The result is *not* a
column of layers (those are kept as a `Vector{CoreScatteringOpticalProperties}`)
— it's a multi-band concatenation along the spectral axis, used when several
bands share the same atmospheric column.

## A worked example

Consider an OCO-2-like O₂ A-band layer at altitude `iz`, with one Rayleigh
contribution and one fine-mode aerosol:

```julia
# inside constructCoreOpticalProperties, conceptually:

# (1) Rayleigh: τ_rayl from molecular cross-section, ϖ_Cab ≈ 1 for pure scattering
rayl = CoreScatteringOpticalProperties(τ_rayl, ϖ_Cab, RaylZ⁺⁺, RaylZ⁻⁺)

# (2) Aerosol fine mode, with δ-M applied via createAero
aer  = createAero(τ_aer_fine, aer_optics_fine, AerZ⁺⁺, AerZ⁻⁺)
#      → modifies (τ, ϖ) for forward-peak truncation per SS2015

# (3) Mix scatterers in the layer
combo = rayl + aer
#      → τ, ϖ, Z⁺⁺, Z⁻⁺ are τϖ-weighted averages

# (4) Add gas absorption (τ_abs is per-wavelength)
layer_optics = combo + CoreAbsorptionOpticalProperties(τ_abs)
#      → τ_λ = τ_abs + τ_scat
#        ϖ_λ = τ_scat·ϖ / τ_λ
#        Z⁺⁺, Z⁻⁺ unchanged
```

After all `Nz` layers are built, `rt_kernel!` consumes them one at a time
from TOA to BOA inside the Fourier-moment loop ([Concepts/04](04_mom_solver.md)).

## What's next

To understand how each ingredient is computed:

- [3a · Gas Absorption](03a_absorption.md) — how ``\tau_\mathrm{abs}`` is built (HITRAN line-by-line, Voigt line shapes).
- [3b · Mie & Rayleigh](03b_scattering.md) — how ``\mathbf{Z}^{++}``, ``\mathbf{Z}^{-+}``, ``\varpi``, ``k`` are built (Mie series, GreekCoefs, NAI-2 vs PCW).
- [3c · Mixing & δ-M Truncation](03c_mixing.md) — how the ingredients combine, and the scattering-vs-total τ split that makes the doubling tractable.

If you only care about what the solver does with these arrays, skip to
[Concepts/04](04_mom_solver.md).

## Code anchors

| Concept | Source |
|---|---|
| Build per-layer (τ, ϖ, Z) | [`src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:11–65`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl#L11-L65) |
| Aerosol δ-M wrapper (`createAero`) | `compEffectiveLayerProperties.jl:67–72` |
| Cabannes vs Rayleigh greek | `compEffectiveLayerProperties.jl:8–9` |
| Mix scatterers (operator `+`) | [`src/CoreRT/types.jl:1063–1093`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/types.jl#L1063-L1093) |
| Stack layers (operator `*`) | `src/CoreRT/types.jl:1096+` |
| `CoreScatteringOpticalProperties` type | `src/CoreRT/types.jl::CoreScatteringOpticalProperties` |
| Gas absorption type | `src/CoreRT/types.jl::CoreAbsorptionOpticalProperties` |

## Hands-on tutorials

Runnable examples with Plotly figures:

- [CoreRT walkthrough](../tutorials/Tutorial_CoreRT.md)
- [Configure a Scene](../tutorials/Tutorial_IO.md)

## References

- Sanghavi et al. (2014), JQSRT **133**:412–433, [doi:10.1016/j.jqsrt.2013.09.004](https://doi.org/10.1016/j.jqsrt.2013.09.004). §2 (layer-averaged optics, Eqs. C.22–C.24 for the forward-only ``\bar{\tau}, \bar{\varpi}_0, \bar{\mathbf{Z}}`` definitions used here).
- Sanghavi & Frankenberg (2023), JQSRT **311**:108791, [doi:10.1016/j.jqsrt.2023.108791](https://doi.org/10.1016/j.jqsrt.2023.108791). §3.2 — separating ``\tau_\mathrm{abs}`` from ``\tau_\mathrm{scat}`` for the constant-`N_doubl` design (more in [Concepts/03c](03c_mixing.md)).
- Crib sheet: `docs/dev_notes/theory_references.md` §F.
