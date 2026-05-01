# 8 · Inelastic Extension (brief)

> **For:** users who need to model the Ring effect or hyperspectral retrievals where rotational/vibrational Raman scattering matters. UV/Vis instruments, O₂ A-band, deep Fraunhofer lines.
>
> **Prev:** [7 · Architecture-Agnostic Code](07_architecture.md) · **Next:** [Manual → Compute Jacobians](../jacobians.md) (loop back to task layer)

This page is intentionally brief — vSmartMOM's Raman support is a parallel
solver path next to the elastic one, not central to most retrievals. For
the full RT formalism see Sanghavi & Frankenberg 2023 (Part II); for the
optical-property derivation see Sanghavi 2022 (Part I).

## In one paragraph

vSmartMOM supports rotational Raman (RRS) and vibrational Raman (VRS,
RVRS) inelastic scattering as a **parallel kernel path** next to the
elastic one. Selection is via the `RS_type` configured in YAML: `noRS`
(default; pure elastic), `RRS`, `VS_0to1`, `VS_1to0`, plus `_plus`
variants for multi-band runs that share spectral structure. The bulk of
retrievals (NIR shortwave for greenhouse gases, polarized aerosol over
ocean) use `noRS` — the inelastic correction is small there. UV/Vis
retrievals (Ring effect, O₂ A-band line shapes, ozone DOAS) need RRS.

## When to enable it

| Application | Mode |
|---|---|
| Ring effect (hyperspectral UV/Vis trace gases) | `RRS` |
| O₂ A-band line-shape correction at high precision | `RRS` (single-scatter approximation often suffices via `rt_run_ss`) |
| Fraunhofer-line ghosting due to N₂/O₂ vibrational Raman | `VS_0to1` and/or `VS_1to0` |
| NIR-only retrievals (XCO₂, XCH₄ from OCO-2/3, GOSAT) | `noRS` (default) |

## How it's wired

The Raman path uses the **same matrix-operator language** as the elastic
path. Sanghavi & Frankenberg 2023 §3.3 shows that the elastic adding
equations (`R₂₀, T₂₀, J₂₀`) carry over essentially unchanged; the
inelastic adding equations (Eqs. 16–21 of SF2023-II) build the inelastic
contributions on top under a **linear-in-inelastic-scattering
approximation** — only one inelastic event per photon path, multiple
elastic scattering before and after. Second-order inelastic has been shown
to be negligible at typical atmospheric conditions.

Architecturally:

- **Mode types** subtype `AbstractRamanType` and live in
  [`src/Inelastic/types.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Inelastic/types.jl). Each carries the precomputed Cabannes/Raman
  greek coefficients and Dunham coefficients for N₂/O₂.
- **Parallel kernel files** at [`src/CoreRT/CoreKernel/elemental_inelastic.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/elemental_inelastic.jl),
  `doubling_inelastic.jl`, `interaction_inelastic.jl`. `rt_kernel!`
  dispatches on `RS_type` to choose which set runs (see
  [`src/CoreRT/CoreKernel/rt_kernel.jl:48–229`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/rt_kernel.jl#L48-L229)). The elastic kernel
  *also* runs alongside — Raman is additive on top.
- **`rt_run_ss`** ([`src/CoreRT/rt_run.jl:364–524`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/rt_run.jl#L364-L524)) is the single-scattering
  approximation that exists alongside `rt_run` for SF2023-II Eq. (32) —
  the inelastic correction `I₁ ≈ I₀ + (I'₁ − I'₀)` that makes the
  linear-in-inelastic approximation usable in absorbing bands. Not a
  debug helper; it's the *production* tool for fast O₂ A-band Raman
  corrections.
- **Linearized Raman is currently elastic-only**. The `_lin.jl` files
  cover the elastic path only; ``\partial \mathbf{R}/\partial \mathbf{x}`` derivatives
  through inelastic scattering are not yet implemented.

## Code anchors

| Concept | Source |
|---|---|
| Mode types (`noRS`, `RRS`, `VS_0to1`, `VS_1to0`, `_plus` variants) | [`src/Inelastic/types.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Inelastic/types.jl) |
| Inelastic optical properties (Cabannes, RRS, VRS, RVRS) | `src/Inelastic/` |
| Inelastic elemental | [`src/CoreRT/CoreKernel/elemental_inelastic.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/elemental_inelastic.jl) |
| Inelastic doubling | [`src/CoreRT/CoreKernel/doubling_inelastic.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/doubling_inelastic.jl) |
| Inelastic interaction | [`src/CoreRT/CoreKernel/interaction_inelastic.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/interaction_inelastic.jl) |
| `_plus` variants (multi-band) | [`src/CoreRT/CoreKernel/elemental_inelastic_plus.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/CoreKernel/elemental_inelastic_plus.jl) |
| Single-scatter inelastic correction | `src/CoreRT/rt_run.jl::rt_run_ss:364–524` |
| Cabannes vs full Rayleigh selector | [`src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:8–9`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl#L8-L9) |

For developers extending the Raman path (adding a new mode, adding new
molecular species), see [Add a Raman Mode](../extending/raman.md).

## References

- **Sanghavi (2022)**, *Raman scattering in the earth's atmosphere, part I: Optical properties*, JQSRT **291**:108328, [doi:10.1016/j.jqsrt.2022.108328](https://doi.org/10.1016/j.jqsrt.2022.108328). **Optical-property derivation; Cabannes vs Rayleigh; RRS/VRS/RVRS phase matrices and cross-sections.**
- **Sanghavi & Frankenberg (2023)**, *Raman scattering in the Earth's atmosphere, Part II: Radiative transfer modeling for remote sensing applications*, JQSRT **311**:108791, [doi:10.1016/j.jqsrt.2023.108791](https://doi.org/10.1016/j.jqsrt.2023.108791). **Inelastic RT formalism in the matrix-operator method; constant-`N_doubl` trick (Eqs. 8–9); single-scatter correction (Eq. 32).**
- Crib sheet: `docs/dev_notes/theory_references.md` §I.
