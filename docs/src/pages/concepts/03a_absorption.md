# 3a · Gas Absorption — building `τ_abs`

> **For:** users who care how spectral gas-absorption optical depth is computed; HITRAN power users; retrieval developers tuning line-shape choices.
>
> **Prev:** [3 · Layer Optical Properties](03_layer_optics.md) · **Next:** [3b · Mie & Rayleigh](03b_scattering.md)

This page covers one of the four arrays the layer-optics pipeline produces:
``\tau_\mathrm{abs}(\nu, z)`` — gas-absorption optical depth as a function of
wavenumber and atmospheric layer.

## What this page produces

For each spectral grid point ``\nu`` (cm⁻¹) and each atmospheric layer ``z``:

```math
\tau_\mathrm{abs}(\nu, z) = \sum_\mathrm{species}\, \sigma_\mathrm{species}(\nu, T_z, p_z)\,N_\mathrm{species}(z),
```

where ``\sigma`` is the absorption cross-section and ``N`` is the column
number density of the absorbing species in the layer. Once computed,
``\tau_\mathrm{abs}`` is added to the scattering optical depth as
``\tau_\lambda = \tau_\mathrm{abs} + \tau_\mathrm{scat}`` in
[Concepts/03c](03c_mixing.md).

## HITRAN line-by-line

vSmartMOM uses the HITRAN line list as the primary source for molecular
absorption parameters. The `Absorption` module exposes:

- `read_hitran(...)` to read a `HitranTable` from a HITRAN-format file.
- `make_hitran_model(...)` to build a `HitranModel` carrying the line list,
  broadening selector, and complex-error-function (CEF) choice.
- `compute_absorption_cross_section(model, grid, p, T)` to evaluate
  ``\sigma(\nu, T, p)`` on a wavenumber grid.

Each line in the HITRAN table contributes through its strength ``S``,
center wavenumber ``\nu_0``, lower-state energy ``E_l``, half-widths for
self- and air-broadening (``\gamma_\mathrm{self}``, ``\gamma_\mathrm{air}``),
pressure-shift coefficient, and a temperature exponent. The sum over lines
is computed on the spectral grid via a batched kernel (see below).

## Line-shape families

Three line shapes are supported:

| Shape | Physical regime | Form |
|---|---|---|
| **Doppler** (Gaussian) | low pressure, hot gas | ``\sigma(\nu) \propto \exp(-(\nu-\nu_0)^2/\gamma_d^2)`` |
| **Lorentz** | high pressure (collision-broadened) | ``\sigma(\nu) \propto \gamma_l / [(\nu-\nu_0)^2 + \gamma_l^2]`` |
| **Voigt** (convolution) | mixed | ``\mathrm{Doppler} \star \mathrm{Lorentz}`` |

Voigt is the universal choice for atmospheric retrievals — it covers the
transition between Doppler-dominated upper atmosphere and
Lorentz-dominated lower atmosphere. It is evaluated through the complex
error function (Faddeyeva function `w(z)`); vSmartMOM provides multiple CEF
implementations including a GPU-friendly Humlicek-style approximation
([`src/Absorption/complex_error_functions.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Absorption/complex_error_functions.jl)).

```
   HITRAN line list  ──→  make_hitran_model  ──┐
                                                │
   Spectral grid ν                              ▼
   Layer p, T, VMR  ─────────────────→  compute_absorption_cross_section
                                                │
                                                ▼
                                        σ(ν, T, p)
                                                │
                                                ▼
                                  τ_abs(ν, z) = σ · N(z)
                                                │
                                                ▼
                                   Layer optics  (Concepts 3)
```

## The GPU line-shape kernel

Because the line-shape sum is a per-grid-point reduction over thousands of
contributing lines, it parallelizes trivially across the spectral grid —
exactly the workload pattern of the RT solver. The Voigt kernel from
[`src/Absorption/compute_absorption_cross_section.jl:229–280`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Absorption/compute_absorption_cross_section.jl#L229-L280):

```julia
@kernel function line_shape_batch!(A, @Const(grid), @Const(ν_arr), @Const(γ_d_arr),
                                   @Const(γ_l_arr), @Const(y_arr), @Const(S_arr),
                                   @Const(istart_arr), @Const(istop_arr),
                                   N_lines, ::Voigt, CEF)
    I = @index(Global, Linear)
    FT = eltype(A)
    acc = zero(FT)
    ν_i = FT(grid[I])
    @inbounds for j in 1:N_lines
        if istart_arr[j] <= I <= istop_arr[j]
            acc += FT(S_arr[j]) * FT(cSqrtLn2divSqrtPi) / FT(γ_d_arr[j]) *
                   real(w(CEF, FT(cSqrtLn2) / FT(γ_d_arr[j]) * (ν_i - FT(ν_arr[j]))
                            + im * FT(y_arr[j])))
        end
    end
    @inbounds A[I] += acc
end
```

`@Const(...)` declares the input arrays as read-only — KernelAbstractions
uses this to enable backend-specific optimizations (constant-memory caching
on CUDA, etc.). `@index(Global, Linear)` extracts the spectral index. The
kernel runs **identically on CPU, CUDA, and Metal** via the same KA
dispatch as the RT kernels. See [Concepts/07](07_architecture.md) for the
backend story.

The window `istart_arr[j] ... istop_arr[j]` is the wavenumber sub-range
where line `j` contributes non-negligibly — pre-computed to avoid evaluating
distant lines.

## Where `τ_abs` enters the layer-optics pipeline

After the kernel produces ``\sigma(\nu, T, p)`` per layer, multiplication by
column number density (and summation across species in the same layer)
gives ``\tau_\mathrm{abs}(\nu, z)``. This is wrapped in a
`CoreAbsorptionOpticalProperties` and added to the layer-optics build chain:

```julia
# inside compEffectiveLayerProperties.jl:58 (paraphrased)
combo = rayleigh + sum(aerosols)                       # τ_scat, ϖ, Z accumulated
final = combo + CoreAbsorptionOpticalProperties(τ_abs) # τ_λ = τ_abs + τ_scat,
                                                        # ϖ_λ = τ_scat·ϖ / τ_λ
```

The `+` overload for `(CoreScatteringOpticalProperties, CoreAbsorptionOpticalProperties)`
is what bumps ``\tau`` and rescales ``\varpi``. Crucially, it does **not**
change the doubling count — see [Concepts/03c](03c_mixing.md) for the
scattering-vs-total τ split that makes this work.

## Code anchors

| Concept | Source |
|---|---|
| HITRAN parser | [`src/Absorption/read_hitran.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Absorption/read_hitran.jl) |
| Absorption types | [`src/Absorption/types.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Absorption/types.jl) |
| Cross-section computation | [`src/Absorption/compute_absorption_cross_section.jl:32–280`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Absorption/compute_absorption_cross_section.jl#L32-L280) |
| Line-shape kernel (`@kernel`) | [`src/Absorption/compute_absorption_cross_section.jl:229–280`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Absorption/compute_absorption_cross_section.jl#L229-L280) |
| Complex error functions | [`src/Absorption/complex_error_functions.jl`](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/main/src/Absorption/complex_error_functions.jl) |
| Constants & TIPS partition function | `src/Absorption/constants/` |

For practical HITRAN data management — downloading line lists, switching
between HITRAN editions, caching — see the [HITRAN Data Management](../Absorption/HITRAN_Data.md)
page in Resources.

## Hands-on tutorials

Runnable examples with Plotly figures:

- [Absorption line-shape & cross-sections](../tutorials/Tutorial_Absorption.md)

## References

- Gordon et al. (2017, 2022, 2024), HITRAN database releases. JQSRT.
- Humlíček (1979, 1982), *Optimized computation of the Voigt and complex probability functions*, JQSRT **21**:309 and **27**:437.
- Faddeyeva, V. N. & Terentiev, N. M. (1961), *Tables of the probability integral for complex argument*. (Original `w(z)` definition.)
- Crib sheet: `docs/dev_notes/theory_references.md` §F.
