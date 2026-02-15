```@meta
EditURL = "Tutorial_CoreRT.jl"
```

# CoreRT Starter Tutorial: Elemental -> Doubling -> Interaction

Goal: understand the execution workflow of the core matrix-operator solver and how it maps to code.

This tutorial is intentionally lightweight. It shows:
1. how to run a minimal RT case,
2. where the core operators are applied, and
3. how to connect runtime steps to equations in the references.

## 1) Minimal run

````@example Tutorial_CoreRT
using vSmartMOM
using vSmartMOM.CoreRT
````

Start from package defaults and make small runtime choices.

````@example Tutorial_CoreRT
params = default_parameters()
params.max_m = 2
params.l_trunc = 20
params.vza = [0.0, 30.0, 60.0]
params.vaz = [0.0, 0.0, 0.0]
````

Build the model and run one spectral band.

````@example Tutorial_CoreRT
model = model_from_parameters(params)
result = rt_run(model, i_band = 1)
````

result[1] is top-of-atmosphere reflected field (SFI branch).

````@example Tutorial_CoreRT
R_SFI = result[1]
size(R_SFI)
````

## 2) What happens internally

The layer loop is implemented in `src/CoreRT/rt_run.jl`.
For each Fourier mode and each atmospheric layer, `rt_kernel!` in
`src/CoreRT/CoreKernel/rt_kernel.jl` applies:

1. `elemental!`:
   builds elemental `r/t/j` operators from layer optical properties.
2. `doubling!`:
   recursively doubles the elemental layer to the target optical thickness.
3. `interaction!`:
   adds the new (doubled) layer to the existing composite atmosphere.

For the first layer (`iz == 1`), the added layer is copied directly into
the composite layer (no interaction yet). For `iz > 1`, adding is applied.

## 3) Equation provenance quick map

- Elemental layer: infinitesimal/elemental operator form in
  Sanghavi et al. (2014), Eq. (19) family.
- Doubling/adding recursions: Sanghavi et al. (2014), Eq. (23)-(28)
  (scalar precursor in Sanghavi et al. (2013), Eq. (22)-(27)).
- D-matrix symmetry (`diag(1,1,-1,-1)` per stream):
  Sanghavi et al. (2014), Eq. (29)-(32).
- SFI single-scattering source terms (`J₀⁺`, `J₀⁻`) in `get_elem_rt_SFI!`:
  Fell (1997), Eq. (1.52)-(1.54), with thin-layer limits Eq. (1.55)-(1.56).
- Truncation updates used by scattering optics:
  Sanghavi and Stephens (2015), Eq. (38a)-(38d).

For a detailed symbol-to-code mapping, see the documentation page:
`vSmartMOM -> Core RT Theory (Doubling/Adding)`.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

