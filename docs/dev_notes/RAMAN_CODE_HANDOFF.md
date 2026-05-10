# How the code runs now — for the Raman code author

This note explains the current unified RT flow so you can check whether the Raman path still works. The branch is `unified-vsmartmom`; it merged analytic Jacobians (linearized RT) with the existing forward RT. **Raman is only used in the forward path; the linearized path is elastic-only (noRS) for now.**

---

## 1. Two separate paths: forward vs linearized

| Path | Purpose | Raman? | Entry point |
|------|--------|--------|-------------|
| **Forward** | Compute R, T (radiance / flux) | **Yes** — you pass `RS_type` | `rt_run(RS_type, model, iBand)` or `rt_run_test(RS_type, model, iBand)` |
| **Linearized** | Compute R, T **and** analytic dR/d(params) | **No** — always `noRS()` | `rt_run(model, lin_model, NAer, NGas, NSurf; i_band=1)` |

So:

- **Your Raman code is only used in the forward path.** The linearized path never takes an `RS_type` and never calls Raman routines.
- To test Raman, you only need to use the **forward** entry points and pass a Raman `RS_type` (e.g. `RRS`, `VS_0to1`).

---

## 2. Forward path (where Raman runs)

**Entry (elastic default):**

```julia
rt_run(model; i_band=1)   # same as rt_run(noRS(), model, 1)
```

**Entry (Raman — your path):**

```julia
rt_run_test(RS_type, model, iBand)   # RS_type = RRS(), VS_0to1(), etc.
```

This dispatches to `rt_run(RS_type, model, iBand)` in **`src/CoreRT/rt_run.jl`**.

**What that function does:**

1. **Setup**  
   Uses `RS_type` to build layers:
   - `make_added_layer(RS_type, FT_dual, arr_type, dims, nSpec)`  
   - `make_composite_layer(RS_type, ...)`  
   For Raman types (`RRS`, `VS_0to1_plus`, `VS_1to0_plus`), `rt_helper_functions.jl` returns `AddedLayerRS` / composite with the extra Raman dimension (`n_Raman`).

2. **Per Fourier moment `m`**  
   - `InelasticScattering.computeRamanZλ!(RS_type, pol_type, Array(qp_μ), m, arr_type)`  
   - `constructCoreOpticalProperties(RS_type, iBand, m, model)`  
   - For each layer `iz`:  
     - If Raman: `RS_type.fscattRayl = expandBandScalars(RS_type, fScattRayleigh[iz])`  
     - `expandOpticalProperties(layer_opt_props[iz], arr_type)`  
     - **`rt_kernel!(RS_type, pol_type, SFI, added_layer, composite_layer, layer_opt, scattering_interface, τ_sum, m, quad_points, I_static, architecture, qp_μN, iz)`**

3. **Surface and postprocess**  
   - `create_surface_layer!(...)`  
   - `interaction!(RS_type, ...)`  
   - `postprocessing_vza!(RS_type, ...)` → fills `R`, `T`, etc.

So in the forward path, **Raman is wired through** `RS_type` into layer creation, optical properties, `rt_kernel!`, interaction, and postprocessing.

**Where the Raman-specific kernel lives:**  
For non-noRS types, `rt_kernel!` is implemented in **`src/CoreRT/CoreKernel/rt_kernel.jl`** (and in multisensor / single-scatter variants). Those implementations use:

- `elemental_inelastic!(RS_type, ...)`  
- `elemental!(...)` (elastic part)  
- `doubling_inelastic!(RS_type, ...)`  
- `interaction!(RS_type, ...)`  

So your Raman logic is in:

- **CoreRT:** `rt_run.jl` (forward driver), `rt_kernel.jl`, `rt_kernel_ss.jl`, `rt_kernel_multisensor.jl`, `rt_helper_functions.jl` (layer types for Raman), and any `interaction!` / `postprocessing_vza!` that dispatch on `RS_type`.
- **Inelastic:** `src/Inelastic/` (types, `computeRamanZλ!`, `getRamanSSProp!`, inelastic helpers, etc.).

The unified code did **not** change the forward `rt_run(RS_type, model, iBand)` signature or the fact that `RS_type` is passed all the way into `rt_kernel!` and down to elemental/doubling/interaction. So the **intended** design is that Raman still runs when you call `rt_run_test(RS_type, model, iBand)` with a Raman type.

---

## 3. Linearized path (no Raman)

**Entry:**

```julia
model, lin_model = model_from_parameters(LinMode(), params)
R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf; i_band=1)
```

Internally this always calls **`rt_run(noRS(), model, lin_model, NAer, NGas, NSurf, i_band)`** (see **`src/CoreRT/rt_run_lin.jl`**). So:

- No `RS_type` argument: the linearized run is **elastic-only**.
- Layer creation uses **`make_added_layer(lin, RS_type, ...)`** and **`make_composite_layer(lin, RS_type, ...)`** with `RS_type = noRS()`. In **`rt_helper_functions_lin.jl`** only **`noRS` / `noRS_plus`** have linearized overloads; there are **no** `AddedLayerRS`-style linearized constructors for Raman.
- The kernel used is **`rt_kernel!(RS_type::noRS, ...)`** in **`src/CoreRT/CoreKernel/rt_kernel_lin.jl`**, which uses:
  - `elemental!` from **`elemental_lin.jl`** (elastic linearized),
  - `lin_added_layer_all_params!`,
  - `doubling_allparams!`,
  - then interaction and postprocessing.

So the **linearized path never touches Raman**. If you only care about “does Raman forward still work?”, you can ignore the linearized path.

---

## 4. Single-scatter (rt_run_ss)

**Forward single-scatter (Raman possible):**

- **`rt_run_ss(model; i_band=1)`** → `rt_run_ss(noRS(), model, i_band)`  
- **`rt_run_test_ss(RS_type, model, iBand)`** → `rt_run_ss(RS_type, model, iBand)` (Raman path)

`rt_run_ss(RS_type, ...)` is implemented in **`src/CoreRT/rt_run_lin.jl`** (same file as the linearized run; name is historical). It:

- Uses **`InelasticScattering.computeRamanZλ!(RS_type, ...)`** when building the run,
- Uses **`rt_kernel_ss!(RS_type, ...)`** (which has Raman branches: `elemental_inelastic!`, etc.),
- Then `create_surface_layer!`, then `postprocessing_vza!`.

So **single-scatter Raman** is still invoked when you call **`rt_run_test_ss(RS_type, model, iBand)`** with a Raman type. Again, the linearized part of the same file is noRS-only and does not affect this.

---

## 5. What to check as the Raman author

1. **Forward multiple scattering with Raman**  
   Call **`rt_run_test(RS_type, model, iBand)`** with your usual Raman `RS_type` and a `model` built from parameters (e.g. `model = model_from_parameters(params)` **without** `LinMode()`).  
   - Confirm layer creation uses `AddedLayerRS` (or your Raman composite type).  
   - Confirm `rt_kernel!` is the one that calls `elemental_inelastic!` and `doubling_inelastic!` (not the noRS-only `rt_kernel_lin.jl`).

2. **Forward single-scatter with Raman**  
   Call **`rt_run_test_ss(RS_type, model, iBand)`** and confirm `rt_kernel_ss!` and `computeRamanZλ!` run as before.

3. **No regressions from unified merge**  
   - In **`rt_run.jl`**, the loop over `m` and `iz` still passes `RS_type` into `rt_kernel!`, `interaction!`, and `postprocessing_vza!`.  
   - In **`rt_helper_functions.jl`**, `make_added_layer` / `make_composite_layer` for `RRS`, `VS_0to1_plus`, `VS_1to0_plus` are unchanged (they still return `AddedLayerRS` and the composite with `n_Raman`).
   So the only way Raman would break is if something in that forward call chain was changed to assume `noRS` or to drop `RS_type`. A quick grep for `RS_type` in **`rt_run.jl`** and **`rt_kernel.jl`** will show that dispatch is still in place.

4. **Linearized vs Raman**  
   You can ignore **`rt_run(model, lin_model, ...)`** and **`rt_kernel_lin.jl`** for Raman; they are explicitly noRS-only. If later someone adds linearized Raman, they would need to add linearized layer constructors for Raman in **`rt_helper_functions_lin.jl`** and a Raman branch in **`rt_kernel_lin.jl`** (and related chain-rule/doubling code).

---

## 6. Short file map (Raman-relevant)

| File | Role |
|------|------|
| **CoreRT/rt_run.jl** | Forward multi-scatter driver; calls `rt_kernel!(RS_type, ...)` and Raman setup. |
| **CoreRT/rt_run_lin.jl** | Defines `rt_run(RS_type, model, lin_model, ...)` (linearized, currently noRS-only) and `rt_run_ss(RS_type, ...)` (forward single-scatter; Raman via `rt_run_test_ss`). |
| **CoreRT/CoreKernel/rt_kernel.jl** | Forward kernel: `noRS` vs Raman; Raman uses `elemental_inelastic!`, `doubling_inelastic!`. |
| **CoreRT/CoreKernel/rt_kernel_ss.jl** | Single-scatter forward kernel with Raman. |
| **CoreRT/CoreKernel/rt_kernel_multisensor.jl** | Multisensor forward kernel with Raman. |
| **CoreRT/tools/rt_helper_functions.jl** | `make_added_layer` / `make_composite_layer` for noRS vs Raman (`AddedLayerRS`, etc.). |
| **CoreRT/tools/rt_helper_functions_lin.jl** | Linearized layers only for noRS (no Raman). |
| **CoreRT/CoreKernel/elemental_lin.jl** | Linearized elemental (noRS only); comments reference Raman equations but code path is elastic. |
| **CoreRT/CoreKernel/rt_kernel_lin.jl** | Linearized kernel; only `rt_kernel!(::noRS, ...)`. |
| **Inelastic/** | Raman types, `computeRamanZλ!`, `getRamanSSProp!`, optical/cross-section helpers. |

If you want, the next step is to add a small “Raman smoke test” (e.g. one call to `rt_run_test(RS_type, model, iBand)` and check that R/T are finite and non-zero) and document it in this repo so the unified code keeps Raman regression-free.
