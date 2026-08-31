# SS 91% error — resolved 2026-04-21

## Root cause

`interaction_ss!` (and `interaction_inelastic_ss!` variants) in
[src/CoreRT/CoreKernel/interaction_ss.jl](../src/CoreRT/CoreKernel/interaction_ss.jl)
rebound the destructured composite-layer fields to **fresh copies** before
passing them to the kernel:

```julia
(; J₀⁺, J₀⁻) = composite_layer   # aliases (correct)
...
J₀⁺ = arr_type(J₀⁺)              # rebinds locals to COPIES
J₀⁻ = arr_type(J₀⁻)
...
kernel!(..., J₀⁺, J₀⁻, ...)      # writes to the copies
```

`Array(A::Array{T,3})` returns a new array (verified: `Array(A) === A` is
false for 3-D). So the kernel's `+=` accumulation targeted a throwaway copy
and the composite layer never saw any contribution past the TOA direct-copy
path in `rt_kernel_ss!` (`iz == 1` branch).

With 12 layers, only the top layer's `j₀⁻` survived — the expected 12-layer
accumulation was reduced to 1 term, explaining the ~1/12 ≈ 0.085 ratio
observed (unified `R_SFI[1,1,50] = 2.348e-4` vs sanghavi `2.752e-3`).

## Fix

Pass `composite_layer.J₀⁺` / `.J₀⁻` (and `.ieJ₀⁺` / `.ieJ₀⁻` for the RRS /
VS variants) directly to the kernel. Removed the dead `J₀⁺/J₀⁻ = arr_type(...)`
rebindings and the now-unused destructures. Matches sanghavi's pattern.

## Verification

After fix, unified produces:
- Per-layer `composite.J₀⁻` max grows 0.00752 → 0.0806 across 12 layers
  (sanghavi's final value: 0.0806 — exact match).
- `R_SFI[1,1,50] = 2.7518e-3` vs sanghavi `2.752e-3`, ratio 0.99994
  (Float32 rounding noise).
- `test/test_forward_ss.jl` — 28/28 tests pass.

## Stale baselines

The captured SS baseline at
`test/benchmarks/baseline_output/sanghavi-unified_3f876c9/phase1b_ss_761-764nm_cpu_output.jld2`
(and `_b4a5ab1/`) was frozen before the fix, so it reflects the broken
1-layer-only accumulation. Needs re-capture on a post-fix commit before it
can serve as a comparison reference.

The Phase 2b per-Stokes tolerance gate (commit eece23a) never
pass/fail-gated on SS outputs — the harness is comparison-only, not a
regression test — so no automated test is broken by the change.

## Hypotheses that turned out wrong

The WIP analysis had fingered the 17-arg `elemental!` method's uppercase
`J₀⁺/J₀⁻` destructure against a lowercase-field `AddedLayer`. That method
**is** buggy (would throw at dispatch), but `rt_run_ss` routes through the
11-arg `CoreScatteringOpticalProperties` overload which correctly uses
lowercase `j₀⁺/j₀⁻`, so the 17-arg method is dead code on the SS path.
That code smell still wants cleanup but wasn't the cause here.
