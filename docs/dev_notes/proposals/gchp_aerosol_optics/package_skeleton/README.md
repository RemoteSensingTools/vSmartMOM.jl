# GCHPAerosolOptics.jl (skeleton)

Couples GCHP (GOCART–TOMAS) aerosol output to vSmartMOM RT. **This is a scaffold** — most
function bodies are `TODO`. Read `../STUDENT_BRIEF.md` first; it explains the design forks,
what already exists on the vSmartMOM `gchp-io` branch, and the staged roadmap.

## Dependency rule
`GCHPAerosolOptics` **depends on** `vSmartMOM`; vSmartMOM never depends back. It reuses
`vSmartMOM.Scattering` (Mie → `AerosolOptics`) rather than reimplementing Mie.

## Getting started
```julia
pkg> activate .
pkg> dev /home/cfranken/code/gitHub/vSmartMOM.jl   # record the local core in THIS env
pkg> instantiate
```
Then port the working `gchp-io` files cited in each stub (`PORT TARGET:` comments) into the
matching module, swapping `..Scattering`/`..Aerosols` for `using vSmartMOM`.

## Layout
```
src/
  GCHPAerosolOptics.jl       module + includes + exports
  io/gchp_scene.jl           GCHP column reader            (PORT: gchp-io GCHPScene.jl)
  representation/
    representation.jl        AerosolRepresentation interface (the Mode A / B fork)
    sectional.jl             Mode A (bins as-they-are; spline) (PORT + new A2)
    lognormal_fit.jl         Mode B (parametric fit)          (NEW — implement the fit)
  mixing/effective_ri.jl     speciation → effective RI; RI(λ) endpoint hook (PORT + core hook)
  optics/
    column_aod.jl            Stage 1 — bins → AOD             (PORT: gchp-io aod_diagnostics)
    layer_optics.jl          Stage 3 — sectional → AerosolOptics (NEW — core missing piece)
  ml/
    scene_dataset.jl         Stage 6 — (input, radiances) dataset
    reduction.jl             EOF/SVD input compression
test/runtests.jl             port gchp-io aerosol tests as the regression baseline
```

## Roadmap (see brief §9)
0 run the AOD example · 1 validate column AOD · 2 add spline + lognormal-fit · 3 build the
per-layer `AerosolOptics` bridge · 4 RI(λ)-across-window (core change) · 5 end-to-end RT
from one column · 6 scale to the ML dataset + EOF reduction.
