# The central design fork (STUDENT_BRIEF §4): how the 15-bin GCHP distribution is turned
# into something optics can be computed from. All representations share one interface so
# the optics/RT code consumes any of them identically.

"""
    AerosolRepresentation

How a per-layer aerosol size distribution is represented for optics:

  - `Sectional`      — use the 15 bins as they are (Mode A1). PORTED from gchp-io.
  - `SplineSectional`— smooth spline through bin centers, then quadrature (Mode A2). NEW.
  - `LogNormalFit`   — fit a (sum of) lognormal(s) per layer (Mode B). NEW (stub on gchp-io).

Speciation is carried per-bin regardless (the hybrid recommended in the brief: parametric
size shape + sectional composition).
"""
abstract type AerosolRepresentation end

"""
    layer_distributions(rep, scene_aerosols, ilev) -> distribution(s) + per-mode effective RI inputs

Contract every representation implements. Returns whatever `layer_aerosol_optics` needs:
for `Sectional`, the per-bin number density + per-bin composition; for `LogNormalFit`, the
fitted `LogNormal` params + the effective composition. Keep the return type uniform enough
that `optics/layer_optics.jl` dispatches cleanly.
"""
function layer_distributions end
