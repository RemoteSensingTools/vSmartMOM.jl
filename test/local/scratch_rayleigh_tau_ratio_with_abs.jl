# Same as scratch_rayleigh_tau_ratio.jl, but adds a tunable absorption
# optical depth on top of the Rayleigh single-layer column.
#
# Hypothesis: a τ_rayl × 0.966 reduction gives ratio > 1 only in the
# presence of significant absorption (large τ_abs) and high albedo.
# Pure Rayleigh (τ_abs = 0) never exceeds 1 (see
# scratch_rayleigh_tau_ratio.jl). Mechanism: at τ_abs ≫ 1 the direct
# beam to the surface is extinguished, so the reflected radiance at
# TOA is dominated by upper-atmosphere multiple scattering; a τ_rayl
# reduction lowers that contribution. But at very high albedo the
# direct-to-surface-to-detector path (proportional to e^{-2·τ_abs/μ})
# is amplified by the τ reduction (larger e^{-2·(τ_abs)/μ} factor
# stays the same but the surface flux gain dominates), tipping the
# ratio above 1.
using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.Architectures: Architectures
using Printf, Statistics, JLD2

params = parameters_from_yaml("benchmarks/natraj.yaml")
params.architecture = Architectures.CPU()
params.vza = [0.0]
params.vaz = [0.0]

FT = Float64
scale = 0.96581   # ϖ_Cabannes equivalent

τ_rayl_list = [0.05, 0.1, 0.25, 0.5]
τ_abs_list  = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
α_list      = [0.0, 0.25, 0.8, 0.99]

function set_τ_rayl!(model, τ_target)
    τ_old = sum(model.optics.τ_rayl[1])
    model.optics.τ_rayl[1] .*= τ_target / τ_old
end

function set_τ_abs!(model, τ_target)
    # Distribute τ_abs uniformly across layers (single layer here, so
    # this is just the total).
    n_layers = size(model.optics.τ_abs[1], 2)
    model.optics.τ_abs[1] .= τ_target / n_layers
end

rows = NamedTuple[]
for τr in τ_rayl_list, τa in τ_abs_list, α in α_list
    params.brdf = CoreRT.AbstractSurfaceType[CoreRT.LambertianSurfaceScalar{FT}(FT(α))]
    model = model_from_parameters(params)

    set_τ_rayl!(model, τr)
    set_τ_abs!(model, τa)
    R_full, _ = rt_run(model)
    I_full    = Array(R_full[1, 1, :])[1]

    model.optics.τ_rayl[1] .*= scale
    R_scl, _  = rt_run(model)
    I_scl     = Array(R_scl[1, 1, :])[1]

    ratio = I_scl / I_full
    push!(rows, (; τr, τa, α, I_full, I_scaled = I_scl, ratio))
end

# Compact table per τr: rows = τa, cols = α
for τr in τ_rayl_list
    println()
    println("══ τ_rayl = $(τr) ════════════════════════════════════════════════")
    print("τ_abs\\α  ")
    for α in α_list; print(@sprintf("%10.2f   ", α)); end
    println()
    for τa in τ_abs_list
        print(@sprintf("%5.2f  ", τa))
        for α in α_list
            r = first(r for r in rows if r.τr == τr && r.τa == τa && r.α == α)
            flag = r.ratio > 1.0 ? "*" : " "
            print(@sprintf(" %8.5f%s   ", r.ratio, flag))
        end
        println()
    end
end
println("\n* marks ratio > 1 (I_scaled exceeds I_full at this τ_rayl,τ_abs,α)")
@save joinpath(@__DIR__, "scratch_rayleigh_tau_ratio_with_abs.jld2") rows τ_rayl_list τ_abs_list α_list scale
