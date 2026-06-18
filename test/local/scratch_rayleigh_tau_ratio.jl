# Does a small τ_rayl reduction flip the I_scaled/I_full ratio above 1
# at high albedo on a pure Rayleigh atmosphere?
#
# Builds a single-layer pure-Rayleigh model from the Natraj-CDS YAML
# (depol = 0, no aerosols, no absorption), then for each (τ_rayl,
# albedo) pair runs RT twice:
#   1. with τ_rayl_full       — the YAML-built value
#   2. with τ_rayl_full × 0.966 — mimicking the Cabannes elastic
#                                  fraction (ϖ_Cabannes ≈ 0.966).
#
# Sweeps τ ∈ {0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1.0}, albedo ∈ {0.0,
# 0.25, 0.8}, SZA = the Natraj YAML value (cos = 1.0 case included via
# extra sweep), VZA = 0° (nadir), VAZ = 0°.

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.Architectures: Architectures
using Printf, Statistics, JLD2

# ── Build a base model from the Natraj benchmark YAML ───────────────
params = parameters_from_yaml("benchmarks/natraj.yaml")
params.architecture = Architectures.CPU()  # tiny problem, CPU is fine
@info "Natraj YAML loaded" depol=params.depol nstreams=params.nstreams

# Override geometry: nadir VZA, principal-plane, SZA fixed (≈23°)
params.vza = [0.0]
params.vaz = [0.0]

FT = Float64
scale = 0.96581   # vSmartMOM's ϖ_Cabannes for this Rayleigh model

τ_list = [0.02, 0.05, 0.1, 0.15, 0.25, 0.5, 1.0]
α_list = [0.0, 0.25, 0.5, 0.8, 0.99]

# Helper: rescale τ_rayl in-place to a given total optical depth
function set_total_τ!(model, τ_target)
    τ_old = sum(model.optics.τ_rayl[1])   # vector × layers — single layer here
    factor = τ_target / τ_old
    model.optics.τ_rayl[1] .*= factor
    return τ_old, factor
end

# ── Sweep ─────────────────────────────────────────────────────────────
rows = NamedTuple[]
for τ in τ_list, α in α_list
    params.brdf = CoreRT.AbstractSurfaceType[CoreRT.LambertianSurfaceScalar{FT}(FT(α))]
    model = model_from_parameters(params)

    # Force total Rayleigh OD to the requested τ (the YAML's two-layer
    # column gives some baseline τ_rayl that we override here).
    τ_before, _ = set_total_τ!(model, τ)
    R_full, _   = rt_run(model)
    I_full      = Array(R_full[1, 1, :])[1]    # nadir VZA, Stokes I

    # Scaled run: τ_rayl × 0.96581
    model.optics.τ_rayl[1] .*= scale
    R_scaled, _ = rt_run(model)
    I_scaled    = Array(R_scaled[1, 1, :])[1]

    ratio = I_scaled / I_full
    push!(rows, (; τ, α, I_full, I_scaled, ratio,
                    ratio_minus_1 = ratio - 1))
    @info @sprintf("τ=%.3f α=%.2f  I_full=%.4e  I_scaled=%.4e  ratio=%.5f  Δ=%+.3e",
                   τ, α, I_full, I_scaled, ratio, ratio - 1)
end

# ── Table ─────────────────────────────────────────────────────────────
println()
println("══════════════════════════════════════════════════════════════════════")
println("I_scaled(τ·0.966) / I_full(τ)  — pure Rayleigh, nadir, SZA = $(params.sza)°")
println("──────────────────────────────────────────────────────────────────────")
print("τ\\α  ")
for α in α_list
    print(@sprintf("%10.2f   ", α))
end
println()
for τ in τ_list
    print(@sprintf("%5.3f", τ))
    for α in α_list
        r = first(r for r in rows if r.τ == τ && r.α == α)
        flag = r.ratio > 1.0 ? "*" : " "
        print(@sprintf(" %8.5f%s   ", r.ratio, flag))
    end
    println()
end
println("══════════════════════════════════════════════════════════════════════")
println("* marks ratio > 1 (I_scaled exceeds I_full)")

# ── Save ──────────────────────────────────────────────────────────────
out = joinpath(@__DIR__, "scratch_rayleigh_tau_ratio.jld2")
@save out rows τ_list α_list scale
@info "Saved" out
