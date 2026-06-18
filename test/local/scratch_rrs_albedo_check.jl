# =================================================================
# RRS perturbation vs albedo: vSmartMOM smoke test
#
# Mirrors the "Inelastic - Elastic" intercomparison plot (vSmartMOM vs
# VLRRS vs SCIATRAN) at SZA=19° in a small O2-A window (13110–13185 cm⁻¹)
# at three Lambertian albedos: 0.0, 0.2, 0.8.
#
# For each albedo we compute (R_RRS + ieR) − R_noRS, i.e. the Ring-effect
# perturbation that the reference panels plot as "Inelastic − Elastic".
#
# Saves results to test/scratch_rrs_albedo_check.jld2 and emits a PNG
# (via global-env Plots.jl) at test/scratch_rrs_albedo_check.png.
#
# Run from the test/ directory:
#     cd test && julia --project=. scratch_rrs_albedo_check.jl
# =================================================================

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Statistics
using JLD2

# ---- Spectral window -----------------------------------------------
# Simulation must cover the RRS source range: rotational Raman shifts for
# N₂/O₂ at 300 K extend to ~±150 cm⁻¹. Display only the central O2-A
# region (13138–13144 cm⁻¹) so per-line filling-in is resolved.
ν_lo        = 13041.0    # ~100 cm⁻¹ buffer below display
ν_hi        = 13241.0    # ~100 cm⁻¹ buffer above display
Δν          = 0.1        # cm⁻¹ — resolves O2 line cores at this pressure
spec        = ν_lo:Δν:ν_hi
nSpec       = length(spec)
disp_lo     = 13138.0    # plotting window low
disp_hi     = 13144.0    # plotting window high
@info "Spectral grid (sim)" ν_lo ν_hi Δν nSpec disp_lo disp_hi

# ---- Build base params ---------------------------------------------
# Start from O2ABandParamsRRS.yaml (includes O2 absorption + Rayleigh +
# tiny aerosol with τ_ref = 0) and override:
#   - spec_bands  -> narrow window
#   - SZA, VZA    -> match reference plot (SZA=19°, nadir VZA=0°)
#   - surfaces    -> swapped per-albedo in the loop
#   - profile_reduction kept at YAML default (=1, ie no reduction);
#     drop to a smaller number if wall-time is too long.
params = parameters_from_yaml("test_parameters/O2ABandParamsRRS.yaml")
params.spec_bands = [collect(spec)]
params.sza        = 19.0
params.vza        = [0.0]
params.vaz        = [0.0]

# Use Stokes_I to keep runtime tractable (the reference plot shows I only).
params.polarization_type = vSmartMOM.Scattering.Stokes_I()

# Single surface placeholder; gets swapped per-albedo below.
FT = Float64
params.brdf = CoreRT.AbstractSurfaceType[CoreRT.LambertianSurfaceScalar{FT}(0.0)]

model = model_from_parameters(params)
iBand = 1
ν     = model.atmosphere.spec_bands[iBand]
ν̃    = mean(ν)
nPol  = CoreRT.polarization_type(model).n

# ---- Optionally zero absorption to match a pure Rayleigh+RRS reference scene
# (matches what VLRRS / SCIATRAN intercomparisons typically use). Toggle via
# the `NO_ABS=1` env var — without it the original O2-A absorption is kept.
if get(ENV, "NO_ABS", "0") == "1"
    @info "NO_ABS=1 → zeroing τ_abs (pure Rayleigh + RRS scene)"
    for τ in model.optics.τ_abs
        fill!(τ, zero(eltype(τ)))
    end
end

# ---- Build RS_type (RRS) and noRS_type templates -------------------
effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

F₀ = zeros(FT, nPol, nSpec); F₀[1, :] .= 1
SIF₀ = zeros(FT, nPol, nSpec)

make_RRS() = begin
    rs = InelasticScattering.RRS(
        n2 = n2, o2 = o2,
        greek_raman = InelasticScattering.GreekCoefs(
            [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
        fscattRayl  = [FT(1)], ϖ_Cabannes = [FT(1)],
        ϖ_λ₁λ₀ = zeros(FT, 1), i_λ₁λ₀ = zeros(Int, 1),
        Z⁻⁺_λ₁λ₀ = zeros(FT, 1, 1), Z⁺⁺_λ₁λ₀ = zeros(FT, 1, 1),
        i_ref = argmin(abs.(ν .- ν̃)),
        n_Raman = 0,
        F₀ = copy(F₀), SIF₀ = copy(SIF₀))
    CoreRT.getRamanSSProp!(rs, 1e7/ν̃, ν)
    rs
end

make_noRS() = begin
    nrs = InelasticScattering.noRS(
        fscattRayl  = [FT(1)], ϖ_Cabannes = [FT(1)],
        bandSpecLim = [], iBand = [1],
        F₀ = zeros(FT, nPol, nSpec))
    nrs.F₀[1, :] .= 1
    nrs
end

# ---- Loop over albedos --------------------------------------------
albedos = [0.0, 0.2, 0.8]
results = Dict{Float64, NamedTuple}()

for α in albedos
    @info "── Running α = $α ──"
    model.surfaces[1] = CoreRT.LambertianSurfaceScalar{FT}(FT(α))

    @info "  noRS run…"
    @time R_noRS, T_noRS, _, _ = CoreRT.rt_run_test(make_noRS(), model, iBand)

    @info "  RRS run…"
    @time R_rrs, T_rrs, ieR, ieT = CoreRT.rt_run_test(make_RRS(), model, iBand)

    # Stokes-I only (vza=1, pol=1)
    I_noRS  = R_noRS[1, 1, :]
    I_total = R_rrs[1, 1, :] .+ ieR[1, 1, :]
    pert    = I_total .- I_noRS

    @info "  α = $α: |pert| max=$(maximum(abs, pert))  mean=$(mean(abs, pert))  median I_noRS=$(median(I_noRS))"
    results[α] = (; ν = collect(ν), I_noRS, I_total, pert)
end

# ---- Save -----------------------------------------------------------
_suffix = get(ENV, "NO_ABS", "0") == "1" ? "_noabs" : ""
out_jld2 = joinpath(@__DIR__, "scratch_rrs_albedo_check$(_suffix).jld2")
@save out_jld2 results albedos ν_lo ν_hi Δν
@info "Saved → $out_jld2"

# ---- Plot via global-env Plots (best-effort) -----------------------
try
    # Look in the default depot's @v#.# environment for Plots.
    pushfirst!(LOAD_PATH, joinpath(homedir(), ".julia", "environments",
                                   "v$(VERSION.major).$(VERSION.minor)"))
    @eval using Plots
    plt = plot(layout = (3, 1), size = (1100, 850),
               link = :x, legend = :topright,
               xlabel = "", ylabel = "Stokes-I I/F reflectance × 100",
               plot_title = "RRS perturbation, Inelastic − Elastic (vSmartMOM, SZA=19°)")
    for (k, α) in enumerate(albedos)
        r = results[α]
        mask = (r.ν .>= disp_lo) .& (r.ν .<= disp_hi)
        plot!(plt[k], r.ν[mask], 100 .* r.pert[mask]; lw=1.2, color=:blue,
              label = "vSmartMOM",
              title = "albedo = $α, SZA = 19")
    end
    xlabel!(plt[3], "Wavenumber [1/cm]")
    out_png = joinpath(@__DIR__, "scratch_rrs_albedo_check$(_suffix).png")
    savefig(plt, out_png)
    @info "Saved → $out_png"
catch err
    @warn "Plotting skipped: $err"
end
