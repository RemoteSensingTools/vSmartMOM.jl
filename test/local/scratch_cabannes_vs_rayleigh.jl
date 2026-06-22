# =================================================================
# Rayleigh vs Cabannes elastic-cross-section comparison.
#
# Two purely-elastic (noRS) runs of the O2 A-band over Lambertian
# surface albedos 0:0.05:1.0:
#
#   Run A — "Rayleigh": full Rayleigh greek coefficients with the
#           total τ_rayl (the conventional noRS default — implicitly
#           includes the inelastic redistribution via the effective
#           depolarization).
#
#   Run B — "Cabannes": elastic-only substitution that mirrors how
#           the RRS code splits Rayleigh into Cabannes + Raman:
#             * greek_rayleigh → greek_cabannes (lower-depol phase
#               matrix, ~3% shift in β/δ)
#             * τ_rayl         → τ_rayl · ϖ_Cabannes (~0.96 elastic
#               fraction; the (1-ϖ_Cab) part is the Raman source we
#               would add back with an RRS run).
#
# The difference R_Rayleigh − R_Cabannes per albedo is exactly the
# elastic bias one would carry if running RRS with the inelastic
# source omitted (or vice versa, if not switching to Cabannes when
# turning RRS on).
#
# Geometry: SZA = 19°, VZA = 0°, VAZ = 0°.
# Polarization: Stokes_IQU.
# Aerosols: none (τ_ref = 0).
# Spectral: (1e7/773):0.02:(1e7/757) cm⁻¹, full O2 A-band.
#
# Run from the test/ directory:
#     cd test && julia --project=. scratch_cabannes_vs_rayleigh.jl
# =================================================================

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.Architectures: Architectures
using vSmartMOM.Scattering: Stokes_I, Stokes_IQU, Stokes_IQUV
using Statistics
using Printf
using JLD2

# ---- Polarization-mode dispatch (ENV-driven) ----------------------
# POL=I    → Stokes_I
# POL=IQUV → Stokes_IQUV
# (default IQU, matches earlier runs)
pol_tag = uppercase(get(ENV, "POL", "IQU"))
pol_type, pol_label = if pol_tag == "I"
    (Stokes_I(), "I")
elseif pol_tag == "IQUV"
    (Stokes_IQUV(), "IQUV")
elseif pol_tag == "IQU"
    (Stokes_IQU(), "IQU")
else
    error("Unknown POL=$pol_tag (use I, IQU, or IQUV)")
end
@info "Polarization mode" pol_tag pol_label nStokes=pol_type.n

# ---- Substitution-mode dispatch (ENV-driven) ----------------------
# MODE=full      → greek_rayleigh → greek_cabannes  AND  τ_rayl → τ_rayl·ϖ_Cab
# MODE=tauonly   → τ_rayl → τ_rayl·ϖ_Cab only;        greek stays full Rayleigh
# MODE=greekonly → greek_rayleigh → greek_cabannes only; τ_rayl unchanged
mode_tag = lowercase(get(ENV, "MODE", "full"))
@assert mode_tag in ("full", "tauonly", "greekonly") "MODE must be 'full', 'tauonly', or 'greekonly'"
@info "Substitution mode" mode_tag
swap_greek  = mode_tag in ("full", "greekonly")
scale_τrayl = mode_tag in ("full", "tauonly")

# Load CUDA before model_from_parameters so the CUDAExt extension
# (ext/vSmartMOMCUDAExt.jl) is registered and `array_type(GPU())`
# returns `CuArray`. Without this, GPU() is selected but kernels still
# get CPU `Array`s and silently run on CPU.
using CUDA
@info "CUDA functional" CUDA.functional() devices=collect(CUDA.devices())

# ---- Spectral window: 13110–13180 cm⁻¹, the band-line zoom we
# settled on for the subset plot. Narrower than the full O2 A-band
# but still includes the strong R-branch lines and the clean
# continuum at 13166–13180 cm⁻¹.
spec  = 13110.0:0.02:13180.0
nSpec = length(spec)
@info "Spectral grid" ν_lo=first(spec) ν_hi=last(spec) Δν=step(spec) nSpec

# ---- Albedo sweep --------------------------------------------------
albedos = collect(0.0:0.05:1.0)
@info "Albedo sweep" albedos

# ---- Build base params --------------------------------------------
# Use the O2-only RRS YAML as the base (single band, O2 absorption,
# small placeholder aerosol that we zero below). Override spec_bands,
# geometry, polarization, single surface.
params = parameters_from_yaml("test_parameters/O2ABandParamsRRS.yaml")
params.spec_bands = [collect(spec)]
params.sza        = 19.0
params.vza        = [0.0]
params.vaz        = [parse(Float64, get(ENV, "VAZ", "0.0"))]
@info "Geometry" sza=params.sza vza=params.vza vaz=params.vaz
params.polarization_type = pol_type
# Disable phase-function truncation — keep the full Rayleigh/Cabannes
# greek coefficients all the way through to the Mie/RT path.
params.truncation = vSmartMOM.Scattering.NoTruncation()
@info "Truncation forced to NoTruncation()"
# Override nstreams from ENV (default keeps YAML value of 3). Public
# contract: stream_l_cap = 2·nstreams - 1, l_trunc = stream_l_cap,
# max_m = stream_l_cap + 1. parameters_from_yaml derives these at
# parse time, so they have to be updated alongside `nstreams` here.
nstreams_override = parse(Int, get(ENV, "NSTREAMS", "0"))
if nstreams_override > 0
    @info "Overriding nstreams" old_nstreams=params.nstreams new_nstreams=nstreams_override old_stream_l_cap=params.stream_l_cap
    params.nstreams     = nstreams_override
    params.stream_l_cap = 2 * nstreams_override - 1
    params.l_trunc      = params.stream_l_cap
    params.max_m        = params.stream_l_cap + 1
    @info "Derived stream caps after override" stream_l_cap=params.stream_l_cap l_trunc=params.l_trunc max_m=params.max_m
end

FT = Float64
params.brdf = CoreRT.AbstractSurfaceType[CoreRT.LambertianSurfaceScalar{FT}(0.0)]
params.architecture = Architectures.GPU()
# Override the long 49-layer profile from O2ABandParamsRRS.yaml with
# the 34-level (33-layer) US Standard atmosphere preset from
# 2BandParameters.yaml. Required because the 49-layer profile's
# original full-level pressures top out at ~956 hPa, so reducing to 20
# uniform half-levels extrapolates past the original grid and trips an
# Interpolations BoundsError. The 34-level preset reduces cleanly.
params.T = [231.62, 244.33, 251.34, 258.09, 264.25, 269.15,
            272.59, 274.07, 273.30, 269.65, 264.27, 258.11,
            251.52, 245.22, 239.20, 234.05, 229.71, 225.70,
            222.70, 220.62, 219.32, 217.93, 216.98, 217.10,
            218.35, 223.33, 234.19, 249.34, 264.12, 277.20,
            280.77, 282.60, 284.40, 285.80]
params.p = [  0.14,   0.22,   0.30,   0.39,   0.53,   0.71,
              0.96,   1.28,   1.70,   2.27,   3.03,   4.03,
              5.44,   7.26,   9.67,  12.90,  17.23,  23.30,
             31.00,  42.07,  56.09,  74.78,  99.69, 131.00,
            176.85, 236.64, 314.58, 418.87, 557.76, 735.00,
            800.12, 849.00, 912.00, 980.00, 1005.0]
params.q = zeros(eltype(params.T), length(params.T))  # dry atmosphere (no H2O VMR)
params.profile_reduction_n = 20

# Zero aerosols (keep one placeholder entry but with τ_ref=0).
if params.scattering_params !== nothing
    for ra in params.scattering_params.rt_aerosols
        ra.τ_ref = 0.0
    end
end

@info "Building Rayleigh model (full greek_rayleigh, τ_rayl as-is)"
model = model_from_parameters(params)
iBand = 1
ν     = model.atmosphere.spec_bands[iBand]

# ---- Optional: zero O2 absorption at one spectral grid point ------
# `ZERO_NU` env var (cm⁻¹) → find the closest sample in ν and set
# τ_abs[iBand][idx, :] = 0 so that ONE wavelength is absorption-free
# (pure Rayleigh+surface). Useful as a clean reference inside an
# otherwise absorbing band.
zero_nu = parse(Float64, get(ENV, "ZERO_NU", "0.0"))
if zero_nu > 0.0
    idx = argmin(abs.(ν .- zero_nu))
    @info "Zeroing O2 absorption at one grid point" requested=zero_nu actual=ν[idx] idx τ_before=sum(model.optics.τ_abs[iBand][idx, :])
    model.optics.τ_abs[iBand][idx, :] .= 0
    @info "After zero" τ_after=sum(model.optics.τ_abs[iBand][idx, :])
end
@info "Model built" Nlayers=size(model.optics.τ_rayl[1], 2) ϖ_Cabannes=model.optics.rayleigh.ϖ_Cabannes[iBand]

# ---- Construct the Cabannes-substituted "shadow" model ------------
# - Replace the greek_rayleigh slot with greek_cabannes so the noRS
#   dispatch (`_rayleigh_greek_source(::noRS, gr, gc) = gr`) picks up
#   the Cabannes phase matrix.
# - Scale τ_rayl by ϖ_Cabannes to drop the inelastic-fraction OD.
# Everything else (τ_abs, profile, aerosols, surfaces vector, sources)
# is shared by reference with the Rayleigh model — only optics
# differs.
ϖ_cab = model.optics.rayleigh.ϖ_Cabannes  # Vector{FT}, per band
cab_greek_slot = swap_greek ? model.optics.rayleigh.greek_cabannes :
                              model.optics.rayleigh.greek_rayleigh
cab_rayleigh = CoreRT.RayleighScattering(
    cab_greek_slot,
    model.optics.rayleigh.greek_cabannes,
    ϖ_cab,
)
scaled_τ_rayl = scale_τrayl ?
    [model.optics.τ_rayl[i] .* ϖ_cab[i] for i in 1:length(model.optics.τ_rayl)] :
    model.optics.τ_rayl
cab_optics = CoreRT.Optics(cab_rayleigh, model.optics.aerosols,
                            model.optics.τ_abs, scaled_τ_rayl)
cab_surfaces = CoreRT.AbstractSurfaceType[CoreRT.LambertianSurfaceScalar{FT}(0.0)]
cab_model = CoreRT.RTModel(model.architecture, model.solver, model.numerics,
                            model.geometry, model.quad_points, model.atmosphere,
                            cab_optics, cab_surfaces, model.sources)
@info "Cabannes shadow model built" τ_rayl_ratio=mean(scaled_τ_rayl[1] ./ model.optics.τ_rayl[1])

# ---- Loop over albedos --------------------------------------------
results = Vector{NamedTuple}(undef, length(albedos))
nVZA  = length(params.vza)
nStokes = CoreRT.polarization_type(model).n
@info "Running $(length(albedos)) × 2 RT runs (Stokes_$nStokes, nVZA=$nVZA, nSpec=$nSpec)…"

t0 = time()
for (k, α) in enumerate(albedos)
    @info @sprintf("── α = %.2f (%d/%d) ──", α, k, length(albedos))
    model.surfaces[1]     = CoreRT.LambertianSurfaceScalar{FT}(FT(α))
    cab_model.surfaces[1] = CoreRT.LambertianSurfaceScalar{FT}(FT(α))

    @info "  Rayleigh run…"
    @time R_ray, T_ray = rt_run(model; i_band = iBand)
    @info "  Cabannes run…"
    @time R_cab, T_cab = rt_run(cab_model; i_band = iBand)

    # R has shape (nVZA, nStokes, nSpec). Pick first VZA and bring to host
    # (R is a CuArray when running on GPU). Keep the full Stokes slab so
    # IQUV runs preserve Q/U/V — slice (nStokes, nSpec).
    S_ray = Array(R_ray[1, :, :])
    S_cab = Array(R_cab[1, :, :])
    diff_S = S_ray .- S_cab  # (nStokes, nSpec)

    I_ray = view(S_ray, 1, :); I_cab = view(S_cab, 1, :)
    diff_I = view(diff_S, 1, :)
    @info @sprintf("  I:    Rayleigh⟨⟩=%.4e  Cabannes⟨⟩=%.4e  Δ⟨⟩=%.3e  |Δ|max=%.3e",
                   mean(I_ray), mean(I_cab), mean(diff_I), maximum(abs, diff_I))

    results[k] = (; α, S_ray, S_cab, diff_S)
end
@info @sprintf("Total elapsed: %.1f s", time() - t0)

# ---- Save (polarization+mode+streams+vaz+zero-tagged filename) ----
ns_suffix   = nstreams_override > 0 ? "_ns$(nstreams_override)" : ""
vaz_suffix  = params.vaz[1] != 0.0 ? "_vaz$(Int(params.vaz[1]))" : ""
zero_suffix = zero_nu > 0.0 ? "_zero$(Int(round(zero_nu)))" : ""
out_jld2 = joinpath(@__DIR__,
    "scratch_cabannes_vs_rayleigh_$(pol_label)_$(mode_tag)$(ns_suffix)$(vaz_suffix)$(zero_suffix).jld2")
@save out_jld2 results albedos ν ϖ_cab pol_label mode_tag zero_nu
@info "Saved → $out_jld2"

# Plotting is handled by scratch_cabannes_vs_rayleigh_plot.jl
# (uses the global v$(VERSION) environment which has Plots).
