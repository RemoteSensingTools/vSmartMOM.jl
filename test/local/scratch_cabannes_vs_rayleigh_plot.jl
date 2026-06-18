# Re-plot Cabannes − Rayleigh comparison from saved JLD2.
#
# Reads scratch_cabannes_vs_rayleigh_$(pol_label).jld2, where pol_label
# is set via ENV["POL"] (I or IQUV; default IQU). Produces, per
# polarization, an absolute-difference subset plot scaled ×100, one
# panel per Stokes component:
#   I     for POL=I
#   I,Q,U for POL=IQU
#   I,Q,U,V for POL=IQUV
# Subset albedos: 0.0, 0.2, 0.4, 0.6, 0.8, 1.0.
pushfirst!(LOAD_PATH, joinpath(homedir(), ".julia", "environments",
                                "v$(VERSION.major).$(VERSION.minor)"))
using Plots, JLD2, Printf, Statistics

pol_label = uppercase(get(ENV, "POL", "IQU"))
mode_tag  = lowercase(get(ENV, "MODE", "full"))
nstreams_tag = parse(Int, get(ENV, "NSTREAMS", "0"))
vaz_tag      = parse(Float64, get(ENV, "VAZ", "0.0"))
zero_tag     = parse(Float64, get(ENV, "ZERO_NU", "0.0"))
kind         = lowercase(get(ENV, "KIND", "diff"))   # "diff" or "ratio"
@assert kind in ("diff", "ratio") "KIND must be 'diff' or 'ratio'"
ns_suffix    = nstreams_tag > 0 ? "_ns$(nstreams_tag)" : ""
vaz_suffix   = vaz_tag != 0.0 ? "_vaz$(Int(vaz_tag))" : ""
zero_suffix  = zero_tag > 0.0 ? "_zero$(Int(round(zero_tag)))" : ""
kind_suffix  = kind == "ratio" ? "_ratio" : ""
@info "Polarization / mode / streams / vaz / zero_nu / kind" pol_label mode_tag nstreams_tag vaz_tag zero_tag kind

jld = joinpath(@__DIR__,
    "scratch_cabannes_vs_rayleigh_$(pol_label)_$(mode_tag)$(ns_suffix)$(vaz_suffix)$(zero_suffix).jld2")
@info "Loading" jld
data = load(jld)
results = data["results"]
albedos = data["albedos"]
ν       = data["ν"]
ϖ_cab   = data["ϖ_cab"]

nStokes = size(results[1].S_ray, 1)
stokes_names = ("I", "Q", "U", "V")[1:nStokes]
@info "Loaded" nAlbedos=length(results) nSpec=length(ν) nStokes ϖ_Cab=round(ϖ_cab[1], digits=5)

# Sign convention: ΔS = S(Cabannes) − S(Rayleigh).
flip_diff(r) = -r.diff_S  # results saved as Rayleigh − Cabannes

# ── Subset of albedos to display ────────────────────────────────────
# - ALL_ALBEDOS=1 → every saved albedo (21 lines)
# - ALBEDO_STEP=Δ → 0:Δ:1.0 step (default 0.2 → 6 lines)
all_albedos = get(ENV, "ALL_ALBEDOS", "0") == "1"
α_step      = parse(Float64, get(ENV, "ALBEDO_STEP", "0.2"))
subset_indices = if all_albedos
    collect(1:length(albedos))
else
    targets = collect(0.0:α_step:1.0)
    inds = [findfirst(α -> isapprox(α, t; atol = 1e-6), albedos) for t in targets]
    @assert all(!isnothing, inds) "Missing one of $targets in saved albedos"
    inds
end
@info "Albedo set" all_albedos α_step n_shown=length(subset_indices)

# ── One subplot per Stokes component ────────────────────────────────
cs_sub = cgrad(:viridis)
ν_lo_disp, ν_hi_disp = 13110.0, 13180.0  # match the saved sim window
disp_mask = (ν .>= ν_lo_disp) .& (ν .<= ν_hi_disp)

title_prefix = kind == "ratio" ? "Cabannes / Rayleigh" : "Cabannes − Rayleigh"
ylabel_fmt   = kind == "ratio" ? (s) -> "$s(Cabannes) / $s(Rayleigh)" :
                                  (s) -> "100·Δ$s"
plt = plot(layout = (nStokes, 1), size = (1200, 350 * nStokes),
           legend = :outertopright, link = :x,
           xlims = (ν_lo_disp, ν_hi_disp),
           left_margin = 12Plots.mm,
           plot_title = "$title_prefix, Stokes_$pol_label, " *
                        "mode=$mode_tag — SZA=19°, VZA=0°, " *
                        @sprintf("VAZ=%.0f°, ns=%d, no aerosols",
                                 vaz_tag, nstreams_tag > 0 ? nstreams_tag : 3),
           plot_titlefontsize = 10,
           titlefontsize     = 9,
           guidefontsize     = 8,
           tickfontsize      = 8,
           legendfontsize    = 6)

for (s, sname) in enumerate(stokes_names)
    sub = plt[s]
    ylabel!(sub, ylabel_fmt(sname))
    for (j, idx) in enumerate(subset_indices)
        r = results[idx]
        col = cs_sub[(j - 1) / max(1, length(subset_indices) - 1)]
        y = if kind == "ratio"
            # I(Cabannes) / I(Rayleigh) = (S_ray − diff_S) / S_ray  (where
            # diff_S was stored as Rayleigh − Cabannes).
            S_ray = r.S_ray[s, disp_mask]
            S_cab = S_ray .- r.diff_S[s, disp_mask]
            S_cab ./ S_ray
        else
            100 .* flip_diff(r)[s, disp_mask]
        end
        plot!(sub, ν[disp_mask], y; lw = 0.9, color = col,
              label = s == 1 ? @sprintf("α = %.1f", r.α) : "")
    end
    # For ratio, draw the y=1 reference line.
    kind == "ratio" && hline!(sub, [1.0]; color = :gray, ls = :dash, lw = 0.7, label = "")
end
xlabel!(plt[nStokes], "Wavenumber [cm⁻¹]")
all_suffix = if all_albedos
    "_allα"
elseif α_step != 0.2
    "_step$(replace(string(α_step), "." => "p"))"
else
    ""
end
out_base = joinpath(@__DIR__,
    "scratch_cabannes_vs_rayleigh_$(pol_label)_$(mode_tag)$(ns_suffix)$(vaz_suffix)$(zero_suffix)$(kind_suffix)$(all_suffix)")
for ext in ("png", "pdf")
    out_file = out_base * "." * ext
    savefig(plt, out_file)
    @info "Saved" out_file
end

# ── Mean-bias summary (Cabannes − Rayleigh, ×100, I only) ───────────
@info "Mean-bias summary (Cabannes − Rayleigh, I only)"
for r in results
    mb = mean(-r.diff_S[1, :])
    mi = mean(r.S_ray[1, :])
    println(@sprintf("  α=%.2f  100·⟨ΔI⟩=% .4e  ⟨I_Rayleigh⟩=%.3e  rel=% .3f%%",
                     r.α, 100*mb, mi, 100*mb/mi))
end
