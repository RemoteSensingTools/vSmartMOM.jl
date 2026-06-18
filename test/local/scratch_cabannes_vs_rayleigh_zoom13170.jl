# Zoom around the absorption-free reference at ν = 13170 cm⁻¹.
# Reads the three zero@13170 JLD2s (full / tauonly / greekonly), one
# per mode, and produces:
#   1. A 3-panel figure (full, tauonly, greekonly) zoomed to a few
#      cm⁻¹ around 13170, 100·ΔI vs ν for 6 subset albedos.
#   2. A scatter of 100·ΔI vs albedo at exactly ν = 13170 cm⁻¹ for the
#      three modes overlaid (one curve per mode).
pushfirst!(LOAD_PATH, joinpath(homedir(), ".julia", "environments",
                                "v$(VERSION.major).$(VERSION.minor)"))
using Plots, JLD2, Printf, Statistics

ν_lo, ν_hi  = 13169.0, 13171.0          # zoom window
ν_zero      = 13170.0                    # the absorption-free pixel
modes       = ["full", "tauonly", "greekonly"]
subset_α    = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

# ── Load all three modes ────────────────────────────────────────────
function load_mode(m)
    jld = joinpath(@__DIR__,
        "scratch_cabannes_vs_rayleigh_I_$(m)_ns8_vaz180_zero13170.jld2")
    d = load(jld)
    return (; results = d["results"], albedos = d["albedos"], ν = d["ν"])
end
data = Dict(m => load_mode(m) for m in modes)
ν = data[modes[1]].ν
albedos = data[modes[1]].albedos
mask    = (ν .>= ν_lo) .& (ν .<= ν_hi)
@info "Loaded" nAlbedos=length(albedos) nSpec_zoom=count(mask) ν_zero

# ── Subset indices ─────────────────────────────────────────────────
sub_idx = [findfirst(α -> isapprox(α, t; atol = 1e-6), albedos) for t in subset_α]
@assert all(!isnothing, sub_idx)

cs = cgrad(:viridis)
col_at(j, n) = cs[(j - 1) / max(1, n - 1)]
flip_diff(r) = -r.diff_S  # Cabannes − Rayleigh

# ── 3-panel zoom figure ────────────────────────────────────────────
plt = plot(layout = (3, 1), size = (1100, 950), legend = :outertopright,
           link = :x, xlims = (ν_lo, ν_hi),
           plot_title = "Cabannes − Rayleigh zoom @ ν=13170 cm⁻¹ " *
                        "(absorption-free pixel)\n" *
                        "Stokes_I, SZA=19°, VZA=0°, VAZ=180°, ns=8, NoTrunc")

for (k, m) in enumerate(modes)
    sub = plt[k]
    ylabel!(sub, "100·ΔI ($m)")
    r = data[m].results
    for (j, idx) in enumerate(sub_idx)
        rj = r[idx]
        y  = 100 .* flip_diff(rj)[1, mask]
        plot!(sub, ν[mask], y;
              lw = 1.2, marker = :circle, ms = 2, color = col_at(j, length(sub_idx)),
              label = k == 1 ? @sprintf("α = %.1f", rj.α) : "")
    end
    # Mark the zeroed pixel
    vline!(sub, [ν_zero]; color = :red, ls = :dash, lw = 1, label = k == 1 ? "ν = 13170" : "")
end
xlabel!(plt[3], "Wavenumber [cm⁻¹]")
out_png1 = joinpath(@__DIR__, "scratch_cabannes_vs_rayleigh_zoom13170.png")
savefig(plt, out_png1)
@info "Saved" out_png1

# ── ΔI at ν=13170 vs albedo, all three modes overlaid ──────────────
idx0 = argmin(abs.(ν .- ν_zero))
@assert ν[idx0] == ν_zero
plt2 = plot(size = (1000, 650), legend = :topleft,
            xlabel = "Lambertian albedo",
            ylabel = "100·ΔI = 100·[I(Cabannes) − I(Rayleigh)] at ν=13170",
            title  = "Per-mode Cabannes − Rayleigh signal at the absorption-free pixel")
markers = Dict("full" => :circle, "tauonly" => :diamond, "greekonly" => :square)
colors  = Dict("full" => :black,  "tauonly" => :red,     "greekonly" => :blue)
for m in modes
    r = data[m].results
    y = [100 * -ri.diff_S[1, idx0] for ri in r]
    plot!(plt2, albedos, y; lw = 2, marker = markers[m], ms = 5,
          color = colors[m], label = m)
end
hline!(plt2, [0]; color = :gray, ls = :dash, label = "")
out_png2 = joinpath(@__DIR__, "scratch_cabannes_vs_rayleigh_zoom13170_vs_albedo.png")
savefig(plt2, out_png2)
@info "Saved" out_png2

# ── Print compact table ────────────────────────────────────────────
@info "100·ΔI at ν=$(ν_zero) cm⁻¹, per mode and albedo"
println(@sprintf("%-5s  %12s  %12s  %12s", "α", "full", "tauonly", "greekonly"))
for (k, α) in enumerate(albedos)
    vals = (m -> 100 * -data[m].results[k].diff_S[1, idx0]).(modes)
    println(@sprintf("%5.2f  %12.4e  %12.4e  %12.4e", α, vals...))
end
