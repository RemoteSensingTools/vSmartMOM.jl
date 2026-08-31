using CairoMakie
using JLD2
using OCORamanWorkflow

function resolve_grid_file()
    if haskey(ENV, "GRID_FILE")
        return ENV["GRID_FILE"]
    elseif haskey(ENV, "GRID_OUTDIR")
        return latest_jld2(ENV["GRID_OUTDIR"])
    else
        error("Set GRID_FILE=/path/to/file.jld2 or GRID_OUTDIR=/path/to/grid/output")
    end
end

grid_file = resolve_grid_file()
fig_dir = ensure_dir(get(ENV, "FIGURE_OUTDIR", joinpath(figures_root(), "quicklook")))

nu, r_rrs, r_nors = jldopen(grid_file, "r") do f
    (f["ν"], f["R_rrs"], f["R_nors"])
end

wavelength_nm = 1.0e7 ./ nu
rrs = first_lastdim_series(r_rrs)
nors = first_lastdim_series(r_nors)
delta = rrs .- nors

fig = Figure(size=(1000, 650))
ax1 = Axis(fig[1, 1], xlabel="wavelength (nm)", ylabel="radiance proxy",
    title="O2 A-band RRS vs noRS")
lines!(ax1, wavelength_nm, rrs, label="RRS")
lines!(ax1, wavelength_nm, nors, label="noRS")
axislegend(ax1, position=:rb)

ax2 = Axis(fig[2, 1], xlabel="wavelength (nm)", ylabel="RRS - noRS",
    title="Raman residual")
lines!(ax2, wavelength_nm, delta)

path = joinpath(fig_dir, splitext(basename(grid_file))[1] * "_quicklook.png")
save(path, fig)
println("Wrote ", path)
print_grid_summary(grid_file)
