# =============================================================================
# plot_ramanXS.jl — dump Cabannes + rotational-Raman cross-sections for N₂
# and O₂ (Phase 6 port from sanghavi).
#
# Extracts the molecular cross-section coefficients that drive the RRS path
# (via getRamanAtmoConstants) and writes four `.dat` files suitable for
# comparison against published RRS cross-section references.
#
# The sanghavi version also produced diagnostic plots via Plots.jl; those
# are dropped in the port. Callers who want plots can load the .dat output
# in any plotting environment. Optional PNG emission behind
# `PLOT_RAMANXS_PLOTS=1` + having `Plots` available in the active env.
# =============================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Statistics
using DelimitedFiles
using Printf

# Where to write cross-section `.dat` outputs. Override with env var.
const RAMANXS_OUTDIR = get(ENV, "RAMANXS_OUTDIR",
                            joinpath(@__DIR__, "ramanXS_output"))

# Load noabs SIF-grid config (noRS + RRS both work on this YAML)
parameters = parameters_from_yaml(joinpath(@__DIR__, "..",
                                  "test_parameters/noabs_parameters2_SIF_grid.yaml"))
parameters.architecture = vSmartMOM.Architectures.CPU()
model = model_from_parameters(parameters)

const FT = CoreRT.float_type(model)

# Zero aerosol, push Rayleigh optical thickness into the last-but-one layer
# (matches sanghavi's pre-conditioning for flat-profile RRS studies).
model.τ_aer[1][1, :, :]     .= 0
model.τ_aer[1][1, :, end]   .= 100
model.τ_rayl[1][:, end - 1] .+= model.τ_rayl[1][:, end]
model.τ_rayl[1][:, end]      .= 0

const iBand = 1
ν  = model.atmosphere.spec_bands[iBand]
ν̃  = FT(0.5) * (model.atmosphere.spec_bands[1][1] +
                model.atmosphere.spec_bands[end][end])
effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)
n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

nPol  = CoreRT.polarization_type(model).n
nSpec = length(ν)
F₀   = zeros(FT, nPol, nSpec); F₀[1, :] .= one(FT)
SIF₀ = zeros(FT, nPol, nSpec)

RS_type_rrs = InelasticScattering.RRS(
    n2          = n2,
    o2          = o2,
    greek_raman = InelasticScattering.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
    fscattRayl  = [FT(1)],
    ϖ_Cabannes  = [FT(1)],
    ϖ_λ₁λ₀      = zeros(FT, 1),
    i_λ₁λ₀      = zeros(Int, 1),
    Z⁻⁺_λ₁λ₀    = zeros(FT, 1, 1),
    Z⁺⁺_λ₁λ₀    = zeros(FT, 1, 1),
    i_ref       = 0,
    n_Raman     = 0,
    F₀          = F₀,
    SIF₀        = SIF₀,
)
CoreRT.getRamanSSProp!(RS_type_rrs, 1e7 / ν̃, ν)

mkpath(RAMANXS_OUTDIR)

# Pure Cabannes cross-section coefficients (scale by ν⁴ for the actual XS).
writedlm(joinpath(RAMANXS_OUTDIR, "effCoeff_Cabannes_N2.dat"), [0 n2.effCoeff.σ_Rayl_coeff])
writedlm(joinpath(RAMANXS_OUTDIR, "effCoeff_Cabannes_O2.dat"), [0 o2.effCoeff.σ_Rayl_coeff])

# Rotational Raman transitions (J → J+2 and J → J-2)
Δν_N2    = vcat(n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2,
                n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2)
kscatt_N2 = vcat(n2.effCoeff.σ_RoRaman_coeff_JtoJp2,
                 n2.effCoeff.σ_RoRaman_coeff_JtoJm2)
Δν_O2    = vcat(o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2,
                o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2)
kscatt_O2 = vcat(o2.effCoeff.σ_RoRaman_coeff_JtoJp2,
                 o2.effCoeff.σ_RoRaman_coeff_JtoJm2)

writedlm(joinpath(RAMANXS_OUTDIR, "effCoeff_RRS_N2.dat"), [Δν_N2 kscatt_N2])
writedlm(joinpath(RAMANXS_OUTDIR, "effCoeff_RRS_O2.dat"), [Δν_O2 kscatt_O2])

@info "plot_ramanXS: wrote Cabannes + RRS cross-section tables" dir = RAMANXS_OUTDIR

# Optional plotting — requires Plots in the active env. Kept as a conditional
# side effect so the script still works on a CPU-only / headless env.
if get(ENV, "PLOT_RAMANXS_PLOTS", "0") == "1"
    try
        @eval begin
            using Plots
            ν₀ = ν̃
            p_n2 = plot([ν₀, ν₀],
                        [n2.effCoeff.σ_Rayl_coeff * ν₀^4, n2.effCoeff.σ_Rayl_coeff * ν₀^4] .* 1e-3,
                        seriestype = :sticks, lw = 2, label = "N₂ σₑ(x10⁻³)", color = :red,
                        yformatter = x -> string(@sprintf("%.2e", x)),
                        xlabel = "wavenumber [cm⁻¹]", ylabel = "Cross sections [cm²]")
            plot!(p_n2, ν₀ .+ n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2,
                  n2.effCoeff.σ_RoRaman_coeff_JtoJp2 .* (ν₀ .+ n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4,
                  seriestype = :sticks, lw = 2, color = :blue, label = "")
            plot!(p_n2, ν₀ .+ n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2,
                  n2.effCoeff.σ_RoRaman_coeff_JtoJm2 .* (ν₀ .+ n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4,
                  seriestype = :sticks, lw = 2, color = :blue, label = "N₂ σᵢₑ")
            savefig(p_n2, joinpath(RAMANXS_OUTDIR, "n2_XS.png"))
        end
        @info "plot_ramanXS: wrote N₂ plot" file = joinpath(RAMANXS_OUTDIR, "n2_XS.png")
    catch err
        @warn "plot_ramanXS: plotting skipped (Plots not available?)" err
    end
end
