# =============================================================================
# pilot_cia.jl — Quick check: do O2–O2 and N2–N2 CIA explain the 2 μm gap?
#
# MODTRAN vs vSmartMOM at 2000 nm scenario
# (H2O=0.05, AOT=0.0001, GNDALT=0, TSZ=168):
#     transm_down_dir  MODTRAN = 0.437   vSmartMOM = 0.975
#     ⇒ slant-τ gap ≈ -ln(0.437)·μ_s - (-ln(0.975)·μ_s)
#                   ≈ 0.638 − 0.0195 ≈ 0.62 nepers
#     ⇒ vertical-τ gap ≈ 0.62 · μ_s ≈ 0.48
#
# This pilot downloads HITRAN-CIA O2-O2 and N2-N2 data, interpolates onto our
# atmospheric profile (see test/test_parameters/ParamsEMIT_MODTRANcomp.yaml)
# and integrates the collision-induced vertical optical depth at 2000 nm.
# If the resulting τ_CIA is of order 0.1–0.4, CIA is a substantial piece of
# the 2 μm gap. If it's ~0.01, CIA isn't the main culprit and we need to
# look elsewhere (baseline aerosol, MT_CKD water continuum, line-absorber
# LUT issues).
#
# Output: printed summary of τ_CIA at 2000, 2100, 2200, 2500 nm plus a
# per-layer breakdown.
# =============================================================================

using Downloads
using Printf
using vSmartMOM
using vSmartMOM.CoreRT

const CIA_DIR = expanduser("~/data/HITRAN_CIA")
# User-supplied HITRAN2024 CIA files, already present at CIA_DIR.
const CIA_FILES = Dict(
    "O2-O2" => joinpath(CIA_DIR, "O2-O2_2024.cia"),
    "O2-N2" => joinpath(CIA_DIR, "O2-N2_2024.cia"),
    "N2-N2" => joinpath(CIA_DIR, "N2-N2_2021.cia"),
)

const YAML_PATH = joinpath(pkgdir(vSmartMOM), "test", "test_parameters",
                           "ParamsEMIT_MODTRANcomp.yaml")

# Wavelengths (nm) to check: UV-Vis (O2-O2 double-transition bands),
# 1.27 μm (O2 a¹Δg), 1.58 μm (O2 a¹Δg x ν_CO2), and the 2 μm gap region.
const LAMBDAS_NM = [
    380.0, 450.0, 477.0, 532.0, 577.0, 630.0,     # UV-Vis O2-O2 double-transitions
    760.0,                                         # O2 A-band (line-dominated too)
    1060.0, 1270.0,                               # O2 a¹Δg bands
    1580.0,                                        # O2 in CO2
    1900.0, 2000.0, 2100.0, 2160.0, 2200.0,       # 2 μm gap region (N2 overtone)
    2300.0, 2500.0, 3000.0,                        # outside CIA bands
]

# ------------------------------------------------------------------
# HITRAN-CIA file parser.
#
# Block structure (from https://hitran.org/cia/):
#   line 1 (header, 100 cols):
#     cols  1–20   chemical formula ("O2-O2               ")
#     cols 21–30   ν_min  (f10.3)
#     cols 31–40   ν_max  (f10.3)
#     cols 41–47   n_pts  (i7)
#     cols 48–54   T_K    (f7.1)
#     cols 55–64   σ_max  (e10.3)
#     cols 65–…    reference / comment / dataset
#   lines 2..n_pts+1:
#     ν  σ  [error_estimate]
# ------------------------------------------------------------------
struct CIABlock
    formula::String
    T::Float64
    ν::Vector{Float64}
    σ::Vector{Float64}       # cm⁵ / molec²
end

function parse_cia_file(path::AbstractString)
    blocks = CIABlock[]
    open(path) do io
        while !eof(io)
            header = readline(io)
            if length(header) < 54 || isempty(strip(header))
                continue
            end
            formula = strip(header[1:20])
            # header cols 21-40 carry ν_min / ν_max (redundant with ν[1]/ν[end])
            n_pts   = parse(Int,     strip(header[41:47]))
            T_K     = parse(Float64, strip(header[48:54]))
            νs      = Vector{Float64}(undef, n_pts)
            σs      = Vector{Float64}(undef, n_pts)
            for i in 1:n_pts
                vals = split(strip(readline(io)))
                νs[i] = parse(Float64, vals[1])
                σs[i] = parse(Float64, vals[2])
            end
            push!(blocks, CIABlock(formula, T_K, νs, σs))
        end
    end
    return blocks
end

"""
    cia_sigma(blocks, ν_query, T_query)

Bi-linear interpolation in (T, ν) over the CIA blocks that contain `ν_query`.
Returns σ in cm⁵/molec² (0.0 if ν_query is outside all blocks at this T).
"""
function cia_sigma(blocks::Vector{CIABlock}, ν_query::Real, T_query::Real)
    # Group blocks into temperature slices and find two bracketing T.
    Ts = sort!(unique([b.T for b in blocks]))
    isempty(Ts) && return 0.0
    T_lo = Ts[1]; T_hi = Ts[end]
    for i in 1:length(Ts) - 1
        if Ts[i] <= T_query <= Ts[i + 1]
            T_lo = Ts[i]; T_hi = Ts[i + 1]; break
        end
    end
    if T_query < Ts[1];   T_lo = T_hi = Ts[1];   end
    if T_query > Ts[end]; T_lo = T_hi = Ts[end]; end

    σ_at_T = function (Ts_pick)
        # Pick the block at this temperature that brackets ν_query.
        for b in blocks
            b.T == Ts_pick || continue
            b.ν[1] <= ν_query <= b.ν[end] || continue
            j = searchsortedfirst(b.ν, ν_query)
            if j <= 1
                return b.σ[1]
            elseif j > length(b.ν)
                return b.σ[end]
            else
                ν1, ν2 = b.ν[j - 1], b.ν[j]
                s1, s2 = b.σ[j - 1], b.σ[j]
                w = (ν_query - ν1) / (ν2 - ν1)
                return (1 - w) * s1 + w * s2
            end
        end
        return 0.0
    end

    σ_lo = σ_at_T(T_lo)
    σ_hi = σ_at_T(T_hi)
    (σ_lo == 0.0 && σ_hi == 0.0) && return 0.0
    if T_hi == T_lo
        return σ_lo
    end
    w = (T_query - T_lo) / (T_hi - T_lo)
    return (1 - w) * σ_lo + w * σ_hi
end

# ------------------------------------------------------------------
# Profile + number-density helpers
# ------------------------------------------------------------------
const k_B = 1.380649e-23     # J/K
const amagat_ref_n = 2.6867811e19  # molec / cm³ at 273.15 K, 1 atm

"""
    number_density(p_hPa, T_K) -> molec/cm³
"""
number_density(p_hPa::Real, T_K::Real) = 1e2 * p_hPa / (k_B * T_K) * 1e-6
# 1e2 converts hPa → Pa, 1e-6 converts m⁻³ → cm⁻³

"""
    layer_thickness_cm(p_top, p_bot, T) -> cm  (hydrostatic, H=R·T/(M·g))
"""
function layer_thickness_cm(p_top_hPa, p_bot_hPa, T_K)
    R = 8.31446  # J/(mol·K)
    M = 28.9647e-3  # kg/mol (dry air)
    g = 9.80665  # m/s²
    H = R * T_K / (M * g)  # m
    Δz_m = H * log(p_bot_hPa / p_top_hPa)
    return Δz_m * 100.0  # cm
end

function load_profile()
    params = parameters_from_yaml(YAML_PATH)
    p_half = Float64.(params.p)   # TOA→BOA half-levels (hPa), length N+1
    T_mid  = Float64.(params.T)   # layer centres (K), length N
    vmr_O2 = params.absorption_params.vmr["O2"]     # scalar 0.209
    vmr_N2 = 1.0 - 0.209 - 3.30e-4 - mean_h2o_vmr(params)  # approx
    return (; p_half, T_mid, vmr_O2, vmr_N2, sza = params.sza)
end

function mean_h2o_vmr(params)
    v = params.absorption_params.vmr["H2O"]
    v isa AbstractVector ? sum(v) / length(v) : v
end

# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
function pilot()
    for (name, path) in CIA_FILES
        isfile(path) || error("Missing CIA file: $path")
    end
    @info "Loading CIA files" CIA_FILES
    blocks_O2O2 = parse_cia_file(CIA_FILES["O2-O2"])
    blocks_O2N2 = parse_cia_file(CIA_FILES["O2-N2"])
    blocks_N2N2 = parse_cia_file(CIA_FILES["N2-N2"])
    for (name, bks) in [("O2-O2", blocks_O2O2), ("O2-N2", blocks_O2N2), ("N2-N2", blocks_N2N2)]
        @info "  $name" n_blocks = length(bks) ν_span = (minimum(b.ν[1] for b in bks), maximum(b.ν[end] for b in bks)) T_span = extrema([b.T for b in bks])
    end

    prof = load_profile()
    μ_s  = cos(deg2rad(Float64(prof.sza)))

    # Vertical CIA optical depth per wavelength, with HITRAN Eq. (3):
    #   k(ν) = k_AA(ν,T)·ρ_A² + k_AB(ν,T)·ρ_A·ρ_B + k_BB(ν,T)·ρ_B²
    @printf("\n  %-8s  %-11s  %-11s  %-11s  %-11s  %-12s  %-12s\n",
            "λ_nm", "τ_O2O2", "τ_O2N2", "τ_N2N2", "τ_CIA_vert", "τ_CIA_slant", "T_dir_slant")
    println("-"^100)

    for λ in LAMBDAS_NM
        ν = 1e7 / λ   # cm⁻¹
        τ_O2O2 = 0.0
        τ_O2N2 = 0.0
        τ_N2N2 = 0.0
        for iz in 1:length(prof.T_mid)
            p_top = prof.p_half[iz]
            p_bot = prof.p_half[iz + 1]
            T     = prof.T_mid[iz]
            p_mid = 0.5 * (p_top + p_bot)
            n_air = number_density(p_mid, T)
            n_O2  = prof.vmr_O2 * n_air
            n_N2  = prof.vmr_N2 * n_air
            Δz    = layer_thickness_cm(p_top, p_bot, T)
            τ_O2O2 += cia_sigma(blocks_O2O2, ν, T) * n_O2 * n_O2 * Δz
            τ_O2N2 += cia_sigma(blocks_O2N2, ν, T) * n_O2 * n_N2 * Δz
            τ_N2N2 += cia_sigma(blocks_N2N2, ν, T) * n_N2 * n_N2 * Δz
        end
        τ_vert  = τ_O2O2 + τ_O2N2 + τ_N2N2
        τ_slant = τ_vert / μ_s
        T_slant = exp(-τ_slant)
        @printf("  %-8.1f  %-11.3e  %-11.3e  %-11.3e  %-11.3e  %-12.3e  %-12.6f\n",
                λ, τ_O2O2, τ_O2N2, τ_N2N2, τ_vert, τ_slant, T_slant)
    end

    @info "Target Δτ (slant) to explain MODTRAN gap at μ_s=$(round(μ_s; digits=4)):"
    @info "   λ=2000 nm: ≈ 0.62   (⇒ Δτ_vert ≈ 0.48)"
    @info "   λ=2100 nm: ≈ 0.21   (⇒ Δτ_vert ≈ 0.16)"
    @info "   λ=2200 nm: ≈ 0.33   (⇒ Δτ_vert ≈ 0.25)"
    @info "N2 2.16 μm overtone CIA band is at 4300–5000 cm⁻¹  (λ = 2.00–2.33 μm)"
    @info "O2 1.27 μm band is at 7545–8355 cm⁻¹            (λ = 1.20–1.32 μm)"
end

if abspath(PROGRAM_FILE) == @__FILE__
    pilot()
end
