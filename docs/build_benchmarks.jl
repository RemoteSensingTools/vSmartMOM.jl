# =============================================================================
# docs/build_benchmarks.jl
#
# Re-runs vSmartMOM against published RT benchmarks and the VLIDORT 2.8.3
# regression-test tables, captures per-cell results (modeled, truth,
# absolute error, relative error %), and writes the markdown page
# `docs/src/pages/benchmarks.md`. CPU only, F64 + F32.
#
# Run manually when reference data, parameters, or the RT solver change:
#
#     julia --project=docs docs/build_benchmarks.jl
#
# Provenance distinction (clearly stated in the page output):
#   - 📚 Natraj 2009 + Siewert 2000 PROBLEM IIA = published, theoretically
#     derived high-accuracy benchmarks. Treat any disagreement as a
#     vSmartMOM defect.
#   - 🧪 VLIDORT 2.8.3 solar_tester (Cases B / C) = saved_results
#     regression tables shipped with VLIDORT — comparing two RT codes
#     against each other, not against analytic theory.
# =============================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering
using Statistics
using Printf

const REPO_ROOT      = normpath(joinpath(@__DIR__, ".."))
const VLIDORT_DIR    = joinpath(REPO_ROOT, "test", "vlidort_baseline")
const VLIDORT_REF    = joinpath(VLIDORT_DIR, "reference_data")
const VLIDORT_CONFIG = joinpath(VLIDORT_DIR, "configs")
const NATRAJ_DIR     = joinpath(REPO_ROOT, "test", "benchmarks")
const OUT_PATH       = joinpath(@__DIR__, "src", "pages", "benchmarks.md")

# -----------------------------------------------------------------------------
# Result records: one Row per (case, stokes-component, float-type, geometry).
# -----------------------------------------------------------------------------
struct Row
    geometry  :: String          # "μ=0.02, ϕ=0°", "vza=10°, raz=0°, TOA-up", …
    modeled   :: Float64
    truth     :: Float64
end

abs_err(r::Row) = abs(r.modeled - r.truth)
function rel_err(r::Row; near_zero_atol::Real = 0.0)
    abs(r.truth) <= near_zero_atol && return NaN
    return abs(r.modeled - r.truth) / abs(r.truth)
end

# -----------------------------------------------------------------------------
# Markdown helpers
# -----------------------------------------------------------------------------
fmt_sci(x; d::Int = 3) = isnan(x) ? "NaN" : isinf(x) ? "∞" : @sprintf("%.*e", d, x)
fmt_pct(x; d::Int = 2) = isnan(x) ? "—" : @sprintf("%.*f%%", d, 100x)
fmt_val(x; d::Int = 4) = @sprintf("%.*f", d, x)

function write_table(io::IO, rows::AbstractVector{Row};
                     near_zero_atol::Real = 0.0)
    println(io, "| Geometry | Modeled | Truth | |Δ| | |Δ|/|truth| |")
    println(io, "|---|---|---|---|---|")
    for r in rows
        re = rel_err(r; near_zero_atol)
        println(io, "| ", r.geometry,
                " | ", fmt_val(r.modeled),
                " | ", fmt_val(r.truth),
                " | ", fmt_sci(abs_err(r)),
                " | ", fmt_pct(re), " |")
    end
end

function summary_line(rows::AbstractVector{Row}; near_zero_atol::Real = 0.0)
    res = [rel_err(r; near_zero_atol) for r in rows]
    valid_idx = findall(!isnan, res)
    valid = res[valid_idx]
    n_total = length(rows)
    n_valid = length(valid)
    if n_valid == 0
        return @sprintf("(n=%d, all values masked as near-zero truth)", n_total)
    end
    abs_errs_valid = [abs_err(rows[i]) for i in valid_idx]
    excluded = n_total - n_valid
    suffix = excluded == 0 ? "" :
             @sprintf("  (excluded %d cells with |truth| < %s)", excluded, fmt_sci(near_zero_atol))
    return @sprintf("n=%d  median rel = %s  max rel = %s  median |Δ| = %s  max |Δ| = %s%s",
                    n_valid,
                    fmt_pct(Statistics.median(valid)),
                    fmt_pct(maximum(valid)),
                    fmt_sci(Statistics.median(abs_errs_valid)),
                    fmt_sci(maximum(abs_errs_valid)),
                    suffix)
end

# =============================================================================
# Case A — Siewert 2000 PROBLEM IIA  (📚 published benchmark)
# =============================================================================

include(joinpath(VLIDORT_REF, "siewert2000_IIA_greek.jl"))
include(joinpath(VLIDORT_REF, "siewert2000_IIA_truth.jl"))

const SIEWERT_YAML    = joinpath(VLIDORT_CONFIG, "siewert2000_IIA.yaml")
const SIEWERT_τ_TOTAL = 1.0
const SIEWERT_ω̃      = 0.973527
const SIEWERT_VZA_COSINES = [cosd(0.0001), cosd(25.841932763), cosd(36.869897646),
                             cosd(45.572995999), cosd(53.130102354), cosd(60.0),
                             cosd(66.421821522), cosd(72.542396876), cosd(78.463040967),
                             cosd(84.260829523), cosd(89.9999)]

function siewert_aerosol_optics()
    greek = Scattering.GreekCoefs(SIEWERT_α, SIEWERT_β, SIEWERT_γ,
                                  SIEWERT_δ, SIEWERT_ϵ, SIEWERT_ζ)
    return Scattering.AerosolOptics(greek_coefs = greek, ω̃ = SIEWERT_ω̃,
                                    k = 1.0, fᵗ = 0.0)
end

function siewert_run_at_az(az_deg::Real, FT::Type)
    params = parameters_from_yaml(SIEWERT_YAML)
    params.architecture = vSmartMOM.Architectures.CPU()
    params.float_type = FT
    fill!(params.vaz, float(az_deg))
    model = model_from_parameters(params)
    model.aerosol_optics[1][1] = siewert_aerosol_optics()
    model.τ_aer[1][1, :] .= SIEWERT_τ_TOTAL
    model.τ_rayl[1] .= 0.0
    return Array(CoreRT.rt_run(model, i_band = 1)[1])
end

function siewert_pick_toa_upwelling(table::AbstractMatrix, vza_cosines)
    out = similar(vza_cosines, Float64)
    for (i, μ) in enumerate(vza_cosines)
        target = -abs(μ)
        idx = argmin(abs.(SIEWERT_TABLE2_COSINES .- target))
        out[i] = table[idx, 1]
    end
    return out
end

const _STOKES_IDX = Dict(:I => 1, :Q => 2, :U => 3, :V => 4)

function run_siewert(FT::Type)
    by_stokes = Dict{Symbol, Vector{Row}}(:I => [], :Q => [], :U => [], :V => [])
    for az_deg in (0.0, 90.0, 180.0)
        L = siewert_run_at_az(az_deg, FT)
        for stokes_sym in (:I, :Q, :U, :V)
            haskey(SIEWERT_TRUTH, (az_deg, stokes_sym)) || continue
            truth_table = SIEWERT_TRUTH[(az_deg, stokes_sym)]
            modeled = π .* L[:, _STOKES_IDX[stokes_sym], 1]
            truth = siewert_pick_toa_upwelling(truth_table, SIEWERT_VZA_COSINES)
            # vSmartMOM is Hovenier, VLIDORT/Siewert tables are Mishchenko →
            # Q/U/V truth needs sign flip for apples-to-apples comparison.
            stokes_sym in (:Q, :U, :V) && (truth = .-truth)
            for (i, μ) in enumerate(SIEWERT_VZA_COSINES)
                vza = round(acosd(μ), digits = 1)
                push!(by_stokes[stokes_sym], Row(
                    "vza=$(vza)°, az=$(round(Int, az_deg))°",
                    Float64(modeled[i]), Float64(truth[i])))
            end
        end
    end
    return by_stokes
end

# =============================================================================
# Natraj 2009  (📚 published benchmark)
# =============================================================================

include(joinpath(NATRAJ_DIR, "natraj_trues.jl"))   # I_trues, Q_trues, U_trues  (16 μ × 7 ϕ)

const NATRAJ_YAML = joinpath(NATRAJ_DIR, "natraj.yaml")
const NATRAJ_TAU  = 0.5
const NATRAJ_λ_NM = 360.0
const NATRAJ_μ    = [0.02, 0.06, 0.10, 0.16, 0.20, 0.28, 0.32, 0.40,
                     0.52, 0.64, 0.72, 0.84, 0.92, 0.96, 0.98, 1.00]
const NATRAJ_ϕs   = collect(0.0:30.0:180.0)

function run_natraj(FT::Type)
    I_modeled = zeros(Float64, length(NATRAJ_μ), length(NATRAJ_ϕs))
    Q_modeled = zeros(Float64, length(NATRAJ_μ), length(NATRAJ_ϕs))
    U_modeled = zeros(Float64, length(NATRAJ_μ), length(NATRAJ_ϕs))
    ν0 = 1e7 / NATRAJ_λ_NM
    for (iϕ, ϕ) in enumerate(NATRAJ_ϕs)
        params = parameters_from_yaml(NATRAJ_YAML)
        params.architecture = vSmartMOM.Architectures.CPU()
        params.float_type = FT
        params.spec_bands = [Float64[ν0, ν0 + 1.0]]
        params.vza = acosd.(NATRAJ_μ)
        params.sza = acosd(0.2)
        params.vaz = repeat([ϕ], length(NATRAJ_μ))
        model = model_from_parameters(params)
        model.τ_rayl[1] .= NATRAJ_TAU
        # vSmartMOM returns radiance factor L = I/F₀; Natraj table is reflectance
        # R = π · L (matches test/test_CoreRT.jl line 71).
        R = π .* Array(CoreRT.rt_run(model, i_band = 1)[1])
        I_modeled[:, iϕ] = R[:, 1, 1]
        Q_modeled[:, iϕ] = R[:, 2, 1]
        U_modeled[:, iϕ] = R[:, 3, 1]
    end
    by_stokes = Dict{Symbol, Vector{Row}}(:I => [], :Q => [], :U => [])
    for (iμ, μ) in enumerate(NATRAJ_μ), (iϕ, ϕ) in enumerate(NATRAJ_ϕs)
        push!(by_stokes[:I], Row("μ=$μ, ϕ=$(round(Int, ϕ))°",
                                  I_modeled[iμ, iϕ], Float64(I_trues[iμ, iϕ])))
        push!(by_stokes[:Q], Row("μ=$μ, ϕ=$(round(Int, ϕ))°",
                                  Q_modeled[iμ, iϕ], Float64(Q_trues[iμ, iϕ])))
        push!(by_stokes[:U], Row("μ=$μ, ϕ=$(round(Int, ϕ))°",
                                  U_modeled[iμ, iϕ], Float64(U_trues[iμ, iϕ])))
    end
    return by_stokes
end


# =============================================================================
# Case B — solar_tester scalar Stokes-I, Task 1   (🧪 VLIDORT regression)
# =============================================================================

include(joinpath(VLIDORT_REF, "solar_tester_atmosphere.jl"))
include(joinpath(VLIDORT_REF, "solar_tester_truth.jl"))

const ST_YAML        = joinpath(VLIDORT_CONFIG, "solar_tester.yaml")
const ST_AERO_G      = 0.8
const ST_AERO_OMEGA  = 0.95
const ST_AERO_TAU    = 0.5
const ST_AERO_NMOM   = 15

function st_aerosol_optics()
    Lp1 = ST_AERO_NMOM + 1
    α = zeros(Float64, Lp1); β = zeros(Float64, Lp1)
    γ = zeros(Float64, Lp1); δ = zeros(Float64, Lp1)
    ϵ = zeros(Float64, Lp1); ζ = zeros(Float64, Lp1)
    β[1] = 1.0
    for L in 1:ST_AERO_NMOM
        β[L + 1] = (2L + 1) * ST_AERO_G ^ L
    end
    return Scattering.AerosolOptics(
        greek_coefs = Scattering.GreekCoefs(α, β, γ, δ, ϵ, ζ),
        ω̃ = ST_AERO_OMEGA, k = 1.0, fᵗ = 0.0)
end

function st_aerosol_extinction()
    h = vcat(60.0, SOLAR_TESTER_HEIGHT_KM)
    n6 = 23 - 6
    parcel = ST_AERO_TAU / (h[n6 + 1] - h[end])
    aerext = zeros(Float64, 23)
    for n in (n6 + 1):23
        aerext[n] = parcel * (h[n] - h[n + 1])
    end
    return aerext
end

function inject_st!(model)
    aerext = st_aerosol_extinction()
    model.aerosol_optics[1][1] = st_aerosol_optics()
    for n in 1:23
        ext = SOLAR_TESTER_MOLEXT[n]
        ssa = SOLAR_TESTER_MOLOMG[n]
        model.τ_rayl[1][:, n] .= ssa * ext
        model.τ_abs[1][:, n]  .= (1 - ssa) * ext
        model.τ_aer[1][1, n]   = aerext[n]
    end
    return model
end

function st_run_at(sza, raz, FT::Type, yaml_path)
    params = parameters_from_yaml(yaml_path)
    params.architecture = vSmartMOM.Architectures.CPU()
    params.float_type = FT
    params.sza = float(sza)
    fill!(params.vaz, float(raz))
    model = model_from_parameters(params)
    inject_st!(model)
    out = CoreRT.rt_run(model, i_band = 1)
    return Array(out[1]), Array(out[2])
end

st_geom_index(i_sza, i_vza, i_raz) = (i_sza - 1) * 9 + (i_vza - 1) * 3 + i_raz

function run_solar_tester_scalar(FT::Type)
    sza   = 35.0
    raz   = SOLAR_TESTER_RAZ_DEG[1]    # 0°
    i_sza = 1; i_raz = 1
    task_idx  = SOLAR_TESTER_TASKS.pp_noDM
    toa_level = 1
    boa_level = length(SOLAR_TESTER_TAU_LEVELS)
    R, T = st_run_at(sza, raz, FT, ST_YAML)

    rows_up = Row[]
    rows_dn = Row[]
    for (i_vza, vza) in enumerate(SOLAR_TESTER_VZA_DEG)
        geom = st_geom_index(i_sza, i_vza, i_raz)
        truth_up = Float64(SOLAR_TESTER_STOKES[geom, toa_level, SOLAR_TESTER_DIR_UP, task_idx])
        truth_dn = Float64(SOLAR_TESTER_STOKES[geom, boa_level, SOLAR_TESTER_DIR_DN, task_idx])
        push!(rows_up, Row("vza=$(vza)°, raz=$(round(Int, raz))°, TOA-up",
                           Float64(R[i_vza, 1, 1]), truth_up))
        push!(rows_dn, Row("vza=$(vza)°, raz=$(round(Int, raz))°, BOA-dn",
                           Float64(T[i_vza, 1, 1]), truth_dn))
    end
    return Dict{Symbol, Vector{Row}}(:I_TOA_up => rows_up, :I_BOA_dn => rows_dn)
end

# =============================================================================
# Case C — solar_tester vector Stokes-IQU, Task 1   (🧪 VLIDORT regression)
# =============================================================================

include(joinpath(VLIDORT_REF, "solar_tester_problemIII_aerosol.jl"))
include(joinpath(VLIDORT_REF, "solar_tester_vector_truth.jl"))

const ST_VEC_YAML = joinpath(VLIDORT_CONFIG, "solar_tester_vector.yaml")

# Reuse st_aerosol_extinction() for layer placement; aerosol Greek coefs come
# from solar_tester_problemIII_aerosol.jl (Problem III).
const ST_VEC_NMOMENTS = 15
function st_vec_aerosol_optics()
    Lp1 = ST_VEC_NMOMENTS + 1
    α = collect(Float64, PROBLEMIII_a2[1:Lp1])
    β = collect(Float64, PROBLEMIII_a1[1:Lp1])
    γ = collect(Float64, PROBLEMIII_b1[1:Lp1])
    δ = collect(Float64, PROBLEMIII_a4[1:Lp1])
    ϵ = .-collect(Float64, PROBLEMIII_b2[1:Lp1])
    ζ = collect(Float64, PROBLEMIII_a3[1:Lp1])
    return Scattering.AerosolOptics(
        greek_coefs = Scattering.GreekCoefs(α, β, γ, δ, ϵ, ζ),
        ω̃ = 0.99999, k = 1.0, fᵗ = 0.0)
end

function inject_st_vec!(model)
    aerext = st_aerosol_extinction()
    model.aerosol_optics[1][1] = st_vec_aerosol_optics()
    for n in 1:23
        ext = SOLAR_TESTER_MOLEXT[n]
        ssa = SOLAR_TESTER_MOLOMG[n]
        model.τ_rayl[1][:, n] .= ssa * ext
        model.τ_abs[1][:, n]  .= (1 - ssa) * ext
        model.τ_aer[1][1, n]   = aerext[n]
    end
    return model
end

function run_solar_tester_vector(FT::Type)
    sza   = 35.0
    raz   = SOLAR_TESTER_VECTOR_RAZ_DEG[1]
    i_sza = 1; i_raz = 1
    task_idx  = 1                         # pp_noDM
    toa_level = 1
    boa_level = length(SOLAR_TESTER_VECTOR_TAU_LEVELS)

    params = parameters_from_yaml(ST_VEC_YAML)
    params.architecture = vSmartMOM.Architectures.CPU()
    params.float_type = FT
    params.sza = float(sza)
    fill!(params.vaz, float(raz))
    model = model_from_parameters(params)
    inject_st_vec!(model)
    out = CoreRT.rt_run(model, i_band = 1)
    R, T = Array(out[1]), Array(out[2])

    out_dict = Dict{Symbol, Vector{Row}}()
    for s_sym in (:I, :Q, :U)
        s_idx = _STOKES_IDX[s_sym]
        rows_up = Row[]; rows_dn = Row[]
        for (i_vza, vza) in enumerate(SOLAR_TESTER_VECTOR_VZA_DEG)
            geom = st_geom_index(i_sza, i_vza, i_raz)
            truth_up = Float64(SOLAR_TESTER_VECTOR_STOKES[s_sym][geom, toa_level, 1, task_idx])
            truth_dn = Float64(SOLAR_TESTER_VECTOR_STOKES[s_sym][geom, boa_level, 2, task_idx])
            s_sym in (:Q, :U) && (truth_up = -truth_up; truth_dn = -truth_dn)
            push!(rows_up, Row("vza=$(vza)°, raz=$(round(Int, raz))°, TOA-up",
                               Float64(R[i_vza, s_idx, 1]), truth_up))
            push!(rows_dn, Row("vza=$(vza)°, raz=$(round(Int, raz))°, BOA-dn",
                               Float64(T[i_vza, s_idx, 1]), truth_dn))
        end
        out_dict[Symbol(s_sym, :_TOA_up)] = rows_up
        out_dict[Symbol(s_sym, :_BOA_dn)] = rows_dn
    end
    return out_dict
end

# =============================================================================
# Page assembly
# =============================================================================

@info "Running benchmark cases (CPU, F64 + F32) — this takes a few minutes…"

results = Dict{Tuple{Symbol, Type}, Dict{Symbol, Vector{Row}}}()
for (case_sym, runner) in (
    (:siewert,        run_siewert),
    (:natraj,         run_natraj),
    (:st_scalar,      run_solar_tester_scalar),
    (:st_vector,      run_solar_tester_vector),
), FT in (Float64, Float32)
    @info "  → $(case_sym) [$(FT)]"
    results[(case_sym, FT)] = runner(FT)
end

mkpath(dirname(OUT_PATH))
open(OUT_PATH, "w") do io
    println(io, "# Benchmarks")
    println(io)
    println(io, "vSmartMOM is validated against two distinct kinds of references with")
    println(io, "very different epistemic standing. This page is **regenerated by**")
    println(io, "`docs/build_benchmarks.jl` whenever reference data, model parameters,")
    println(io, "or the RT solver itself change. Last regenerated: ",
            Base.Libc.strftime("%Y-%m-%d", time()),
            ", on Float64 + Float32 CPU.")
    println(io)
    println(io, "## Two kinds of references")
    println(io)
    println(io, "- **📚 Published benchmarks** — Natraj 2009 polarised Rayleigh tables")
    println(io, "  and Siewert 2000 PROBLEM IIA — are theoretically-derived, high-")
    println(io, "  accuracy references. The numerical procedures used to compute them")
    println(io, "  are documented in the original papers and have themselves been")
    println(io, "  independently cross-validated. **Treat any disagreement as a")
    println(io, "  vSmartMOM defect** to track down.")
    println(io, "- **🧪 VLIDORT 2.8.3 regression tables** — solar_tester scalar (Case B,")
    println(io, "  Stokes-I only) and solar_tester vector (Case C, Stokes-IQU) — are")
    println(io, "  *saved_results* tables shipped with the VLIDORT 2.8.3 distribution.")
    println(io, "  These compare two independently-implemented RT codes against each")
    println(io, "  other, **not** against an analytic ground truth. Disagreement at")
    println(io, "  the level of 10⁻³–10⁻² could be either code's numerical behavior at")
    println(io, "  the chosen NSTREAMS / convergence settings; we report it for")
    println(io, "  transparency, not as a correctness benchmark.")
    println(io)
    println(io, "## Error metric")
    println(io)
    println(io, "Each row reports the modeled value, the truth value, the absolute")
    println(io, "error `|Δ| = |modeled − truth|`, and the relative error")
    println(io, "`|Δ|/|truth|·100%`. Cells where `|truth|` falls below a per-case")
    println(io, "*near-zero threshold* are reported as `—` for the relative error")
    println(io, "and excluded from the median/max in the summary line — at small")
    println(io, "|truth| the relative error is dominated by the absolute noise floor")
    println(io, "and is not a meaningful accuracy measure. The threshold and the")
    println(io, "number of excluded cells are stated in each summary line.")
    println(io)
    println(io, "Where applicable, thresholds match the in-tree regression tests:")
    println(io, "for Natraj Q/U we use `|truth| < 0.01`, the same filter applied in")
    println(io, "[`test/test_CoreRT.jl`](../../../test/test_CoreRT.jl) (`@test maximum(Q_deltas[Q_modeled .≥ 0.01]) < 0.008`).")
    println(io)

    # Explicit per-case (stokes_label => near_zero_atol) ordering.
    # Dict iteration is unstable; using a Vector of Pairs preserves order.
    # near_zero_atol = 0 means no filtering. Natraj Q/U use 0.01 to match the
    # in-tree test_CoreRT.jl filter (Q_modeled .≥ 0.01); the relative error at
    # smaller |truth| is dominated by the absolute noise floor.
    case_meta = (
        (:siewert, :published, "Case A — Siewert 2000 PROBLEM IIA",
         """A single-layer aerosol slab (`τ = 1.0`, `ω̃ = 0.973527`) with the
         PROBLEM IIA Greek phase-matrix coefficients shipped with VLIDORT
         2.8.3. TOA upwelling Stokes vector at 11 viewing-zenith cosines and
         3 azimuths (0°, 90°, 180°). Q/U/V truth signs flipped to vSmartMOM's
         Hovenier convention.""",
         [:I => 0.0, :Q => 0.0, :U => 0.0, :V => 0.0]),
        (:natraj, :published,
         "Natraj 2009 — polarised Rayleigh τ=0.5",
         """Pure Rayleigh, optical depth `τ = 0.5`, sza = `acos(0.2)` ≈ 78.46°,
         depolarization `0` (Natraj's idealization). 16 viewing-zenith cosines
         from 0.02 to 1.00 across 7 azimuths from 0° to 180°. Modeled
         reflectance `R = π · (I/F₀)`. **`GaussLegQuad` + `NoTruncation()`**:
         half-space Gauss-Legendre on `[0, 1]` (5–50× more accurate than
         RadauQuad on this Rayleigh-only setup) and Rayleigh has only
         `β₀, β₁, β₂` non-zero so δ-fit forward-peak truncation is
         meaningless. Q/U use a near-zero filter of `|truth| < 0.01`
         (mirrors the `|modeled| ≥ 0.01` filter applied in
         [`test/test_CoreRT.jl`](../../../test/test_CoreRT.jl)).""",
         [:I => 0.0, :Q => 0.01, :U => 0.01]),
        (:st_scalar, :vlidort_regression, "Case B — solar_tester scalar (Task 1)",
         """23-layer atmosphere from VLIDORT's `2p8p3_solar_tester.f90`,
         Stokes-I only, isotropic Henyey-Greenstein aerosol (`g=0.8`,
         `ω̃=0.95`, `τ=0.5`) with NMOMENTS=15 to match VLIDORT's
         NSTREAMS=8. SZA=35°, RAZ=0°, 3 viewing zenith angles. TOA upwelling
         and BOA downwelling Stokes-I.""",
         [:I_TOA_up => 0.0, :I_BOA_dn => 0.0]),
        (:st_vector, :vlidort_regression, "Case C — solar_tester vector Stokes-IQU (Task 1)",
         """Same 23-layer atmosphere as Case B but with vector RT
         (Stokes-IQU) and the PROBLEM-III aerosol Greek coefficients shipped
         with VLIDORT. Q/U truth signs flipped to vSmartMOM's Hovenier
         convention.""",
         [:I_TOA_up => 0.0, :I_BOA_dn => 0.0,
          :Q_TOA_up => 0.0, :Q_BOA_dn => 0.0,
          :U_TOA_up => 0.0, :U_BOA_dn => 0.0]),
    )
    provenance_label(p::Symbol) = p === :published ? "📚 published benchmark" : "🧪 VLIDORT regression table"

    for (case_sym, prov, title, descr, stokes_order) in case_meta
        println(io, "## ", title)
        println(io)
        println(io, "*", provenance_label(prov), "*")
        println(io)
        println(io, descr)
        println(io)
        for FT in (Float64, Float32)
            println(io, "### Float", FT === Float64 ? "64" : "32")
            println(io)
            by_stokes = results[(case_sym, FT)]
            for (s_label, atol) in stokes_order
                haskey(by_stokes, s_label) || continue
                rows = by_stokes[s_label]
                println(io, "**Stokes ", s_label, "**  ",
                        summary_line(rows; near_zero_atol = atol))
                println(io)
                write_table(io, rows; near_zero_atol = atol)
                println(io)
            end
        end
    end

    println(io, "## Reproducing")
    println(io)
    println(io, "```bash")
    println(io, "julia --project=docs docs/build_benchmarks.jl")
    println(io, "```")
    println(io)
    println(io, "Reference data lives under [`test/vlidort_baseline/reference_data/`](../../../test/vlidort_baseline/reference_data/)")
    println(io, "(Siewert + VLIDORT-shipped) and [`test/benchmarks/natraj_trues.jl`](../../../test/benchmarks/natraj_trues.jl)")
    println(io, "(Natraj 2009 transcribed). Configurations are at")
    println(io, "[`test/vlidort_baseline/configs/`](../../../test/vlidort_baseline/configs/)")
    println(io, "and [`test/benchmarks/natraj.yaml`](../../../test/benchmarks/natraj.yaml).")
end

@info "✅ Benchmarks page written" path = OUT_PATH
