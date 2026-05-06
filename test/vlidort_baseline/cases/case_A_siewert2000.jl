using Test
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering

# Loads SIEWERT_α/β/γ/δ/ϵ/ζ (ϵ already sign-flipped per VLIDORT mapping)
# and SIEWERT_TRUTH (Dict (az_deg, stokes_sym) → 22×7 matrix).
include(joinpath(REF_DIR, "siewert2000_IIA_greek.jl"))
include(joinpath(REF_DIR, "siewert2000_IIA_truth.jl"))

const SIEWERT_YAML = joinpath(CONFIG_DIR, "siewert2000_IIA.yaml")
const SIEWERT_τ_TOTAL = 1.0
const SIEWERT_ω̃     = 0.973527

# Cosines used by VLIDORT tables (column 1 of each table). All eight tables
# share the same cosine grid (-1.0, -0.9, …, -0.0, 0.0, …, 1.0).
const SIEWERT_TABLE_COSINES = SIEWERT_TABLE2_COSINES

"Hand-built `AerosolOptics` carrying the Siewert PROBLEM_IIA Greek coefficients."
function siewert_aerosol_optics()
    greek = Scattering.GreekCoefs(SIEWERT_α, SIEWERT_β, SIEWERT_γ,
                                  SIEWERT_δ, SIEWERT_ϵ, SIEWERT_ζ)
    return Scattering.AerosolOptics(greek_coefs=greek, ω̃=SIEWERT_ω̃,
                                    k=1.0, fᵗ=0.0)
end

"Run vSmartMOM at all 11 viewing angles for one azimuth. Returns L = I/F₀
shape `(11, 4, 1)`."
function siewert_run_at(az_deg::Real)
    params = parameters_from_yaml(SIEWERT_YAML)
    fill!(params.vaz, float(az_deg))
    model = model_from_parameters(params)
    model.aerosol_optics[1][1] = siewert_aerosol_optics()
    model.τ_aer[1][1, :] .= SIEWERT_τ_TOTAL
    model.τ_rayl[1] .= 0.0
    return CoreRT.rt_run(model, i_band=1)[1]
end

"For each VZA cosine, find the table row at -|μ| (TOA upwelling) and read
column 1 (τ-level=0)."
function pick_toa_upwelling(table::AbstractMatrix, vza_cosines::AbstractVector)
    out = similar(vza_cosines, Float64)
    for (i, μ) in enumerate(vza_cosines)
        target = -abs(μ)
        idx = argmin(abs.(SIEWERT_TABLE_COSINES .- target))
        out[i] = table[idx, 1]
    end
    return out
end

@testset "VLIDORT baseline: Case A — Siewert 2000 Problem IIA" begin
    vza_cosines = [cosd(0.0001), cosd(25.841932763), cosd(36.869897646),
                   cosd(45.572995999), cosd(53.130102354), cosd(60.0),
                   cosd(66.421821522), cosd(72.542396876), cosd(78.463040967),
                   cosd(84.260829523), cosd(89.9999)]

    # Per-(az, Stokes) results: max & median rel error over 11 viewing angles
    # at TOA upwelling. Excludes near-zero entries (|truth| < 1e-9).
    results = NamedTuple{(:az, :stokes, :max_rel, :median_rel),
                         Tuple{Float64,Symbol,Float64,Float64}}[]

    stokes_idx = Dict(:I => 1, :Q => 2, :U => 3, :V => 4)

    for az_deg in (0.0, 90.0, 180.0)
        @info "running Siewert IIA at azimuth $az_deg"
        L = siewert_run_at(az_deg)
        for stokes_sym in (:I, :Q, :U, :V)
            haskey(SIEWERT_TRUTH, (az_deg, stokes_sym)) || continue
            truth_table = SIEWERT_TRUTH[(az_deg, stokes_sym)]
            modeled = π .* L[:, stokes_idx[stokes_sym], 1]
            truth = pick_toa_upwelling(truth_table, vza_cosines)
            re = [abs(truth[i]) < 1e-9 ? NaN :
                  abs(modeled[i] - truth[i]) / abs(truth[i])
                  for i in eachindex(modeled)]
            valid = filter(!isnan, re)
            isempty(valid) && continue
            push!(results, (az=az_deg, stokes=stokes_sym,
                            max_rel=maximum(valid),
                            median_rel=Statistics.median(valid)))
        end
    end

    println("\nCase A — Siewert 2000 IIA, TOA-upwelling rel-error summary:")
    println("    az    stokes    median       max")
    for r in results
        println(string("    ", lpad(round(r.az, digits=1), 5), "    ",
                       lpad(string(r.stokes), 4), "    ",
                       lpad(string(round(r.median_rel, sigdigits=3)), 8), "    ",
                       lpad(string(round(r.max_rel, sigdigits=3)), 8)))
    end

    # Per-Stokes tolerance gates. With μ₀ = 0.6 (matching the f90 driver
    # override of the .cfg SZA), observed maxima are ~5e-6 for I/U/V and
    # ~3e-3 for Q (Q max-rel hits a near-zero crossing). Gates set ~10×
    # observed to flag regressions while tolerating numerical-noise drift.
    for r in results
        rtol = r.stokes == :I ? 1e-4 :
               r.stokes == :Q ? 1e-2 :
               r.stokes == :U ? 1e-4 : 1e-4
        @test r.max_rel < rtol
    end
end
