# =============================================================================
# pilot_lut_resolution.jl — Line-by-line Δν convergence study
#
# Goal: decide the wavenumber grid spacing for the rebuilt HITRAN LUTs.
# For two representative bands (O2 A-band and H2O 1.94 μm), compute direct
# line-by-line cross-sections at several Δν values and report how the
# band-integrated absorption converges toward a fine reference.
#
# One (p, T) — 500 hPa, 250 K — that is a reasonable column-mean, though
# low-pressure stratospheric lines will be narrower and need finer Δν.
#
# Output: printed table of {band, Δν, n_points, σ_mean, Tmean, max_error}
#         and a compact JLD2 file with spectra for follow-on plotting.
# =============================================================================

using vSmartMOM
using vSmartMOM.Absorption
using vSmartMOM.Architectures
using Printf
using Statistics
using JLD2

const OUTFILE = expanduser("~/pilot_lut_resolution_results.jld2")
const ARCH    = Architectures.GPU()   # switch to CPU() if no GPU
const CEF     = HumlicekWeidemann32SDErrorFunction()
const WING_CUTOFF = 40.0
# Test at two (p, T) conditions: tropospheric and stratospheric. Line widths
# shrink at low pressure, so stratospheric conditions are the stressful case
# for Δν selection.
const PT_CASES = [
    (label = "trop_500hPa_250K",  p = 500.0, T = 250.0),
    (label = "strat_10hPa_220K",  p = 10.0,  T = 220.0),
]

# Reference full-atmosphere column densities (molec / cm²). Used to turn σ
# into a band-mean transmittance exp(-σ·N) that is physically meaningful.
const N_O2   = 4.5e24
const N_H2O  = 5.0e22

# Bands to scan (in cm⁻¹).
# `name` is the HITRAN molecule name used by `fetch_hitran`; `mol`/`iso` are
# the mol/iso IDs used by `read_hitran` to filter the .par file.
const BANDS = Dict(
    "O2_Aband"   => (name = "O2",  mol = 7, iso = 1, ν_lo = 13000.0, ν_hi = 13200.0, N_col = N_O2),
    "H2O_1p94um" => (name = "H2O", mol = 1, iso = 1, ν_lo = 5050.0,  ν_hi = 5250.0,  N_col = N_H2O),
)

# Δν values to test (cm⁻¹). First value is the reference.
const DELTA_NU = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]

"""
    compute_band_spectrum(hitran_table, ν_grid; p, T, arch)

Return the cross-section σ(ν) [cm²/molec] on `ν_grid` at (p, T).
"""
function compute_band_spectrum(hitran_tbl, ν_grid; p, T, arch = ARCH)
    model = make_hitran_model(hitran_tbl, Voigt();
                              wing_cutoff = WING_CUTOFF,
                              CEF         = CEF,
                              architecture = arch)
    return collect(compute_absorption_cross_section(model, collect(ν_grid), p, T))
end

"""
    band_metrics(σ_test, σ_ref, ν_test, ν_ref; N_col)

Compare `σ_test` on `ν_test` to `σ_ref` on `ν_ref` (finer grid).

Reports:
  σ_mean_test, σ_mean_ref
  Tmean_test, Tmean_ref — mean in-band transmittance exp(-σ·N_col)
  |ΔTmean|              — band-mean transmittance error
  max_point_T, mean_point_T — pointwise transmittance errors
"""
function band_metrics(σ_test, σ_ref, ν_test, ν_ref; N_col)
    σ_mean_test = mean(σ_test)
    σ_mean_ref  = mean(σ_ref)

    T_test = exp.(-σ_test .* N_col)
    T_ref  = exp.(-σ_ref  .* N_col)

    Tmean_test = mean(T_test)
    Tmean_ref  = mean(T_ref)

    # For pointwise comparison, interpolate σ_ref linearly onto ν_test.
    σ_ref_on_test = similar(σ_test)
    @inbounds for (i, ν) in enumerate(ν_test)
        j = searchsortedfirst(ν_ref, ν)
        if j <= 1
            σ_ref_on_test[i] = σ_ref[1]
        elseif j > length(ν_ref)
            σ_ref_on_test[i] = σ_ref[end]
        else
            ν1, ν2 = ν_ref[j - 1], ν_ref[j]
            s1, s2 = σ_ref[j - 1], σ_ref[j]
            w = (ν - ν1) / (ν2 - ν1)
            σ_ref_on_test[i] = (1 - w) * s1 + w * s2
        end
    end
    T_ref_on_test = exp.(-σ_ref_on_test .* N_col)

    max_point_T  = maximum(abs.(T_test .- T_ref_on_test))
    mean_point_T = mean(abs.(T_test   .- T_ref_on_test))
    dTmean       = abs(Tmean_test - Tmean_ref)

    return (; σ_mean_test, σ_mean_ref,
              Tmean_test, Tmean_ref, dTmean,
              max_point_T, mean_point_T)
end

function main()
    results = Dict{String, Any}()

    for (band_name, b) in BANDS
        println("\n" * "="^72)
        println("Band: $band_name   [ν=$(b.ν_lo) … $(b.ν_hi) cm⁻¹,  N_col=$(b.N_col)]")
        println("="^72)

        # Fetch & parse HITRAN once per band, with a buffer on each side for
        # wing_cutoff (default 40 cm⁻¹). `fetch_hitran` pulls all isotopologues
        # of the named molecule; `read_hitran` then filters to mol/iso.
        numin, numax = b.ν_lo - WING_CUTOFF, b.ν_hi + WING_CUTOFF
        @info "Fetching HITRAN" molecule = b.name numin numax
        path = fetch_hitran(b.name; numin = numin, numax = numax)
        tbl  = read_hitran(path; mol = b.mol, iso = b.iso,
                           ν_min = numin, ν_max = numax)
        @info "  lines retained" n_lines = length(tbl.νᵢ)

        results[band_name] = Dict{String, Any}()

        for pt in PT_CASES
            println("\n--- (p, T) = ($(pt.p) hPa, $(pt.T) K) ---")

            # Reference (finest grid).
            Δν_ref = DELTA_NU[1]
            ν_ref  = b.ν_lo:Δν_ref:b.ν_hi
            t0 = time()
            σ_ref = compute_band_spectrum(tbl, ν_ref; p = pt.p, T = pt.T)
            t_ref = time() - t0
            @info "  reference σ computed" Δν_ref n_points = length(ν_ref) t_ref

            pt_result = Dict{String, Any}(
                "nu_ref"    => collect(ν_ref),
                "sigma_ref" => σ_ref,
                "dnu_ref"   => Δν_ref,
                "p_hPa"     => pt.p,
                "T_K"       => pt.T,
            )

            @printf("\n %-8s  %-10s  %-12s  %-12s  %-12s  %-12s  %-8s\n",
                    "Δν", "n_pts", "σ_mean(test)", "Tmean(test)", "Tmean(ref)", "|ΔTmean|", "t (s)")
            println("-"^92)

            for Δν in DELTA_NU
                ν_test = b.ν_lo:Δν:b.ν_hi
                t0 = time()
                σ_test = compute_band_spectrum(tbl, ν_test; p = pt.p, T = pt.T)
                t_run  = time() - t0

                m = band_metrics(σ_test, σ_ref, collect(ν_test), collect(ν_ref);
                                 N_col = b.N_col)
                @printf(" %-8.4f  %-10d  %-12.4e  %-12.6f  %-12.6f  %-12.3e  %-8.2f\n",
                        Δν, length(ν_test), m.σ_mean_test, m.Tmean_test,
                        m.Tmean_ref, m.dTmean, t_run)

                key = @sprintf("dnu_%.4f", Δν)
                pt_result[key] = Dict(
                    "dnu"          => Δν,
                    "n_points"     => length(ν_test),
                    "sigma_mean"   => m.σ_mean_test,
                    "Tmean_test"   => m.Tmean_test,
                    "Tmean_ref"    => m.Tmean_ref,
                    "dTmean"       => m.dTmean,
                    "max_point_T"  => m.max_point_T,
                    "mean_point_T" => m.mean_point_T,
                    "t_run"        => t_run,
                )
            end

            results[band_name][pt.label] = pt_result
        end
    end

    @info "Saving pilot results" OUTFILE
    @save OUTFILE results
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
