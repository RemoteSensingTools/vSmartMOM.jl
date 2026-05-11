## =============================================================================
## create_HITRAN_LUTs.jl
##
## One-time script to pre-compute broadband absorption cross-section LUTs for
## every HITRAN molecule used in the EMIT / MODTRAN comparison.  Each LUT is a
## 3-D (ν, p, T) B-spline interpolation table stored as a JLD2 file.  Once
## created, adding the paths to the YAML `LUTfiles` section replaces the ~30 min
## line-by-line HITRAN computation with a ~10 s LUT load.
##
## Usage:
##   julia --project=/home/sanghavi/code/github/vSmartMOM.jl \
##         /home/sanghavi/code/github/vSmartMOM.jl/test/benchmarks/create_HITRAN_LUTs.jl
##
## After completion, update ParamsEMIT_MODTRANcomp.yaml:
##
##   LUTfiles:
##     O2:  "/home/sanghavi/data/EMIT_MODTRANcomp/LUTs/O2_emit_lut.jld2" # for o2 a-band in OCO-2/3 simulations: "/net/fluo/data1/ABSCO_CS_Database/v5.2_final/o2_v52_v2.jld2"
##     CO2: "/home/sanghavi/data/EMIT_MODTRANcomp/LUTs/CO2_emit_lut.jld2"
##     H2O: "/home/sanghavi/data/EMIT_MODTRANcomp/LUTs/H2O_emit_lut.jld2"
##     CH4: "/home/sanghavi/data/EMIT_MODTRANcomp/LUTs/CH4_emit_lut.jld2"
##     O3:  "/home/sanghavi/data/EMIT_MODTRANcomp/LUTs/O3_emit_lut.jld2"
##     N2O: "/home/sanghavi/data/EMIT_MODTRANcomp/LUTs/N2O_emit_lut.jld2"
##     CO:  "/home/sanghavi/data/EMIT_MODTRANcomp/LUTs/CO_emit_lut.jld2"
## =============================================================================

using CUDA
CUDA.device!(1)

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.Absorption: read_hitran, make_interpolation_model,
                            save_interpolation_model, load_interpolation_model,
                            absorption_cross_section, make_hitran_model,
                            Voigt, HumlicekWeidemann32SDErrorFunction
using Dates

# ---- Configuration ----------------------------------------------------------

const LUT_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/LUTs"

# Spectral grid: covers the EMIT simulation range (380–2510 nm → 3984–26316 cm⁻¹)
# at 0.1 cm⁻¹ resolution, matching the simulation grid in ParamsEMIT_MODTRANcomp.yaml.
const ν_min  = Float64(1e7 / 2510.0)    # ≈ 3984 cm⁻¹
const ν_max  = Float64(1e7 / 380.0)     # ≈ 26316 cm⁻¹
const ν_step = 0.1
const ν_grid = ν_min:ν_step:ν_max

# Pressure grid: 0.01 to 1100 hPa in 20 hPa steps.
# Covers the AFGL midlatitude winter profile (3.6e-5 to 1018 hPa) with margin.
const p_grid = 0.01:20.0:1100.0

# Temperature grid: 150 to 330 K in 5 K steps.
# Covers the AFGL range (199–313 K) with margin for interpolation.
const t_grid = 150.0:5.0:330.0

# Molecules to process: (name, wing_cutoff)
# The q-based H₂O in model_from_parameters uses wing_cutoff=150;
# the YAML-based H₂O and all other molecules use wing_cutoff=40.
# We create TWO H₂O LUTs (wing_cutoff=150 and 40) for completeness.
const MOLECULES = [
    ("H2O",  150, "H2O_wc150_emit_lut.jld2"),  # for q-based H₂O
    ("H2O",   40, "H2O_emit_lut.jld2"),         # for YAML-based variable H₂O
    ("CO2",   40, "CO2_emit_lut.jld2"),
    ("CH4",   40, "CH4_emit_lut.jld2"),
    ("O3",    40, "O3_emit_lut.jld2"),
    ("N2O",   40, "N2O_emit_lut.jld2"),
    ("CO",    40, "CO_emit_lut.jld2"),
    ("O2",    40, "O2_emit_lut.jld2"),
]

# ---- Main -------------------------------------------------------------------

function create_all_luts()

    mkpath(LUT_DIR)

    @info "LUT creation parameters" ν_range=(ν_min, ν_max) ν_step n_ν=length(ν_grid) n_p=length(p_grid) n_T=length(t_grid) n_total_pts=length(ν_grid)*length(p_grid)*length(t_grid)

    for (mol, wc, fname) in MOLECULES

        out_path = joinpath(LUT_DIR, fname)

        if isfile(out_path)
            @info "LUT already exists, skipping" mol out_path
            continue
        end

        @info "=== Creating LUT for $mol (wing_cutoff=$wc) ==="
        t0 = time()

        # Load HITRAN line list
        hitran_data = read_hitran(artifact(mol), iso=1)
        @info "  HITRAN lines loaded" mol n_lines=length(hitran_data.νᵢ)

        # Create the 3-D interpolation model.
        # This calls compute_absorption_cross_section for every (p, T) pair
        # on the GPU, then fits a B-spline interpolator.
        lut = make_interpolation_model(
                  hitran_data, Voigt(),
                  ν_grid, p_grid, t_grid;
                  wing_cutoff   = wc,
                  CEF           = HumlicekWeidemann32SDErrorFunction(),
                  architecture  = vSmartMOM.Architectures.GPU(),
                  vmr           = 0)

        # Save to disk
        save_interpolation_model(lut, out_path)
        elapsed = (time() - t0) / 60.0
        @info "  LUT saved" mol path=out_path elapsed_min=round(elapsed; digits=1)
    end

    # ---- Validation ---------------------------------------------------------
    @info "=== Validating LUTs against HITRAN line-by-line ==="

    # Test at a mid-troposphere (p, T) point
    p_test = 500.0
    T_test = 250.0
    ν_test = collect(ν_grid)

    for (mol, wc, fname) in MOLECULES
        out_path = joinpath(LUT_DIR, fname)
        isfile(out_path) || continue

        lut = load_interpolation_model(out_path)

        hitran_data = read_hitran(artifact(mol), iso=1)
        hitran_model = make_hitran_model(hitran_data, Voigt();
                           wing_cutoff  = wc,
                           CEF          = HumlicekWeidemann32SDErrorFunction(),
                           architecture = vSmartMOM.Architectures.GPU(),
                           vmr          = 0)

        σ_lut    = absorption_cross_section(lut, ν_test, p_test, T_test)
        σ_hitran = Array(absorption_cross_section(hitran_model, ν_test, p_test, T_test))

        # Relative error where cross-section is non-negligible
        mask = abs.(σ_hitran) .> 1e-30
        if any(mask)
            rel_err = abs.(σ_lut[mask] .- σ_hitran[mask]) ./ abs.(σ_hitran[mask])
            @info "  $mol (wc=$wc): max_rel_err = $(round(maximum(rel_err); sigdigits=3)), " *
                  "mean_rel_err = $(round(mean(rel_err); sigdigits=3)), " *
                  "n_nonzero = $(count(mask))"
        else
            @info "  $mol (wc=$wc): no absorption lines in test range"
        end
    end

    @info "=== All done ==="
    @info "Add the following to ParamsEMIT_MODTRANcomp.yaml under absorption.LUTfiles:"
    for (mol, wc, fname) in MOLECULES
        wc == 150 && continue    # q-based H₂O LUT needs code change, not YAML
        println("    $mol: \"$(joinpath(LUT_DIR, fname))\"")
    end
    println()
    @info "For the q-based H₂O (wing_cutoff=150), use:" H2O_wc150=joinpath(LUT_DIR, "H2O_wc150_emit_lut.jld2")
end

using Statistics: mean
create_all_luts()
