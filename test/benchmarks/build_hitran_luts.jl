# =============================================================================
# build_hitran_luts.jl — Rebuild the HITRAN absorption LUTs
#
# Produces one JLD2 per species under `~/data/HITRAN_LUTs/` (override via
# `HITRAN_LUT_DIR`). Spec per plans/HANDOFF_LUT_REBUILD_2026-04-23.md:
#
#   Δν         = 0.01 cm⁻¹  (pilot confirmed Δν≤0.01 caps ⟨ΔT⟩ ≤ 3.3e-5)
#   ν range    = 2857.14–35714.29 cm⁻¹ (280–3500 nm)
#   p grid     = 15 log-spaced points from 0.01 to 1080 hPa
#   T grid     = 10 points 180–360 K
#   storage    = Float32
#   interp     = BSpline(Linear())  (no cubic over/undershoot)
#   broadening = Voigt, wing_cutoff = 40 cm⁻¹, HW32SD CEF
#   arch       = GPU by default (override via `HITRAN_LUT_ARCH=cpu`)
#
# Species list: 7 already-rebuilt + 15 missing that MODTRAN activates by
# default in the 280–3500 nm window (see plan doc for rationale).
#
# ENV overrides:
#   HITRAN_LUT_DIR       default ~/data/HITRAN_LUTs
#   HITRAN_LUT_SPECIES   comma-separated subset, e.g. "NO2,SO2"; default = all
#   HITRAN_LUT_ARCH      "gpu" (default) or "cpu"
#   HITRAN_LUT_FORCE     "1" to overwrite existing output files
# =============================================================================

using vSmartMOM
using vSmartMOM.Absorption
using vSmartMOM.Architectures
using Printf
using Dates
using JLD2

# ------------------------------------------------------------------
# Build spec
# ------------------------------------------------------------------
const Δν_CM     = 0.01
const ν_LO_CM   = 1e7 / 3500.0           # ≈ 2857.14
const ν_HI_CM   = 1e7 /  280.0           # ≈ 35714.29
const ν_GRID    = ν_LO_CM:Δν_CM:ν_HI_CM  # ~3.286M points

# Pressure & temperature grids.
#
# `scale(model.itp, ν_grid, p_grid, t_grid)` (the runtime path in
# `compute_absorption_cross_section.jl`) requires AbstractRange with uniform
# spacing, so we use linear spacing here. 15 pts × 10 pts is what was agreed
# in the plan. Step = 77.14 hPa / 20 K.
const P_GRID = range(0.01, 1080.01; length = 15)   # hPa, Δp ≈ 77.14
const T_GRID = range(180.0, 360.0;  length = 10)   # K,   ΔT = 20

const WING_CUTOFF = 40.0
const CEF         = HumlicekWeidemann32SDErrorFunction()
const BROADENING  = Voigt()

const DEFAULT_LUT_DIR = expanduser("~/data/HITRAN_LUTs")

# ------------------------------------------------------------------
# Species list.
#   - `name`: HITRAN molecule name for `fetch_hitran`
#   - `mol`:  HITRAN molecule id for `read_hitran` filter
#   - `isos`: isotopologues to include; use `:all` (every iso in the .par)
#             or a vector of local iso indices. For most radiatively
#             important transitions iso=1 dominates; :all is safer.
# ------------------------------------------------------------------
const SPECIES = (
    # --- existing 7 (rebuilt in the new format) ---
    (name = "H2O", mol =  1, isos = :all),
    # CO2: iso=1 only. HITRAN CO2 has >9 isotopologues whose iso-column uses
    # letter codes ("A"=10, "B"=11, ...); `read_hitran` only parses digits,
    # so those rows collapse to iso=0 and then crash `mol_weight`. iso=1
    # (12C16O2) is ~98% of atmospheric CO2 anyway.
    (name = "CO2", mol =  2, isos = [1]),
    (name = "O3",  mol =  3, isos = :all),
    (name = "N2O", mol =  4, isos = :all),
    (name = "CO",  mol =  5, isos = :all),
    (name = "CH4", mol =  6, isos = :all),
    (name = "O2",  mol =  7, isos = :all),
    # --- MODTRAN-equivalent additions ---
    (name = "NO",   mol =  8, isos = :all),
    (name = "SO2",  mol =  9, isos = :all),
    (name = "NO2",  mol = 10, isos = :all),
    (name = "NH3",  mol = 11, isos = :all),
    (name = "HNO3", mol = 12, isos = :all),
    (name = "HF",   mol = 14, isos = :all),
    (name = "HCl",  mol = 15, isos = :all),
    (name = "ClO",  mol = 18, isos = :all),
    (name = "OCS",  mol = 19, isos = :all),
    (name = "H2CO", mol = 20, isos = :all),
    (name = "HOCl", mol = 21, isos = :all),
    (name = "HCN",  mol = 23, isos = :all),
    (name = "H2O2", mol = 25, isos = :all),
    (name = "C2H2", mol = 26, isos = :all),
    (name = "C2H6", mol = 27, isos = :all),
)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
function parse_arch(val::AbstractString)
    lowercase(val) in ("gpu", "cuda") && return Architectures.GPU()
    lowercase(val) == "cpu" && return Architectures.CPU()
    throw(ArgumentError("HITRAN_LUT_ARCH must be 'gpu' or 'cpu', got $val"))
end

function parse_species_filter(val::AbstractString, all_species)
    want = Set(strip.(split(val, ",")))
    isempty(want) && return all_species
    return filter(s -> s.name in want, all_species)
end

"""
    load_hitran_table(species; ν_lo=ν_LO_CM, ν_hi=ν_HI_CM)

Download + parse HITRAN lines for `species`. Returns a concatenated
`HitranTable` that merges all requested isotopologues.
"""
function load_hitran_table(species; ν_lo = ν_LO_CM, ν_hi = ν_HI_CM)
    # Extend by wing_cutoff on each side so Voigt wings that reach into the
    # grid are captured.
    numin = ν_lo - WING_CUTOFF
    numax = ν_hi + WING_CUTOFF
    @info "  fetch_hitran" molecule = species.name numin numax
    path = fetch_hitran(species.name; numin = numin, numax = numax)
    # iso = -1 → keep every iso present in the file; read_hitran still
    # honours the ν_min / ν_max bounds.
    iso_filter = species.isos === :all ? -1 : first(species.isos)  # read_hitran only supports single iso
    tbl = read_hitran(path; mol = species.mol, iso = iso_filter,
                      ν_min = numin, ν_max = numax)
    @info "  lines retained" mol = species.mol iso = iso_filter n_lines = length(tbl.νᵢ)
    return tbl
end

"""
    build_one_species!(species, out_dir; arch, force)

Build a single LUT. Skips if output exists unless `force` is true.
Writes `{out_dir}/{species.name}.jld2` with one top-level `itp_model`.
"""
function build_one_species!(species, out_dir; arch, force::Bool)
    path = joinpath(out_dir, species.name * ".jld2")
    if isfile(path) && !force
        @info "[skip] already built" path
        return path
    end

    @info "[build] $(species.name)  (mol=$(species.mol))"
    t0 = time()
    tbl = load_hitran_table(species)
    if length(tbl.νᵢ) == 0
        @warn "  zero lines in range — writing empty sentinel" species = species.name
    end
    @info "  building interpolation model" n_ν = length(ν_GRID) n_p = length(P_GRID) n_T = length(T_GRID) arch
    itp_model = make_interpolation_model(tbl, BROADENING, ν_GRID, P_GRID, T_GRID;
                                         wing_cutoff  = WING_CUTOFF,
                                         CEF          = CEF,
                                         architecture = arch,
                                         interp_type  = :linear,
                                         storage_type = Float32)
    save_interpolation_model(itp_model, path)
    t_tot = time() - t0
    @info "  saved" path size_MB = round(filesize(path) / 1024^2; digits = 1) t_tot
    return path
end

# ------------------------------------------------------------------
# Main driver
# ------------------------------------------------------------------
function main()
    out_dir = get(ENV, "HITRAN_LUT_DIR", DEFAULT_LUT_DIR)
    arch    = parse_arch(get(ENV, "HITRAN_LUT_ARCH", "gpu"))
    force   = get(ENV, "HITRAN_LUT_FORCE", "0") == "1"
    want    = get(ENV, "HITRAN_LUT_SPECIES", "")

    mkpath(out_dir)
    species_list = want == "" ? collect(SPECIES) : parse_species_filter(want, collect(SPECIES))

    @info "=== HITRAN LUT build ==="
    @info "  Δν            = $(Δν_CM) cm⁻¹"
    @info "  ν range       = $(round(ν_LO_CM; digits=2)) – $(round(ν_HI_CM; digits=2)) cm⁻¹ ($(length(ν_GRID)) points)"
    @info "  p grid        = $(length(P_GRID)) points: $(first(P_GRID))…$(last(P_GRID)) hPa"
    @info "  T grid        = $(length(T_GRID)) points: $(first(T_GRID))…$(last(T_GRID)) K"
    @info "  architecture  = $arch"
    @info "  out_dir       = $out_dir"
    @info "  # species     = $(length(species_list))"
    @info "  storage       = Float32, interp=linear"

    t_start = time()
    for (k, species) in enumerate(species_list)
        @info "--- species $k / $(length(species_list)): $(species.name) ---"
        try
            build_one_species!(species, out_dir; arch = arch, force = force)
        catch e
            @error "  FAILED on $(species.name)" exception = (e, catch_backtrace())
            # keep going with the remaining species
        end
    end
    @info "=== done in $(round(time() - t_start; digits = 1)) s ==="
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
