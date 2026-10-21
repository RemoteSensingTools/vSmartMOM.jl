#!/usr/bin/env julia

# High-pressure numerical convergence of isolated-Voigt line cutoffs for the
# ExoEarth bands. H2O=25 cm^-1 is tested but interpreted separately because
# MT_CKD defines the complementary continuum beyond the local-line treatment.

using vSmartMOM
using vSmartMOM.Absorption
using vSmartMOM.Architectures
using Printf

const CASES = [
    ("O2_A",  "O2",  7, 12820.0, 13520.0),
    ("CO2_16", "CO2", 2, 5880.0, 6455.0),
    ("CH4_16", "CH4", 6, 5880.0, 6455.0),
    ("H2O_16", "H2O", 1, 5880.0, 6455.0),
    ("CO2_204", "CO2", 2, 4760.0, 5130.0),
    ("H2O_204", "H2O", 1, 4760.0, 5130.0),
    ("CH4_23", "CH4", 6, 4165.0, 4547.0),
    ("H2O_23", "H2O", 1, 4165.0, 4547.0),
]
const CUTOFFS = [25, 50, 100, 150, 250, 500, 1000]
const PRESSURES_HPA = [2000.0, 5000.0, 9000.0]
const TEMPERATURES_K = Dict(2000.0 => 228.863118,
                            5000.0 => 250.40959,
                            9000.0 => 270.9169715)
const GRID_STEP = parse(Float64, get(ENV, "GRID_STEP_CM", "0.5"))
const ONLY = get(ENV, "ONLY_CASE", "all")
const PAD_CM = maximum(CUTOFFS)

println("case,species,p_hpa,T_K,cutoff_cm,mean_rel_to_1000,l1_rel_to_1000,max_abs_rel_to_1000")
for (tag, species, mol, lo, hi) in CASES
    ONLY == "all" || ONLY == tag || continue
    lines = read_hitran(artifact(species), mol=mol, iso=1,
                        ν_min=lo - PAD_CM, ν_max=hi + PAD_CM)
    grid = lo:GRID_STEP:hi
    for p in PRESSURES_HPA
        T = TEMPERATURES_K[p]
        ref_model = make_hitran_model(lines, Voigt(); wing_cutoff=last(CUTOFFS),
            CEF=HumlicekWeidemann32SDErrorFunction(), architecture=GPU())
        ref = Float64.(collect(compute_absorption_cross_section(ref_model, grid, p, T)))
        ref_mean = sum(ref) / length(ref)
        ref_l1 = sum(abs, ref)
        for cutoff in CUTOFFS
            model = make_hitran_model(lines, Voigt(); wing_cutoff=cutoff,
                CEF=HumlicekWeidemann32SDErrorFunction(), architecture=GPU())
            σ = Float64.(collect(compute_absorption_cross_section(model, grid, p, T)))
            mean_rel = ref_mean == 0 ? 0.0 : (sum(σ) / length(σ) - ref_mean) / ref_mean
            l1_rel = ref_l1 == 0 ? 0.0 : sum(abs.(σ .- ref)) / ref_l1
            scale = maximum(ref)
            max_abs_rel = scale == 0 ? 0.0 : maximum(abs.(σ .- ref)) / scale
            @printf("%s,%s,%.1f,%.6f,%d,%.9e,%.9e,%.9e\n",
                    tag, species, p, T, cutoff, mean_rel, l1_rel, max_abs_rel)
        end
    end
end
