# =============================================================================
# prototype_inelastic.jl — smallest RRS demo (Phase 6 port from sanghavi).
#
# Builds a model from O2Parameters.yaml, constructs an RRS type placeholder
# with the Phase 1b pattern, and runs both noRS and RRS forward RT via
# rt_run_test(). Demonstrates:
#   - noRS ↔ RRS comparison (elastic + inelastic outputs)
#   - iBand-vector form rt_run_test(RS_type, model, [1, 1]) for multi-band-like
#     concatenation (here just band 1 twice for schema coverage)
#   - The minimal RRS kwargs needed to satisfy the type constructor
#     (these placeholder values are overwritten by getRamanSSProp!).
# =============================================================================

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering
using Statistics

# Load YAML files into parameter struct
parameters = parameters_from_yaml(joinpath(pkgdir(vSmartMOM), "test",
                                           "test_parameters", "O2Parameters.yaml"))
# Create model struct (precomputes optical properties) from parameters
model      = model_from_parameters(parameters)

const iBand = 1
const FT    = CoreRT.float_type(model)

# Compute all Raman properties
ν  = model.atmosphere.spec_bands[iBand]
ν̃  = mean(ν)
# Find central reference index for RRS:
i_ref = argmin(abs.(ν .- ν̃))
# Effective temperature for Raman calculations
effT = (model.profile.vcd_dry' * model.profile.T) / sum(model.profile.vcd_dry)

# N₂, O₂ molecular constants at the effective temperature
n2, o2 = InelasticScattering.getRamanAtmoConstants(ν̃, effT)

nPol  = CoreRT.polarization_type(model).n
nSpec = length(ν)
F₀   = zeros(FT, nPol, nSpec); F₀[1, :] .= one(FT)
SIF₀ = zeros(FT, nPol, nSpec)

# Minimal RRS placeholder — getRamanSSProp! fills in the real values.
RS_type = InelasticScattering.RRS(
    n2          = n2,
    o2          = o2,
    greek_raman = InelasticScattering.GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)]),
    fscattRayl  = [FT(1)],
    ϖ_Cabannes  = [FT(1)],
    ϖ_λ₁λ₀      = zeros(FT, 1),
    i_λ₁λ₀      = zeros(Int, 1),
    Z⁻⁺_λ₁λ₀    = zeros(FT, 1, 1),
    Z⁺⁺_λ₁λ₀    = zeros(FT, 1, 1),
    i_ref       = i_ref,
    n_Raman     = 0,
    F₀          = F₀,
    SIF₀        = SIF₀,
)

# Compute Raman SSA properties (populates ϖ_λ₁λ₀, Z matrices in-place)
CoreRT.getRamanSSProp!(RS_type, 1e7 / ν̃, ν)

# noRS forward RT (elastic-only, for comparison)
RnoRS, TnoRS, _, _ = CoreRT.rt_run_test(InelasticScattering.noRS{FT}(F₀ = F₀, SIF₀ = SIF₀), model, iBand)

# RRS forward RT (elastic + inelastic)
R_rrs, T_rrs, ieR, ieT = CoreRT.rt_run_test(RS_type, model, iBand)

# Multi-band-like form (both entries = band 1 — schema exerciser)
RnoRS_multi, TnoRS_multi, _, _ = CoreRT.rt_run_test(InelasticScattering.noRS{FT}(F₀ = F₀, SIF₀ = SIF₀), model, [1, 1])

@info "prototype_inelastic outputs" size_R = size(R_rrs) size_ieR = size(ieR)
