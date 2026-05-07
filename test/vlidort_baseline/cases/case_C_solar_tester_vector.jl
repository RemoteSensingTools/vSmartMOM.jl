using Test
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering

include(joinpath(REF_DIR, "solar_tester_atmosphere.jl"))
include(joinpath(REF_DIR, "solar_tester_problemIII_aerosol.jl"))
include(joinpath(REF_DIR, "solar_tester_vector_truth.jl"))

const SOLAR_TESTER_VECTOR_YAML = joinpath(CONFIG_DIR, "solar_tester_vector.yaml")

# Aerosol parameters mirroring V2p8p3_solar_tester.f90:282-322 (vector path).
const VEC_AEROSOL_OMEGA      = 0.99999
const VEC_AEROSOL_TAU_TOTAL  = 0.5
# NMOMENTS truncation to match VLIDORT internal cap min(2·NSTREAMS−1, n_moms).
# NSTREAMS=8 ⇒ NMOMENTS = 15. We provide β/α/γ/δ/ϵ/ζ of length 16 (L=0..15).
const VEC_AEROSOL_NMOMENTS   = 15

"Build the Problem III `AerosolOptics` from the parsed `PROBLEMIII_*`
arrays.

    a1   → β   (B[1,1] phase function)
    b1   → γ   (B[1,2] = B[2,1])
    a2   → α   (B[2,2])
    a3   → ζ   (B[3,3])
    b2   → −ϵ  (B[3,4]; sign-flipped per VLIDORT GREEKMAT(.,12) = −b2)
    a4   → δ   (B[4,4])

Important: VLIDORT's vector solar_tester driver later applies
`SMASK(2:3)=-1` when converting `ProblemIII.Moms` into `GREEKMAT`.
The committed `PROBLEMIII_b1` values are the pre-SMASK file values, not the
final VLIDORT `GREEKMAT(1,2)` values. Passing `b1` through here keeps the
aerosol in the same internal convention as vSmartMOM's Rayleigh
`get_greek_rayleigh`; the VLIDORT truth Q/U values are sign-flipped at
comparison."
function solar_tester_vector_aerosol_optics()
    Lp1 = VEC_AEROSOL_NMOMENTS + 1   # L = 0..VEC_AEROSOL_NMOMENTS
    @assert length(PROBLEMIII_a1) >= Lp1 "PROBLEMIII has $(length(PROBLEMIII_a1)) moments, need ≥$(Lp1)"
    α = collect(Float64, PROBLEMIII_a2[1:Lp1])
    β = collect(Float64, PROBLEMIII_a1[1:Lp1])
    γ = collect(Float64, PROBLEMIII_b1[1:Lp1])
    δ = collect(Float64, PROBLEMIII_a4[1:Lp1])
    ϵ = .-collect(Float64, PROBLEMIII_b2[1:Lp1])
    ζ = collect(Float64, PROBLEMIII_a3[1:Lp1])
    greek = Scattering.GreekCoefs(α, β, γ, δ, ϵ, ζ)
    return Scattering.AerosolOptics(greek_coefs=greek, ω̃=VEC_AEROSOL_OMEGA,
                                    k=1.0, fᵗ=0.0)
end

"VLIDORT-style aerosol extinction per layer: parcel × (h[n−1] − h[n])."
function solar_tester_vector_aerosol_extinction()
    h = vcat(60.0, SOLAR_TESTER_HEIGHT_KM)
    n6 = 23 - 6
    parcel = VEC_AEROSOL_TAU_TOTAL / (h[n6 + 1] - h[end])
    aerext = zeros(Float64, 23)
    for n in (n6 + 1):23
        aerext[n] = parcel * (h[n] - h[n + 1])
    end
    return aerext
end

"Override per-layer optical properties on a built RTModel with the
vector solar_tester atmosphere + Problem III aerosol."
function inject_solar_tester_vector_optics!(model)
    aerext = solar_tester_vector_aerosol_extinction()
    model.aerosol_optics[1][1] = solar_tester_vector_aerosol_optics()
    @assert size(model.τ_rayl[1], 2) == 23 "expected 23 layers, got $(size(model.τ_rayl[1], 2))"
    for n in 1:23
        ext = SOLAR_TESTER_MOLEXT[n]
        ssa = SOLAR_TESTER_MOLOMG[n]
        model.τ_rayl[1][:, n] .= ssa * ext
        model.τ_abs[1][:, n]  .= (1 - ssa) * ext
        model.τ_aer[1][1, n]   = aerext[n]
    end
    return model
end

"Run vSmartMOM at one (SZA, RAZ). Returns `(R, T)` — TOA upwelling
reflectance and BOA downwelling transmittance with shape `(n_vza, 3, n_spec)`
for Stokes_IQU. `spec`, when supplied, overrides quadrature/float/arch."
function solar_tester_vector_run_at(sza::Real, raz::Real; spec = nothing)
    params = parameters_from_yaml(SOLAR_TESTER_VECTOR_YAML)
    spec === nothing || apply_overrides!(params, spec)
    params.sza = float(sza)
    fill!(params.vaz, float(raz))
    model = model_from_parameters(params)
    inject_solar_tester_vector_optics!(model)
    out = CoreRT.rt_run(model, i_band=1)
    return Array(out[1]), Array(out[2])
end

"Geometry index in 1..36 for (i_sza, i_vza, i_raz). RAZ inner, VZA mid, SZA outer."
geom_index_vec(i_sza, i_vza, i_raz) = (i_sza - 1) * 9 + (i_vza - 1) * 3 + i_raz

for spec in axis_specs()
    @testset "VLIDORT baseline: Case C — solar_tester vector Stokes-3 (Task 1) [$(spec_tag(spec))]" begin
        # vSmartMOM internal convention is kept consistent between Rayleigh and
        # imported aerosol γ. Truth Q/U from VLIDORT's convention is sign-flipped
        # at comparison.
        #
        # VLIDORT vector solar_tester cfg/reference data: USER_RELAZMS = 10°,
        # 90°, 170°; geometry index 1 is raz=10°.
        sza = 35.0
        raz = SOLAR_TESTER_VECTOR_RAZ_DEG[1]
        i_sza = 1
        i_raz = 1
        task_idx = 1
        toa_level = 1
        boa_level = length(SOLAR_TESTER_VECTOR_TAU_LEVELS)
        s = tol_scale(spec)

        @info "Case C: running solar_tester vector [$(spec_tag(spec))]" sza raz
        R, T = solar_tester_vector_run_at(sza, raz; spec=spec)

        # Per-Stokes representative scale = max truth magnitude in this set.
        # `atol = 100·eps(FT)·scale` is the FT noise floor.
        stokes_map = [(1, :I), (2, :Q), (3, :U)]
        truth_scales = Dict(s_sym => max(maximum(abs.(SOLAR_TESTER_VECTOR_STOKES[s_sym][:, toa_level, 1, task_idx])),
                                          maximum(abs.(SOLAR_TESTER_VECTOR_STOKES[s_sym][:, boa_level, 2, task_idx])))
                            for (_, s_sym) in stokes_map)
        atol_per_stokes = Dict(s_sym => 100 * eps(spec.float) * truth_scales[s_sym]
                               for (_, s_sym) in stokes_map)

        function compare(modeled, truth, label, sza, vza, raz, geom, stokes_sym)
            # Regularized rel-err: |err| / (|truth| + atol). Naturally
            # bounded when |truth| → 0 (saturates at |err|/atol).
            atol = atol_per_stokes[stokes_sym]
            re = abs(modeled - truth) / (abs(truth) + atol)
            @info "  geom=$geom $label $stokes_sym" sza vza raz modeled=round(modeled, sigdigits=6) truth=round(truth, sigdigits=6) rel_err=round(re, sigdigits=3)
            return re
        end
        rel_up = Dict{Symbol, Vector{Float64}}(:I => Float64[], :Q => Float64[], :U => Float64[])
        rel_dn = Dict{Symbol, Vector{Float64}}(:I => Float64[], :Q => Float64[], :U => Float64[])

        vza_list = SOLAR_TESTER_VECTOR_VZA_DEG
        for (i_vza, vza) in enumerate(vza_list)
            geom = geom_index_vec(i_sza, i_vza, i_raz)
            for (s_idx, s_sym) in stokes_map
                truth_up = SOLAR_TESTER_VECTOR_STOKES[s_sym][geom, toa_level, 1, task_idx]
                truth_dn = SOLAR_TESTER_VECTOR_STOKES[s_sym][geom, boa_level, 2, task_idx]
                # Convert VLIDORT truth to the vSmartMOM internal convention for
                # comparison. The ProblemIII.Moms b1 file values are pre-SMASK;
                # the VLIDORT driver applies SMASK before solving.
                if s_sym in (:Q, :U)
                    truth_up = -truth_up
                    truth_dn = -truth_dn
                end
                re_u = compare(R[i_vza, s_idx, 1], truth_up, "TOA-up", sza, vza, raz, geom, s_sym)
                re_d = compare(T[i_vza, s_idx, 1], truth_dn, "BOA-dn", sza, vza, raz, geom, s_sym)
                isnan(re_u) || push!(rel_up[s_sym], re_u)
                isnan(re_d) || push!(rel_dn[s_sym], re_d)
            end
        end

        println("\nCase C [$(spec_tag(spec))] — solar_tester vector Task 1, sza=$sza raz=$raz:")
        for s_sym in (:I, :Q, :U)
            if !isempty(rel_up[s_sym])
                println("  TOA-up $s_sym: median=$(round(Statistics.median(rel_up[s_sym]), sigdigits=3))  max=$(round(maximum(rel_up[s_sym]), sigdigits=3))")
            end
            if !isempty(rel_dn[s_sym])
                println("  BOA-dn $s_sym: median=$(round(Statistics.median(rel_dn[s_sym]), sigdigits=3))  max=$(round(maximum(rel_dn[s_sym]), sigdigits=3))")
            end
        end

        # FT-objective gate: √(floor² + ft_margin²) with ft_margin=K·√(eps(FT)).
        # Quadrature combination of physics-model truncation + FT roundoff.
        # K=300 calibrated empirically to ~30-doubling F32 noise floor.
        # F64: ~4.5e-6 (negligible); F32: ~0.1 (dominates).
        ft_margin = 300 * sqrt(eps(spec.float))
        for s_sym in (:I, :Q, :U)
            @test maximum(rel_up[s_sym]) < sqrt((1e-3)^2 + ft_margin^2)
            @test maximum(rel_dn[s_sym]) < sqrt((2e-3)^2 + ft_margin^2)
        end
        @info "Case C maxima" TOA_I_max=maximum(rel_up[:I]) BOA_I_max=maximum(rel_dn[:I]) Q_TOA_max=maximum(rel_up[:Q]) Q_BOA_max=maximum(rel_dn[:Q]) U_TOA_max=maximum(rel_up[:U]) U_BOA_max=maximum(rel_dn[:U])
    end
end
