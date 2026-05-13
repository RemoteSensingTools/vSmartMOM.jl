# =========================================================================
# v0.7 Phase A — ThermalEmission per-layer Planck volume source
#
# Unit-test coverage:
#   T-A0a  Construction (default, B_layer kwarg, T_layers+ν_grid convenience)
#   T-A0b  prepare_source produces correct PreparedThermalEmission
#   T-A0c  SourceSet composition (SolarBeam + ThermalEmission)
#   T-A1   contribute!(PreparedThermalEmission, ...) at m=0 writes
#          2π·(1-ϖ)·B·(1-exp(-dτ/μ)) into I-rows of j₀±, leaves Q/U/V at 0
#   T-A2   contribute!(...) at m>0 is a no-op (TIR weight rule)
#   T-A3   contribute!(NoSource, ...) is a no-op
#   T-A4   contribute!(SourceSet{(SolarBeam, ThermalEmission)}, ...) propagates
#          to the ThermalEmission member at m=0; no-op at m>0
#   T-A5   prepare_source rejects mismatched B_layer column count
#
# The TIR weight rule (overrides REFACTOR_SPEC_v6 §2.11) is enforced at
# injection time: thermal is isotropic, so it fires only at m=0, and the
# explicit 2π factor undoes the uniform 0.5/π weight applied downstream by
# `postprocessing_vza!`. The test cases below pin down this contract so a
# future regression cannot silently break it.
# =========================================================================

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.SolarModel: planck_spectrum_wn
using Test

@testset "v0.7 sources — ThermalEmission volume Planck source" begin

    @testset "T-A0a: ThermalEmission construction" begin
        # Default (placeholder, no thermal payload)
        te_default = ThermalEmission()
        @test te_default isa ThermalEmission
        @test te_default.B_layer === nothing

        # Pre-computed B_layer kwarg
        B = ones(Float64, 4, 7)            # 4 layers × 7 spectral pts
        te_pre = ThermalEmission(B_layer = B)
        @test te_pre.B_layer === B

        # T_layers + ν_grid convenience constructor → Planck-shaped B_layer
        T_layers = [200.0, 250.0, 300.0]
        ν_grid = collect(800.0:50.0:1200.0)        # cm⁻¹
        te_conv = ThermalEmission(T_layers, ν_grid)
        @test te_conv isa ThermalEmission
        @test size(te_conv.B_layer) == (length(T_layers), length(ν_grid))
        # Spot-check: row iz must equal planck_spectrum_wn(T_layers[iz], ν_grid)
        for iz in eachindex(T_layers)
            B_ref = planck_spectrum_wn(T_layers[iz], collect(float.(ν_grid)))
            @test te_conv.B_layer[iz, :] ≈ B_ref
        end
        # Source AD-mode default is analytic
        @test source_ad_mode(te_conv) === AnalyticSourceJacobian()
    end

    @testset "T-A0b: prepare_source(ThermalEmission, ...) precision conversion" begin
        nLayers, nSpec = 3, 5
        B = randn(Float64, nLayers, nSpec) .+ 100.0          # plausible Planck values
        te = ThermalEmission(B_layer = B)

        prep64 = prepare_source(te, Float64, 3, nSpec, identity)
        @test prep64 isa PreparedThermalEmission
        @test size(prep64.B_layer) == (nLayers, nSpec)
        @test eltype(prep64.B_layer) === Float64
        @test prep64.B_layer == B                            # exact, not just ≈

        prep32 = prepare_source(te, Float32, 3, nSpec, identity)
        @test eltype(prep32.B_layer) === Float32
        @test prep32.B_layer ≈ Float32.(B)

        # Default ThermalEmission() → 1×nSpec zero placeholder
        prep_zero = prepare_source(ThermalEmission(), Float64, 3, nSpec, identity)
        @test prep_zero isa PreparedThermalEmission
        @test all(iszero, prep_zero.B_layer)
        @test size(prep_zero.B_layer) == (1, nSpec)

        @test source_ad_mode(prep64) === AnalyticSourceJacobian()
        @test has_thermal_emission(prep64) == true
        @test has_thermal_emission(NoSource()) == false
    end

    @testset "T-A0c: SourceSet composition with ThermalEmission" begin
        nLayers, nSpec = 3, 5
        B = ones(Float64, nLayers, nSpec)
        combo = SolarBeam() + ThermalEmission(B_layer = B)
        @test combo isa SourceSet
        @test length(combo) == 2

        prep_combo = prepare_sources(combo, Float64, 3, nSpec, identity)
        @test prep_combo isa SourceSet
        @test length(prep_combo) == 2
        # `has_thermal_emission` is a *prepared-set* predicate (mirrors the
        # `has_surface_sif`/`has_solar_beam` convention) — call it after
        # prepare_sources, not on the raw composition.
        @test has_thermal_emission(prep_combo) == true
        # Predicate is false when only SolarBeam present
        prep_solo = prepare_sources(SolarBeam(), Float64, 3, nSpec, identity)
        @test has_thermal_emission(prep_solo) == false
    end

    @testset "T-A1, T-A2: contribute! at m=0 (write) and m>0 (no-op)" begin
        # v0.7 Phase A.2a — thermal now writes to `added_layer.j₀_by_src[:thermal]`
        # slot (a `SourceSlot` carrier), not the legacy `j₀⁺/j₀⁻` fields. Build a
        # minimal `added_layer` stub that exposes the same `.j₀_by_src.thermal`
        # access pattern the contribute! function expects.
        FT = Float64
        nSpec = 4
        Nquad = 3
        nStokes = 3
        NquadN = Nquad * nStokes
        thermal_slot = CoreRT.SourceSlot(
            j₀⁺     = zeros(FT, NquadN, 1, nSpec),
            j₀⁻     = zeros(FT, NquadN, 1, nSpec),
            dbl_j₁⁺ = zeros(FT, NquadN, 1, nSpec),
            dbl_j₁⁻ = zeros(FT, NquadN, 1, nSpec),
            expk    = ones(FT, nSpec),
        )
        added = (j₀_by_src = (thermal = thermal_slot,),)

        struct _PolStubA; n::Int; end
        pol = _PolStubA(nStokes)

        # quad_points stub carrying qp_μ
        qp_μ = FT[0.2, 0.6, 0.9]
        struct _QPStubA{T}; qp_μ::T; end
        qp = _QPStubA(qp_μ)

        # 2 layers × 4 spectral pts; non-uniform Planck values so we can
        # check that iz indexing works correctly.
        B = FT[10.0 20.0 30.0 40.0;
                5.0 15.0 25.0 35.0]
        prep_te = PreparedThermalEmission{FT, typeof(B)}(B)

        # Per-spectral-point ϖ and dτ vectors
        ϖ_λ = FT[0.0, 0.2, 0.5, 0.8]                   # absorption fractions = [1.0, 0.8, 0.5, 0.2]
        dτ_λ = FT[0.1, 0.5, 1.0, 2.0]
        arch = vSmartMOM.Architectures.CPU()

        # T-A1: contribute! at m=0 writes thermal source into the :thermal slot's
        # I-rows of j₀±. Legacy `j₀⁺/j₀⁻` fields (not exposed on this stub) are
        # untouched — the slot is a separate buffer pool.
        CoreRT.contribute!(prep_te, added, ϖ_λ, dτ_λ, 1, 0, pol, qp, arch)

        # Expected: for each stream iμ, slot.j₀±[i_I, 1, n] = 2π·(1-ϖ_λ[n])·B[1,n]·(1 - exp(-dτ_λ[n]/μ[iμ]))
        for iμ in 1:Nquad
            i_I = (iμ - 1) * nStokes + 1
            for n in 1:nSpec
                expected = 2π * (1 - ϖ_λ[n]) * B[1, n] *
                           (1 - exp(-dτ_λ[n] / qp_μ[iμ]))
                @test thermal_slot.j₀⁺[i_I, 1, n] ≈ expected
                @test thermal_slot.j₀⁻[i_I, 1, n] ≈ expected
            end
            # Q/U/V rows must remain zero (thermal is unpolarized)
            for iS in 2:nStokes
                iS_row = (iμ - 1) * nStokes + iS
                @test all(thermal_slot.j₀⁺[iS_row, 1, :] .== 0)
                @test all(thermal_slot.j₀⁻[iS_row, 1, :] .== 0)
            end
        end

        # T-A1 with iz=2 picks the second row of B (different magnitudes).
        # The slot is zeroed at the start of contribute!, so we don't need
        # to reset by hand — but a second call with iz=2 fully overwrites
        # the iz=1 contribution.
        CoreRT.contribute!(prep_te, added, ϖ_λ, dτ_λ, 2, 0, pol, qp, arch)
        i_I_stream1 = 1
        n = 1
        expected_iz2 = 2π * (1 - ϖ_λ[n]) * B[2, n] *
                        (1 - exp(-dτ_λ[n] / qp_μ[1]))
        @test thermal_slot.j₀⁺[i_I_stream1, 1, n] ≈ expected_iz2
        @test thermal_slot.j₀⁺[i_I_stream1, 1, n] != 0
        @test !(thermal_slot.j₀⁺[i_I_stream1, 1, n] ≈ 2π * (1 - ϖ_λ[n]) * B[1, n] *
                                                       (1 - exp(-dτ_λ[n] / qp_μ[1])))

        # T-A2: contribute! at m>0 is a no-op (TIR weight rule).
        # Re-zero the slot manually (contribute! at m>0 doesn't enter the
        # zero-then-write path) and verify it stays zero across multiple m.
        thermal_slot.j₀⁺ .= 0; thermal_slot.j₀⁻ .= 0
        for m_test in 1:4
            CoreRT.contribute!(prep_te, added, ϖ_λ, dτ_λ, 1, m_test, pol, qp, arch)
        end
        @test all(thermal_slot.j₀⁺ .== 0)
        @test all(thermal_slot.j₀⁻ .== 0)
    end

    @testset "T-A3: contribute! NoSource is a no-op" begin
        FT = Float64
        Nquad = 2; nStokes = 3; NquadN = Nquad * nStokes
        nSpec = 3
        # Build a slot we can check stays untouched; NoSource doesn't read it
        # but we want to confirm it isn't accidentally zeroed either.
        thermal_slot = CoreRT.SourceSlot(
            j₀⁺     = ones(FT, NquadN, 1, nSpec),   # nonzero sentinel
            j₀⁻     = ones(FT, NquadN, 1, nSpec),
            dbl_j₁⁺ = zeros(FT, NquadN, 1, nSpec),
            dbl_j₁⁻ = zeros(FT, NquadN, 1, nSpec),
            expk    = ones(FT, nSpec),
        )
        added = (j₀_by_src = (thermal = thermal_slot,),)
        struct _PolStubB; n::Int; end
        pol = _PolStubB(nStokes)
        qp = (qp_μ = FT[0.4, 0.8],)
        arch = vSmartMOM.Architectures.CPU()

        CoreRT.contribute!(NoSource(), added,
                            FT[0.1, 0.2, 0.3], FT[1.0, 1.0, 1.0],
                            1, 0, pol, qp, arch)
        # NoSource doesn't touch the slot — sentinel ones remain.
        @test all(thermal_slot.j₀⁺ .== 1)
        @test all(thermal_slot.j₀⁻ .== 1)
    end

    @testset "T-A4: contribute! SourceSet routes to ThermalEmission member" begin
        FT = Float64
        Nquad = 2; nStokes = 3; NquadN = Nquad * nStokes
        nSpec = 3
        thermal_slot = CoreRT.SourceSlot(
            j₀⁺     = zeros(FT, NquadN, 1, nSpec),
            j₀⁻     = zeros(FT, NquadN, 1, nSpec),
            dbl_j₁⁺ = zeros(FT, NquadN, 1, nSpec),
            dbl_j₁⁻ = zeros(FT, NquadN, 1, nSpec),
            expk    = ones(FT, nSpec),
        )
        added = (j₀_by_src = (thermal = thermal_slot,),)
        struct _PolStubC; n::Int; end
        pol = _PolStubC(nStokes)
        qp_μ = FT[0.4, 0.8]
        qp = (qp_μ = qp_μ,)
        arch = vSmartMOM.Architectures.CPU()

        B = FT[100.0 200.0 300.0]                              # 1 layer × 3 spec
        prep_combo = prepare_sources(
            SolarBeam() + ThermalEmission(B_layer = B),
            FT, nStokes, nSpec, identity)
        ϖ_λ = FT[0.2, 0.4, 0.6]
        dτ_λ = FT[1.0, 1.0, 1.0]

        CoreRT.contribute!(prep_combo, added, ϖ_λ, dτ_λ, 1, 0, pol, qp, arch)

        # Thermal member fires into the :thermal slot — SolarBeam member is a
        # no-op for the volume contribute! seam (its SFI math lives inside
        # elemental!).
        for iμ in 1:Nquad
            i_I = (iμ - 1) * nStokes + 1
            for n in 1:nSpec
                expected = 2π * (1 - ϖ_λ[n]) * B[1, n] *
                           (1 - exp(-dτ_λ[n] / qp_μ[iμ]))
                @test thermal_slot.j₀⁺[i_I, 1, n] ≈ expected
                @test thermal_slot.j₀⁻[i_I, 1, n] ≈ expected
            end
        end

        # And at m=1 the whole SourceSet is a no-op (thermal gated, solar SFI
        # no-op here). Re-zero manually since m>0 contribute! skips the
        # zero-then-write path.
        thermal_slot.j₀⁺ .= 0; thermal_slot.j₀⁻ .= 0
        CoreRT.contribute!(prep_combo, added, ϖ_λ, dτ_λ, 1, 1, pol, qp, arch)
        @test all(thermal_slot.j₀⁺ .== 0)
        @test all(thermal_slot.j₀⁻ .== 0)
    end

    @testset "T-A5: prepare_source rejects shape mismatch on B_layer" begin
        # B_layer second-axis must equal nSpec
        B_bad = ones(Float64, 3, 4)                # 4 spec columns
        @test_throws ErrorException prepare_source(
            ThermalEmission(B_layer = B_bad), Float64, 3, 5, identity)
    end

    @testset "T-A6: thermal RT is SZA-independent (per-source doubling correctness)" begin
        # End-to-end correctness check for the v0.7 Phase A.2a refactor.
        #
        # PHYSICS PROPERTY UNDER TEST. Thermal volume emission is isotropic
        # and self-generated — it has no incident-direction dependence. So the
        # emergent thermal radiance at TOA, viewed at fixed VZA, must NOT
        # depend on the solar zenith angle μ₀ that the model was configured
        # with. Pre-Phase-A.2a, the thermal contribution lived in the legacy
        # solar j₀± and was multiplied by `expk_solar = exp(-dτ/μ₀)` during
        # doubling — so changing SZA changed the result. After A.2a, thermal
        # lives in its own slot with `expk = ones` and the doubling math
        # collapses to the Fortran TIR recipe (rt_doubling.f90:191-197).
        #
        # TEST DESIGN.
        #   1. Load PureRayleighParameters (1-band Rayleigh scattering, no
        #      surface reflection, full Stokes IQUV).
        #   2. Build the model, derive an isothermal B_layer from the
        #      atmospheric T profile + ν grid.
        #   3. Run `rt_run` thermal-only at SZA = 30°  → R_30
        #   4. Re-run at SZA = 60°                     → R_60
        #   5. Assert that R_30 ≈ R_60 to a tight tolerance (the doubling
        #      bug would produce a > 10× difference at thick layers).
        params = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        params.architecture = vSmartMOM.Architectures.CPU()

        # Pre-compute thermal B_layer matching the model's spectral grid.
        # NOTE 1. The test YAML uses a visible band (~19417 cm⁻¹). At realistic
        # atmospheric T ~ 250 K the Planck radiance is ~ 10⁻⁴⁰ at this ν —
        # well below Float64 ε. We therefore use a *synthetic* B_layer with
        # non-trivially nonzero magnitudes; the math correctness of the
        # SZA-independence property does not depend on whether B is
        # physically realistic at this band.
        # NOTE 2. PureRayleigh has ϖ = 1 (pure scattering, no absorption) at
        # every layer, so (1 − ϖ_λ) ≡ 0 inside `contribute!(thermal)` and the
        # thermal source would vanish. We inject per-layer τ_abs so each
        # layer absorbs ~ 10% of its scattering depth → ϖ ≈ 0.9 and the
        # thermal contribution is non-zero.
        params.sza = 30.0
        model_30 = model_from_parameters(params)
        FT_mod = CoreRT.float_type(model_30)
        spec_band = collect(CoreRT.get_spec_bands(model_30)[1])
        nSpec = length(spec_band)
        nLayers = length(model_30.profile.T)
        B_layer = zeros(FT_mod, nLayers, nSpec)
        for iz in 1:nLayers
            # Synthetic B with non-uniform per-layer magnitude — captures the
            # iz indexing path while staying in numerically stable range.
            B_layer[iz, :] .= FT_mod(0.1) + FT_mod(0.01) * (iz - 1)
        end
        thermal_only = ThermalEmission(B_layer = B_layer)

        # Inject some absorption so the thermal source has a non-zero (1−ϖ)
        # factor. `model.τ_abs[iBand]` is shape `(nSpec, nLayers)` and is
        # writable post-construction (see rt_run.jl docstrings).
        function _inject_τ_abs!(m)
            for iB in 1:length(m.τ_abs)
                m.τ_abs[iB] .= FT_mod(0.05)
            end
        end
        _inject_τ_abs!(model_30)

        # Solar-only baseline for comparison (default model.sources = SolarBeam())
        _R_solar_baseline = rt_run(model_30)[1]

        # Thermal-only RT at SZA = 30
        R_30 = rt_run(model_30; sources = thermal_only)[1]

        # Same configuration at SZA = 60
        params.sza = 60.0
        model_60 = model_from_parameters(params)
        _inject_τ_abs!(model_60)
        R_60 = rt_run(model_60; sources = thermal_only)[1]

        # Property check: emergent thermal radiance is independent of SZA.
        # Tolerance: the thermal Planck radiance at this band (~19417 cm⁻¹ for
        # T ~ 220-285 K) is tiny but non-zero. Use a relative tolerance against
        # the max R magnitude; the doubling bug would scale the result by
        # `exp(-τ/μ_30) / exp(-τ/μ_60)` which is ~10× for atmospheric τ ~ 1.
        max_R = maximum(abs, R_30)
        max_R > 0 || @warn "Thermal R_30 is identically zero — band/T mismatch?"
        rel_diff = maximum(abs, R_30 .- R_60) / max(max_R, eps(FT_mod))
        @test rel_diff < 1e-6

        # Sanity: thermal-only must differ from solar-only (different physics
        # entirely; this catches a wiring error where thermal isn't reaching
        # the kernel at all).
        @test R_30 != _R_solar_baseline
        # Sanity: thermal contribution at this visible band should be small
        # but positive in I (B(T ≈ 250K, ν = 19417 cm⁻¹) is tiny, but the
        # Rayleigh-scattered thermal emission is nonzero).
        @test maximum(abs, R_30[:, 1, :]) > 0
    end

    @testset "T-A7: thick isothermal magnitude — R/B = 1 across m_max (slot-leak regression)" begin
        # MAGNITUDE check that catches the per-source slot-leak bug.
        #
        # For an optically thick isothermal column with no scattering
        # (ϖ = 0) the emergent thermal radiance at TOA must equal
        # `B(T, ν)` to machine precision — independent of how many
        # Fourier moments the solver runs. Before the slot-reset fix in
        # `rt_kernel.jl`, the per-source `slot.j₀±` from m=0 leaked into
        # m>0 doublings and produced `R/B = 1 + 2·m_max` (vaz=0 phase
        # gathers all moments coherently): m_max=0 → 1, m_max=1 → 3,
        # m_max=2 → 5, …
        #
        # This test pins R/B = 1 at three different forced m_max values.
        #
        # Build PureRayleighParameters in IR; force τ_abs huge per layer
        # so the column is fully opaque; T_iso so the only relevant
        # quantity is B(T_iso, ν).
        params = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        params.architecture = vSmartMOM.Architectures.CPU()
        params.vza = [0.0]
        params.vaz = [0.0]
        T_iso = 250.0

        for m_force in (0, 1, 2)
            model = model_from_parameters(params)
            ν = collect(CoreRT.get_spec_bands(model)[1])
            # Force the Fourier-moment count for this Pkg.test invocation.
            model.solver.m_max_bands[1] = m_force
            # Fully opaque per-layer absorption so the atmosphere is a
            # blackbody at T_iso to numerical precision.
            fill!(model.τ_abs[1], 100.0)

            FT_mod = CoreRT.float_type(model)
            T_layers = fill(FT_mod(T_iso), length(model.profile.T))
            thermal = ThermalEmission(T_layers, ν)

            R, _T_out = rt_run(model; sources = thermal)
            B_expected = vSmartMOM.SolarModel.planck_spectrum_wn(T_iso, Float64.(ν))[1]
            ratio = R[1, 1, 1] / B_expected
            @test isapprox(ratio, 1.0; atol = 1e-3)
        end
    end

end
