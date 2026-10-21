using Test
using vSmartMOM
using vSmartMOM.Absorption
using vSmartMOM.CoreRT
import AtmosphericAbsorption

function _synthetic_h2o_model(::Type{FT}) where {FT<:AbstractFloat}
    lines = AtmosphericAbsorption.LineDatabase(
        mol=Int32[1],
        iso=Int32[1],
        ν0=FT[1000],
        S=FT[1e-21],
        E_lower=FT[0],
        g_upper=FT[1],
        γ_air=FT[0.03],
        γ_self=FT[0.30],
        n_air=FT[0.7],
        δ_air=FT[-0.01],
        n_self=FT[0.5],
        δ_self=FT[0.04],
        molar_mass=FT[18],
        meta=AtmosphericAbsorption.SourceMetadata(
            "synthetic H2O self-broadening regression", 296.0, 1013.25),
    )
    flat_partition = AtmosphericAbsorption.TabulatedPF(
        FT[100, 300, 500], FT[1, 1, 1])
    return AtmosphericAbsorption.LineByLineModel(
        lines, flat_partition;
        profile=AtmosphericAbsorption.Voigt(),
        wing_cutoff=FT(10),
        vmr=zero(FT),
        architecture=AtmosphericAbsorption.CPU(),
    )
end

function _synthetic_h2o_profile(::Type{FT}) where {FT<:AbstractFloat}
    r_h2o = FT[0.2, 0.5] # deliberately wet so r and r/(1+r) are distinct
    vcd_dry = FT[1e22, 2e22]
    return CoreRT.AtmosphericProfile(
        FT[260, 280],
        FT[600, 900],
        FT[0, 0],
        FT[400, 800, 1000],
        r_h2o,
        vcd_dry,
        r_h2o .* vcd_dry,
        Dict{String,FT}(),
        FT[4000, 2000],
    )
end

@testset "H2O line self-broadening uses moist mole fraction" begin
    for FT in (Float32, Float64)
        model = _synthetic_h2o_model(FT)
        profile = _synthetic_h2o_profile(FT)
        grid = collect(FT, FT(999):FT(0.01):FT(1001))
        n_spec = length(grid)
        n_layers = length(profile.T)
        rtol = FT === Float32 ? 5f-6 : 2e-13

        # vSmartMOM's wrapper must expose AA's per-call broadener override.
        x_test = FT(0.25)
        σ_wrapped = absorption_cross_section(
            model, grid, FT(800), FT(270); vmr=x_test)
        σ_direct = AtmosphericAbsorption.compute_cross_section(
            model, grid, FT(800), FT(270); vmr=x_test)
        @test σ_wrapped == σ_direct

        expected = zeros(FT, n_spec, n_layers)
        expected_dot = zeros(FT, n_layers, n_spec, n_layers)
        wrong_dry_ratio = similar(expected)
        for iz in eachindex(profile.T)
            r = profile.vmr_h2o[iz]
            x = r / (one(FT) + r)
            σ = collect(absorption_cross_section(
                model, grid, profile.p_full[iz], profile.T[iz]; vmr=x))
            σ_wrong = collect(absorption_cross_section(
                model, grid, profile.p_full[iz], profile.T[iz]; vmr=r))
            expected[:, iz] .= σ .* profile.vcd_dry[iz] .* r
            expected_dot[iz, :, iz] .= σ .* profile.vcd_dry[iz]
            wrong_dry_ratio[:, iz] .=
                σ_wrong .* profile.vcd_dry[iz] .* r
        end

        # Forward construction path.
        τ_forward = zeros(FT, n_spec, n_layers)
        CoreRT.compute_h2o_absorption_profile!(
            τ_forward, model, grid, profile)
        @test τ_forward ≈ expected rtol=rtol
        @test !isapprox(τ_forward, wrong_dry_ratio; rtol=FT(1e-3))

        # Linearized construction must use the identical broadened cross
        # section for its forward optical depth and layer-local tangent.
        τ_linearized = zeros(FT, n_spec, n_layers)
        τ̇ = zeros(FT, n_layers, n_spec, n_layers)
        CoreRT.compute_h2o_absorption_profile!(
            τ_linearized, τ̇, 1, model, grid, profile)
        @test τ_linearized ≈ expected rtol=rtol
        @test τ̇ ≈ expected_dot rtol=rtol

        # BatchContext updates share _compute_band_absorption!; pin its H₂O
        # branch to the same result without requiring external HITRAN data.
        τ_update = zeros(FT, n_spec, n_layers)
        CoreRT._compute_band_absorption!(
            τ_update, Any[], model, grid, String[], profile, true)
        @test τ_update ≈ expected rtol=rtol
    end
end
