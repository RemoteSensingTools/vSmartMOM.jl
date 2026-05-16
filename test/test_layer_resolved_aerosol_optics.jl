using Test
using vSmartMOM
using vSmartMOM.CoreRT

function _layer_resolved_aerosol_config()
    return Dict(
        "radiative_transfer" => Dict(
            "spec_bands" => ["[13000.0, 13000.1]"],
            "surface" => ["LambertianSurfaceScalar(0.1)"],
            "quadrature_type" => "GaussLegQuad()",
            "polarization_type" => "Stokes_I()",
            "max_m" => 1,
            "Δ_angle" => 0.0,
            "l_trunc" => 12,
            "depol" => 0.0,
            "float_type" => "Float64",
            "architecture" => "CPU()",
        ),
        "geometry" => Dict(
            "sza" => 30.0,
            "vza" => [0.0],
            "vaz" => [0.0],
            "obs_alt" => 1000.0,
        ),
        "atmospheric_profile" => Dict(
            "T" => [250.0, 285.0],
            "p" => [100.0, 550.0, 1000.0],
            "profile_reduction" => -1,
        ),
        "scattering" => Dict(
            "aerosols" => [Dict(
                "τ_ref" => 0.02,
                "phase_function" => "SyntheticPolarizedHenyeyGreensteinPhaseFunction(g=0.45, polarization_fraction=0.0)",
                "ssa" => 0.92,
                "p₀" => 550.0,
                "σp" => 200.0,
            )],
            "r_max" => 1.0,
            "nquad_radius" => 8,
            "λ_ref" => 0.55,
            "decomp_type" => "NAI2()",
        ),
    )
end

function _with_layer_resolved_aerosol_optics(model)
    n_layers = size(model.τ_rayl[1], 2)
    base_optics = model.aerosol_optics[1][1]
    layered = LayerResolvedAerosolOptics([deepcopy(base_optics) for _ in 1:n_layers])

    aerosol_optics = [[layered]]
    aerosols = AerosolState(aerosol_optics, model.τ_aer)
    optics = Optics(model.optics.rayleigh, aerosols, model.τ_abs, model.τ_rayl)

    return RTModel(model.architecture, model.solver, model.numerics,
                   model.geometry, model.quad_points, model.atmosphere,
                   optics, model.surfaces, model.sources)
end

@testset "LayerResolvedAerosolOptics" begin
    params = parameters_from_dict(_layer_resolved_aerosol_config())
    model = model_from_parameters(params)
    layered_model = _with_layer_resolved_aerosol_optics(model)

    @test layered_model.aerosol_optics[1][1] isa LayerResolvedAerosolOptics
    @test length(layered_model.aerosol_optics[1][1]) == size(model.τ_rayl[1], 2)

    R_ref, T_ref = rt_run(model)
    R_layered, T_layered = rt_run(layered_model)

    @test R_layered == R_ref
    @test T_layered == T_ref

    short_layers = LayerResolvedAerosolOptics([deepcopy(model.aerosol_optics[1][1])])
    bad_aerosols = AerosolState([[short_layers]], model.τ_aer)
    bad_optics = Optics(model.optics.rayleigh, bad_aerosols, model.τ_abs, model.τ_rayl)
    bad_model = RTModel(model.architecture, model.solver, model.numerics,
                        model.geometry, model.quad_points, model.atmosphere,
                        bad_optics, model.surfaces, model.sources)

    if size(model.τ_rayl[1], 2) != 1
        @test_throws ArgumentError rt_run(bad_model)
    end
end
