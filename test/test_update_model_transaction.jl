using Test
using vSmartMOM
using vSmartMOM.CoreRT
using Interpolations
import AtmosphericAbsorption
using Distributions
using vSmartMOM.Scattering: Aerosol

function transaction_test_parameters(reduction)
    params = read_parameters(Dict(
        "radiative_transfer" => Dict(
            "spec_bands" => ["[6199.8, 6200.0, 6200.2]",
                             "[6199.8, 6200.0, 6200.2]"],
            "surface" => ["LambertianSurfaceScalar(0.2)",
                          "LambertianSurfaceScalar(0.2)"],
            "nstreams" => 3,
            "polarization_type" => "Stokes_I()",
            "truncation" => "NoTruncation()",
            "depol" => -1,
            "float_type" => "Float64",
            "architecture" => "CPU()"),
        "geometry" => Dict(
            "sza" => 35.0,
            "vza" => [25.0],
            "vaz" => [40.0],
            "obs_alt" => [0]),
        "atmospheric_profile" => Dict(
            "T" => [270.0, 280.0],
            "p" => [100.0, 500.0, 1000.0],
            "profile_reduction" => reduction)))

    ν = range(6199.0, 6201.0; length=3)
    p = range(100.0, 1100.0; length=3)
    make_lut(T) = begin
        table = [1e-24 * (1 + 0.001pp + 0.002tt + i)
                 for i in eachindex(ν), pp in p, tt in T]
        itp = interpolate(table, BSpline(Linear()))
        vSmartMOM.Absorption.InterpolationModel(itp, 2, 1, ν, p, T)
    end

    # The first band accepts the failing 350/360 K candidate; the second band
    # throws. This exercises failure after earlier trial work has completed.
    wide_lut = make_lut(range(200.0, 400.0; length=3))
    narrow_lut = make_lut(range(200.0, 300.0; length=3))
    params.absorption_params = CoreRT.AbsorptionParameters(
        [["CO2"], ["CO2"]], [String[], String[]],
        Dict{String,Any}("CO2" => 4e-4),
        AtmosphericAbsorption.Voigt(),
        AtmosphericAbsorption.HumlicekWeideman32(),
        10.0, [Any[wide_lut], Any[narrow_lut]],
        Any[:disabled, :disabled], String[], "")
    return params
end

profile_snapshot(profile) = (
    T=copy(profile.T), p_full=copy(profile.p_full), q=copy(profile.q),
    p_half=copy(profile.p_half), vmr_h2o=copy(profile.vmr_h2o),
    vcd_dry=copy(profile.vcd_dry), vcd_h2o=copy(profile.vcd_h2o),
    vmr=deepcopy(profile.vmr), Δz=copy(profile.Δz))

geometry_snapshot(geometry) = (
    sza=geometry.sza, vza=copy(geometry.vza), vaz=copy(geometry.vaz),
    obs_alt=copy(geometry.obs_alt), sensor_levels=copy(geometry.sensor_levels),
    sensor_altitudes=copy(geometry.sensor_altitudes),
    include_toa=geometry.include_toa, include_boa=geometry.include_boa,
    toa_altitude=geometry.toa_altitude)

@testset "transactional BatchContext updates" begin
    for reduction in (-1, 1)
        @testset "profile_reduction=$reduction" begin
            params = transaction_test_parameters(reduction)
            ctx = BatchContext(params)

            # Construction owns both its live profile and remembered inputs.
            @test ctx.model.profile.T !== params.T
            @test ctx.current_T !== params.T
            @test ctx.current_T !== ctx.model.profile.T
            @test ctx.params !== params
            @test ctx.params.absorption_params.luts !== params.absorption_params.luts
            for ib in 1:ctx.n_bands
                lut = params.absorption_params.luts[ib][1]
                @test ctx.params.absorption_params.luts[ib][1].itp.coefs === lut.itp.coefs
                @test ctx.absorption_models[ib][1].itp.coefs === lut.itp.coefs
            end

            old_profile = profile_snapshot(ctx.model.profile)
            old_current = (T=copy(ctx.current_T), p=copy(ctx.current_p_half),
                           q=copy(ctx.current_q), vmr=deepcopy(ctx.current_vmr))
            old_τ_abs = map(copy, ctx.model.τ_abs)
            old_τ_rayl = map(copy, ctx.model.τ_rayl)
            old_τ_aer = map(copy, ctx.model.optics.aerosols.τ_aer)
            old_ϖ = copy(ctx.model.optics.rayleigh.ϖ_Cabannes)
            old_cabannes_β = map(g -> copy(g.β),
                                 ctx.model.optics.rayleigh.greek_cabannes)
            old_rayleigh_β = map(g -> copy(g.β),
                                 ctx.model.optics.rayleigh.greek_rayleigh)
            old_geometry = geometry_snapshot(ctx.model.geometry)
            old_radiance = [copy(rt_run(ctx.model; i_band=i).toa)
                            for i in 1:ctx.n_bands]

            @test_throws BoundsError update_model!(ctx; T=[350.0, 360.0])

            @test profile_snapshot(ctx.model.profile) == old_profile
            @test ctx.current_T == old_current.T
            @test ctx.current_p_half == old_current.p
            @test ctx.current_q == old_current.q
            @test ctx.current_vmr == old_current.vmr
            @test ctx.model.τ_abs == old_τ_abs
            @test ctx.model.τ_rayl == old_τ_rayl
            @test ctx.model.optics.aerosols.τ_aer == old_τ_aer
            @test ctx.model.optics.rayleigh.ϖ_Cabannes == old_ϖ
            @test map(g -> g.β, ctx.model.optics.rayleigh.greek_cabannes) == old_cabannes_β
            @test map(g -> g.β, ctx.model.optics.rayleigh.greek_rayleigh) == old_rayleigh_β
            @test geometry_snapshot(ctx.model.geometry) == old_geometry
            @test [rt_run(ctx.model; i_band=i).toa for i in 1:ctx.n_bands] == old_radiance

            # A valid retry must match a fresh model, and later caller mutation
            # must not alter either remembered or live state.
            retry_T = [275.0, 285.0]
            update_model!(ctx; T=retry_T)
            fresh_params = transaction_test_parameters(reduction)
            fresh_params.T .= retry_T
            fresh = model_from_parameters(fresh_params)
            @test ctx.model.τ_abs == fresh.τ_abs
            @test ctx.model.τ_rayl == fresh.τ_rayl
            @test ctx.model.profile.T == fresh.profile.T

            live_T = copy(ctx.model.profile.T)
            remembered_T = copy(ctx.current_T)
            retry_T .+= 50
            @test ctx.model.profile.T == live_T
            @test ctx.current_T == remembered_T
        end
    end
end

@testset "Aerosol updates commit all bands together" begin
    params = transaction_test_parameters(-1)
    template = parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
    params.scattering_params = deepcopy(template.scattering_params)
    params.scattering_params.r_max = 1.0
    params.scattering_params.nquad_radius = 20
    ctx = BatchContext(params)
    original = params.scattering_params.rt_aerosols[1]
    replacement = Aerosol(LogNormal(log(0.12), 0.4), 1.4, 0.001)
    snapshot() = (τ=map(copy, ctx.model.τ_aer), optics=map(copy, ctx.model.aerosol_optics),
        k=copy(ctx.k_ref), loading=copy(ctx.current_τ_ref), dist=copy(ctx.current_profile_dist),
        m=copy(ctx.model.solver.m_max_bands), l=copy(ctx.model.solver.l_max),
        n=copy(ctx.model.solver.n_fourier_moments_bands))
    before = snapshot()
    @test_throws MethodError update_aerosol_loading!(ctx, 1;
        τ_ref=0.3, profile_dist=:invalid_distribution)
    @test snapshot() == before

    # Inject a late staging failure: the first band's trial work succeeds,
    # then the second scratch slice cannot accept its candidate. Live arrays,
    # optical objects, loading and solver bounds must all remain unchanged.
    scratch = ctx.scratch_τ_aer[2]
    ctx.scratch_τ_aer[2] = zeros(1, 2, ctx.Nz)
    @test_throws DimensionMismatch update_aerosol_loading!(ctx, 1; τ_ref=0.3)
    @test snapshot() == before
    @test_throws DimensionMismatch update_aerosol_microphysics!(ctx, 1, replacement; τ_ref=0.3)
    @test snapshot() == before
    ctx.scratch_τ_aer[2] = scratch

    update_aerosol_loading!(ctx, 1; τ_ref=0.3)
    update_aerosol_microphysics!(ctx, 1, replacement)
    fresh_params = copy_parameters(params; share_luts=true)
    fresh_params.scattering_params.rt_aerosols[1].aerosol = replacement
    fresh_params.scattering_params.rt_aerosols[1].τ_ref = 0.3
    fresh = model_from_parameters(fresh_params)
    @test ctx.current_τ_ref == [0.3]
    @test ctx.model.τ_aer == fresh.τ_aer
    @test ctx.model.solver.m_max_bands == fresh.solver.m_max_bands
    @test ctx.model.solver.l_max == fresh.solver.l_max
    for ib in 1:ctx.n_bands
        @test rt_run(ctx.model; i_band=ib).toa ≈ rt_run(fresh; i_band=ib).toa rtol=2e-12
    end
    @test original.τ_ref != 0.3
end
