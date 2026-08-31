using Test
using vSmartMOM
using vSmartMOM.CoreRT
using Distributions: LogNormal

# Retrieval-selective Jacobian regression hierarchy:
#   1. Independent central finite differences for every active physical class
#      used by OCO_RRS_synth: pressure, aerosol tau_ref/z0, layer absorption,
#      three Legendre coefficients, SIF760, and mSIF.
#   2. Exact 16-layer OCO global/local layout and zero-fill behavior.
#   3. Bit-for-bit equality between compact propagation and corresponding
#      columns of the historical full native Jacobian.
#   4. Direct assertions that the Fourier-invariant aerosol/gas caches carry
#      compact derivative dimensions before phase mixing and MOM allocation.

struct _MockOCOJacobianModel
    profile
    spec_bands
    surfaces
end

@testset "Float32 normalized column tangent is overflow-safe" begin
    column = Float32[1.2e24, 3.4e24, 8.1e24, 1.7e25]
    column_dot = Float32[0, 0, 0, 1.3e22]
    fraction, fraction_dot = CoreRT._normalized_column_fraction_tangent(
        column, column_dot)
    @test all(isfinite, fraction)
    @test all(isfinite, fraction_dot)
    @test sum(fraction) ≈ 1f0 rtol=2f-7
    @test sum(fraction_dot) ≈ 0f0 atol=2f-10

    reference_column = Float64.(column)
    reference_dot = Float64.(column_dot)
    total = sum(reference_column)
    reference = (reference_dot .* total .-
                 reference_column .* sum(reference_dot)) ./ total^2
    @test Float64.(fraction_dot) ≈ reference rtol=5e-6 atol=1e-12
end

function _selective_fd_metrics(analytic, finite_difference; threshold=1e-10)
    absolute = abs.(analytic .- finite_difference)
    mask = abs.(finite_difference) .> threshold
    relative = any(mask) ? maximum(absolute[mask] ./ abs.(finite_difference[mask])) : 0.0
    return (max_abs=maximum(absolute), max_rel=relative)
end

function _test_selective_fd(label, analytic, finite_difference;
                            rtol, atol, threshold=1e-10)
    metrics = _selective_fd_metrics(analytic, finite_difference; threshold)
    @info "selective Jacobian finite difference" parameter=label metrics...
    @test analytic ≈ finite_difference rtol=rtol atol=atol
end

vSmartMOM.CoreRT.get_spec_bands(model::_MockOCOJacobianModel) = model.spec_bands
vSmartMOM.CoreRT.get_surface(model::_MockOCOJacobianModel, i_band) =
    model.surfaces[i_band]
vSmartMOM.CoreRT.n_aerosols(::_MockOCOJacobianModel) = 3

function _mock_oco_profile()
    FT = Float64
    Nz = 16
    # Four coarse upper layers have centers above 10 km; the remaining twelve
    # reproduce the active/fixed split of the synthetic OCO retrieval.
    Δz = FT.([20_000, 10_000, 5_000, 5_000, fill(800, 12)...])
    p_half = collect(range(FT(0.1), FT(1000), length=Nz + 1))
    p_full = (p_half[1:end-1] .+ p_half[2:end]) ./ 2
    vmr = Dict{String,Union{Real,Vector}}("CO2" => fill(400e-6, Nz))
    return CoreRT.AtmosphericProfile(fill(250.0, Nz), p_full, zeros(Nz), p_half,
        zeros(Nz), ones(Nz), zeros(Nz), vmr, Δz)
end


@testset "Selective aerosol/gas Jacobians vs central finite differences" begin
    params = parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
    params.architecture = vSmartMOM.Architectures.CPU()
    z₀ = 2.0
    σ₀ = 0.4
    aerosol = only(params.scattering_params.rt_aerosols)
    aerosol.profile = LogNormal(log(z₀), σ₀)

    model, lin_model = model_from_parameters(LinMode(), params;
        compute_aerosol_microphysics_jacobians=false,
        compute_h2o_jacobians=false)
    Nz = length(model.profile.p_full)
    ngas = size(lin_model.τ̇_abs[1], 1)
    native = ParameterLayout(n_aerosols=1, n_gases=ngas, n_surface=1)
    iz = Nz
    gas_native_column = gas_layer_index(native, 1, iz, Nz)

    # A synthetic absorption scale isolates the selected layer-gas column
    # without requiring an external CO2 LUT. Existing absorption remains in
    # the forward state; x changes only the bottom-layer additive opacity.
    κ = collect(range(0.01, 0.02; length=size(model.τ_abs[1], 1)))
    lin_model.τ̇_abs[1][iz, :, iz] .= κ

    keys = [
        ParameterKey(:atmosphere, :surface_pressure),
        ParameterKey(:aerosol, :tau_ref; component=1),
        ParameterKey(:aerosol, :z0; component=1),
        ParameterKey(:gas, :vmr; component="synthetic", layer=iz),
        ParameterKey(:surface, :albedo; component=1, band=1),
    ]
    names = ["surface_pressure", "aerosol_tau_ref", "aerosol_z0",
             "synthetic_gas_layer$(iz)", "surface_albedo"]
    layout = ActiveParameterLayout(
        keys, names, keys, names,
        [1, aerosol_range(native, 1)[1], aerosol_range(native, 1)[6],
         gas_native_column], 4;
        surface_columns=5:5)
    plan = JacobianPlan(OCO_RRS_synth(), keys, names, [layout])
    selected = rt_run(model, PlannedRTModelLin(lin_model, plan); i_band=1)

    # BatchContext changes aerosol loading/profile while reusing identical
    # forward Mie optics, keeping these central differences independent of the
    # linearized optical-property builder and inexpensive enough for CI.
    context = BatchContext(params)
    τₒ = aerosol.τ_ref
    hτ = τₒ * 1e-3
    update_aerosol_loading!(context, 1;
        τ_ref=τₒ + hτ, profile_dist=LogNormal(log(z₀), σ₀))
    Rτp = rt_run_toa(context.model)
    update_aerosol_loading!(context, 1;
        τ_ref=τₒ - hτ, profile_dist=LogNormal(log(z₀), σ₀))
    Rτm = rt_run_toa(context.model)
    fd_τ = (Rτp .- Rτm) ./ (2hτ)
    _test_selective_fd("tau_ref", selected.toa_jacobian[:, :, :, 2], fd_τ;
                       rtol=8e-3, atol=2e-8, threshold=1e-8)

    hz = 1e-3
    update_aerosol_loading!(context, 1;
        τ_ref=τₒ, profile_dist=LogNormal(log(z₀ + hz), σ₀))
    Rzp = rt_run_toa(context.model)
    update_aerosol_loading!(context, 1;
        τ_ref=τₒ, profile_dist=LogNormal(log(z₀ - hz), σ₀))
    Rzm = rt_run_toa(context.model)
    fd_z = (Rzp .- Rzm) ./ (2hz)
    _test_selective_fd("z0", selected.toa_jacobian[:, :, :, 3], fd_z;
                       rtol=1e-2, atol=2e-8, threshold=1e-9)

    hx = 1e-5
    plus = deepcopy(model)
    minus = deepcopy(model)
    plus.τ_abs[1][:, iz] .+= hx .* κ
    minus.τ_abs[1][:, iz] .-= hx .* κ
    fd_gas = (rt_run_toa(plus) .- rt_run_toa(minus)) ./ (2hx)
    _test_selective_fd("layer absorption", selected.toa_jacobian[:, :, :, 4], fd_gas;
                       rtol=5e-4, atol=2e-8, threshold=1e-8)
end


@testset "Selective pressure/surface/SIF Jacobians vs central finite differences" begin
    params = parameters_from_yaml("test_parameters/JacobianTestRayleigh.yaml")
    params.architecture = vSmartMOM.Architectures.CPU()
    coeff = [0.24, 0.015, -0.004]
    params.brdf = CoreRT.AbstractSurfaceType[
        CoreRT.LambertianSurfaceLegendre(coeff)]
    model, lin_model = model_from_parameters(LinMode(), params;
        compute_aerosol_microphysics_jacobians=false,
        compute_h2o_jacobians=false)

    keys = [
        ParameterKey(:atmosphere, :surface_pressure),
        ParameterKey(:surface, :P0; component=1, band=1),
        ParameterKey(:surface, :P1; component=1, band=1),
        ParameterKey(:surface, :P2; component=1, band=1),
        ParameterKey(:source, :SIF760; component="SIF", band=1),
        ParameterKey(:source, :mSIF; component="SIF", band=1),
    ]
    names = ["surface_pressure", "surface_P0", "surface_P1", "surface_P2",
             "SIF760", "mSIF"]
    layout = ActiveParameterLayout(keys, names, keys, names, [1], 1;
                                   surface_columns=2:4, sif_columns=5:6)
    plan = JacobianPlan(OCO_RRS_synth(), keys, names, [layout])

    ν = collect(model.atmosphere.spec_bands[1])
    F₀ = zeros(Float64, CoreRT.polarization_type(model).n, length(ν))
    F₀[1, :] .= collect(range(1.9, 2.1; length=length(ν)))
    sif760 = 0.03
    msif = 1e-4
    make_sources(s, m) = CoreRT.SolarBeam(F₀=F₀) +
        CoreRT.SurfaceSIF(SIF760=s, mSIF=m, wavenumber_cm1=ν)
    sources = make_sources(sif760, msif)
    selected = rt_run(model, PlannedRTModelLin(lin_model, plan);
                      i_band=1, sources)
    @test selected.toa ≈ rt_run_toa(model; sources) rtol=2e-12 atol=2e-14

    hp = 0.01
    pp = deepcopy(params); pp.p[end] += hp
    pm = deepcopy(params); pm.p[end] -= hp
    fd_pressure = (rt_run_toa(model_from_parameters(pp); sources) .-
                   rt_run_toa(model_from_parameters(pm); sources)) ./ (2hp)
    _test_selective_fd("surface pressure", selected.toa_jacobian[:, :, :, 1],
                       fd_pressure; rtol=2e-4, atol=2e-8, threshold=1e-8)

    hc = 1e-5
    for icoeff in eachindex(coeff)
        cp = copy(coeff); cp[icoeff] += hc
        cm = copy(coeff); cm[icoeff] -= hc
        pplus = deepcopy(params)
        pminus = deepcopy(params)
        pplus.brdf = CoreRT.AbstractSurfaceType[
            CoreRT.LambertianSurfaceLegendre(cp)]
        pminus.brdf = CoreRT.AbstractSurfaceType[
            CoreRT.LambertianSurfaceLegendre(cm)]
        fd_surface = (rt_run_toa(model_from_parameters(pplus); sources) .-
                      rt_run_toa(model_from_parameters(pminus); sources)) ./ (2hc)
        _test_selective_fd("surface P$(icoeff - 1)",
                           selected.toa_jacobian[:, :, :, 1 + icoeff],
                           fd_surface; rtol=3e-4, atol=2e-8, threshold=1e-8)
    end

    hsif = 1e-6
    fd_sif = (rt_run_toa(model; sources=make_sources(sif760 + hsif, msif)) .-
              rt_run_toa(model; sources=make_sources(sif760 - hsif, msif))) ./ (2hsif)
    _test_selective_fd("SIF760", selected.toa_jacobian[:, :, :, 5], fd_sif;
                       rtol=2e-8, atol=2e-9, threshold=1e-8)

    hslope = 1e-8
    fd_slope = (rt_run_toa(model; sources=make_sources(sif760, msif + hslope)) .-
                rt_run_toa(model; sources=make_sources(sif760, msif - hslope))) ./ (2hslope)
    _test_selective_fd("mSIF", selected.toa_jacobian[:, :, :, 6], fd_slope;
                       rtol=2e-7, atol=2e-8, threshold=1e-8)
end

@testset "OCO_RRS_synth plan" begin
    profile = _mock_oco_profile()
    surfaces = [CoreRT.LambertianSurfaceLegendre([0.2, 0.0, 0.0]) for _ in 1:3]
    model = _MockOCOJacobianModel(
        profile, [Float64[1.0] for _ in 1:3], surfaces)
    params = (absorption_params=(
        variable_molecules=[String[], ["CO2"], ["CO2"]],),)
    lin_model = (τ̇_abs=[zeros(32, 1, 16) for _ in 1:3],)

    plan = jacobian_plan(OCO_RRS_synth(), params, model, lin_model)
    @test n_global(plan) == 30
    @test parameter_names(plan)[1] == "surface_pressure"
    @test parameter_names(plan)[2] == "co2_vmr_layer05"
    @test parameter_names(plan)[13] == "co2_vmr_layer16"
    @test parameter_names(plan)[14] == "sulfate_tau_ref"
    @test parameter_names(plan)[17] == "sulfate_z0"
    @test parameter_names(plan)[29:30] == ["SIF760", "mSIF"]
    @test n_total(band_layout(plan, 1)) == 12
    @test n_total(band_layout(plan, 2)) == 22
    @test n_total(band_layout(plan, 3)) == 22
    @test n_layer_params(band_layout(plan, 1)) == 7
    @test n_layer_params(band_layout(plan, 2)) == 19
    @test length(surface_range(band_layout(plan, 1))) == 3
    @test length(sif_range(band_layout(plan, 1))) == 2
    @test isempty(sif_range(band_layout(plan, 2)))

    local_jac = reshape(collect(1.0:24.0), 2, 12)
    global_jac = globalize_jacobian(local_jac, band_layout(plan, 1))
    @test size(global_jac) == (2, 30)
    @test global_jac[:, local_to_global(band_layout(plan, 1))] == local_jac
    unused = setdiff(1:30, local_to_global(band_layout(plan, 1)))
    @test all(iszero, global_jac[:, unused])
end

@testset "Compact tangent propagation equals selected full columns" begin
    params = parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
    model_full, lin_full = model_from_parameters(LinMode(), params)
    model_fixed, lin_fixed = model_from_parameters(LinMode(), params;
        compute_aerosol_microphysics_jacobians=false,
        compute_h2o_jacobians=false)

    @test model_fixed.τ_aer[1] == model_full.τ_aer[1]
    @test all(iszero, lin_fixed.lin_aerosol_optics[1][1].k̇)

    ngas = size(lin_full.τ̇_abs[1], 1)
    full_layout = ParameterLayout(n_aerosols=1, n_gases=ngas, n_surface=1)
    full = rt_run(model_full, lin_full, 1, ngas, 1)

    keys = [
        ParameterKey(:atmosphere, :surface_pressure),
        ParameterKey(:aerosol, :tau_ref; component=1),
        ParameterKey(:aerosol, :profile_location; component=1),
        ParameterKey(:surface, :albedo; component=1, band=1),
    ]
    names = ["surface_pressure", "aerosol_tau_ref",
             "aerosol_profile_location", "surface_albedo"]
    layout = ActiveParameterLayout(keys, names, keys, names, [1, 2, 7], 3;
                                   surface_columns=4:4)
    plan = JacobianPlan(OCO_RRS_synth(), keys, names, [layout])
    cache = CoreRT.build_m_invariant_cache_lin(
        1, model_fixed, lin_fixed; active_layout=layout)
    @test size(cache.aerosol[1][1][1].τ̇, 2) == 2
    @test size(cache.lin_gas[1][1].τ̇, 2) == 0
    selected = rt_run(model_fixed, PlannedRTModelLin(lin_fixed, plan); i_band=1)
    native_columns = [1, 2, 7, n_total(full_layout)]

    @test selected.toa == full.toa
    @test selected.boa == full.boa
    @test size(selected.toa_jacobian, 4) == 4
    @test selected.toa_jacobian == full.toa_jacobian[:, :, :, native_columns]
    @test selected.boa_jacobian == full.boa_jacobian[:, :, :, native_columns]
end
