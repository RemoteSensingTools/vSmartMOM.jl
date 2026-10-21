# Surface-pressure Jacobians under the two vertical-grid conventions.
#
# FixedPressureGrid (default): only the bottom interface p_half[end] = p_surf
# moves. FixedAltitudeGrid: every interface is fixed in geometric altitude, so
# every interface pressure scales with p_surf. In both cases the analytic
# tangent must equal a centred finite difference of the forward model rebuilt
# with `surface_pressure_interfaces(convention, p_half, p_surf ± h)`.
#
# Run from test/:  julia --project=. test_vertical_grid_psurf.jl
using Test
using LinearAlgebra
using Distributions
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Absorption
import AtmosphericAbsorption

const _VG_GRIDS = (FixedPressureGrid(), FixedAltitudeGrid())

include("vertical_grid_psurf_fixture.jl")
include(joinpath("reference", "generate_vertical_grid_golden.jl"))  # read_golden


# A plan selecting surface pressure must see the requested convention.
struct _VGPressureFlavor <: AbstractJacobianFlavor end
CoreRT.requires_aerosol_microphysics_jacobians(::_VGPressureFlavor) = false
function CoreRT.jacobian_plan(flavor::_VGPressureFlavor, params, model, lin)
    keys = [ParameterKey(:atmosphere, :surface_pressure),
            ParameterKey(:surface, :P0; component=1, band=1)]
    names = ["psurf", "albedo"]
    layout = ActiveParameterLayout(keys, names, keys, names, [1], 1;
                                   surface_columns=2:2)
    JacobianPlan(flavor, keys, names, [layout])
end

@testset "vertical-grid surface-pressure Jacobians" begin

@testset "surface_pressure_interfaces" begin
    for FT in (Float32, Float64)
        p = FT[0.14, 50, 250, 500, 750, 1000]
        for grid in _VG_GRIDS
            same = surface_pressure_interfaces(grid, p, p[end])
            @test same == p && same !== p
            @test reinterpret(UInt8, same) == reinterpret(UInt8, p)
            @test eltype(same) === FT
        end
        fixed = surface_pressure_interfaces(FixedPressureGrid(), p, FT(990))
        @test fixed[1:end-1] == p[1:end-1] && fixed[end] == FT(990)
        @test_throws ArgumentError surface_pressure_interfaces(
            FixedPressureGrid(), p, FT(740))
        scaled = surface_pressure_interfaces(FixedAltitudeGrid(), p, FT(970))
        @test scaled[end] == FT(970)
        @test scaled ./ p ≈ fill(FT(0.97), length(p)) rtol=4eps(FT)
        @test_throws ArgumentError surface_pressure_interfaces(
            FixedAltitudeGrid(), p, FT(0))
        @test_throws ArgumentError surface_pressure_interfaces(
            FixedAltitudeGrid(), p, FT(NaN))
    end
end

@testset "profile tangents ($(nameof(typeof(grid))))" for grid in _VG_GRIDS
    T = [230.0, 250.0, 270.0, 285.0]
    q = [0.002, 0.004, 0.006, 0.01]
    p = [50.0, 250.0, 500.0, 750.0, 1000.0]
    vmr = Dict{String,Any}("CO2" => 4e-4)
    profile_at(ps) = first(CoreRT.prepare_observer_profile(
        T, surface_pressure_interfaces(grid, p, ps), q, vmr, 0.0, -1))
    base = profile_at(1000.0)
    tangents = CoreRT.psurf_profile_tangents(grid, base)
    h = 0.05
    plus, minus = profile_at(1000.0 + h), profile_at(1000.0 - h)
    for (field, dot) in ((:vcd_dry, :vcd_dry_dot), (:vcd_h2o, :vcd_h2o_dot),
                         (:Δz, :Δz_dot), (:p_full, :p_full_dot))
        fd = (getfield(plus, field) - getfield(minus, field)) / (2h)
        @test getfield(tangents, dot) ≈ fd rtol=1e-7 atol=1e-9 * maximum(abs, getfield(base, field))
    end
    if grid isa FixedAltitudeGrid
        # Geometric altitudes are fixed: the scaled grid reproduces them.
        @test CoreRT.half_level_altitudes(plus) ≈ CoreRT.half_level_altitudes(base) rtol=1e-12
        @test all(!iszero, tangents.p_full_dot) && all(iszero, tangents.Δz_dot)
    else
        @test all(iszero, tangents.p_full_dot[1:end-1]) && tangents.p_full_dot[end] == 0.5
    end
end

@testset "optics tangents vs rebuild: $(nameof(typeof(grid))) $kind water=$water" for
        grid in _VG_GRIDS, kind in (:line, :absco, :regular), water in (false, true)
    FT = Float64
    aerosols = kind === :absco && !water
    params = _vg_parameters(FT, kind; water, aerosols)
    model, lin = _vg_lin(params, grid)
    h = 0.05
    # Gas absorption: the cross-section response dσ/dp of every layer that
    # moves, plus the molecular-column response.
    fd_abs = _vg_fd(m -> m.τ_abs[1], params, grid, h)
    @test _vg_relerr(lin.τ̇_abs_psurf[1], fd_abs) < 2e-6
    if grid isa FixedAltitudeGrid
        # Every layer responds, and each layer separately matches.
        for iz in axes(fd_abs, 2)
            @test norm(fd_abs[:, iz]) > 0
            @test _vg_relerr(lin.τ̇_abs_psurf[1][:, iz], fd_abs[:, iz]) < 2e-6
        end
    else
        @test all(iszero, lin.τ̇_abs_psurf[1][:, 1:end-1])
    end
    fd_rayl = _vg_fd(m -> m.τ_rayl[1], params, grid, h)
    if grid isa FixedAltitudeGrid
        @test _vg_relerr(lin.τ̇_rayl_psurf[1], fd_rayl) < 1e-8
    else
        # The forward depolarization uses the column-weighted temperature,
        # which shifts slightly with the bottom-layer weight; the existing
        # FixedPressureGrid tangent holds it fixed (relative effect ~2e-6).
        # Under FixedAltitudeGrid all weights scale together, so it is exact.
        @test _vg_relerr(lin.τ̇_rayl_psurf[1], fd_rayl) < 1e-5
    end
    if aerosols
        fd_aer = _vg_fd(m -> m.τ_aer[1], params, grid, h)
        # Normal (pressure-form) aerosol: genuine redistribution.
        @test _vg_relerr(lin.τ̇_aer_psurf[1][2, :, :], fd_aer[2, :, :]) < 1e-6
        if grid isa FixedAltitudeGrid
            # LogNormal (altitude-form) aerosol is fixed in altitude.
            @test all(iszero, lin.τ̇_aer_psurf[1][1, :, :])
            @test maximum(abs, fd_aer[1, :, :]) < 1e-10 * maximum(model.τ_aer[1][1, :, :])
        else
            @test _vg_relerr(lin.τ̇_aer_psurf[1][1, :, :], fd_aer[1, :, :]) < 1e-6
        end
    end
end

@testset "CIA and MT_CKD tangents ($(nameof(typeof(grid))))" for grid in _VG_GRIDS
    # Through the constructor: a synthetic O2-O2 CIA file in the band.
    mktempdir() do directory
        params = _vg_parameters(Float64, :absco; cia_path=_vg_write_cia(
            joinpath(directory, "O2-O2.cia")))
        no_cia = _vg_parameters(Float64, :absco)
        model, lin = _vg_lin(params, grid)
        _, lin_no_cia = _vg_lin(no_cia, grid)
        h = 0.05
        fd_total = _vg_fd(m -> m.τ_abs[1], params, grid, h)
        fd_line = _vg_fd(m -> m.τ_abs[1], no_cia, grid, h)
        cia_tangent = lin.τ̇_abs_psurf[1] - lin_no_cia.τ̇_abs_psurf[1]
        @test norm(fd_total - fd_line) > 0
        @test _vg_relerr(cia_tangent, fd_total - fd_line) < 1e-6
        @test _vg_relerr(lin.τ̇_abs_psurf[1], fd_total) < 2e-6
    end
    # MT_CKD through the constructor's helper on the actual continuum code.
    T = [230.0, 250.0, 270.0, 285.0]
    q = [0.002, 0.004, 0.006, 0.01]
    p = [50.0, 250.0, 500.0, 750.0, 1000.0]
    profile_at(ps) = first(CoreRT.prepare_observer_profile(
        T, surface_pressure_interfaces(grid, p, ps), q, Dict{String,Any}(), 0.0, -1))
    band = Absorption.MTCKDBand([1e-27], [2e-27], [0.5], 1013.0, 296.0)
    ν = [5000.0]
    function continuum(profile)
        τ = zeros(1, length(profile.T))
        Absorption.compute_τ_h2o_continuum!(τ, band, ν, profile, profile.vmr_h2o)
        return τ
    end
    base = profile_at(1000.0)
    h = 0.01
    fd = (continuum(profile_at(1000.0 + h)) - continuum(profile_at(1000.0 - h))) / (2h)
    analytic = zeros(1, length(base.T))
    CoreRT._add_binary_psurf_tangent!(analytic, continuum(base), grid,
        CoreRT.psurf_profile_tangents(grid, base), base)
    @test analytic ≈ fd rtol=1e-8
end

@testset "radiance Jacobian vs FD: $(nameof(typeof(grid))) $kind external=$external" for
        grid in _VG_GRIDS, kind in (:line, :absco), external in (false, true)
    params = _vg_parameters(Float64, kind; water=true, aerosols=true)
    model, lin = _vg_lin(params, grid; external_solar=external)
    analytic = rt_run(model, lin, 2, size(lin.τ̇_abs[1], 1), 1)
    h = 0.05
    forward(m) = external ? rt_run_toa(m) : rt_run(m).toa
    fd = _vg_fd(forward, params, grid, h; external_solar=external)
    pressure_column = analytic.toa_jacobian[:, :, :, 1]
    @test norm(fd) > 0
    @test pressure_column ≈ fd rtol=5e-5
end

@testset "reference state is convention independent" begin
    params = _vg_parameters(Float64, :absco; water=true, aerosols=true)
    m_p, l_p = _vg_lin(params, FixedPressureGrid())
    m_z, l_z = _vg_lin(params, FixedAltitudeGrid())
    @test m_p.τ_abs == m_z.τ_abs && m_p.τ_rayl == m_z.τ_rayl && m_p.τ_aer == m_z.τ_aer
    @test l_p.τ̇_abs == l_z.τ̇_abs && l_p.τ̇_aer == l_z.τ̇_aer
    @test l_p.τ̇_abs_psurf != l_z.τ̇_abs_psurf
    r_p = rt_run(m_p, l_p, 2, size(l_p.τ̇_abs[1], 1), 1)
    r_z = rt_run(m_z, l_z, 2, size(l_z.τ̇_abs[1], 1), 1)
    @test r_p.toa == r_z.toa
    @test r_p.toa_jacobian[:, :, :, 2:end] == r_z.toa_jacobian[:, :, :, 2:end]
end

@testset "default convention is FixedPressureGrid, bit for bit" begin
    # Golden values written by the UNMODIFIED 7acab570 code (see
    # reference/generate_vertical_grid_golden.jl; Julia 1.12.5, x86-64 curry).
    # On the generating platform the default convention reproduces them bit
    # for bit; set VSMARTMOM_EXACT_GOLDEN=1 to require that. Elsewhere (other
    # Julia versions or CPUs: Mie/Voigt arithmetic, SIMD reductions) allow
    # last-digit differences.
    exact = get(ENV, "VSMARTMOM_EXACT_GOLDEN", "0") == "1"
    same(a, b) = exact ? a == b : isapprox(a, b; rtol=1e-12)
    golden = read_golden(joinpath(@__DIR__, "reference",
                                  "vertical_grid_fixed_pressure_7acab570.txt"))
    for kind in (:line, :absco)
        params = _vg_parameters(Float64, kind; water=true, aerosols=true)
        # 7acab570 derived the automatic (depol: -1) Rayleigh depolarization
        # from the constant Raman-module anisotropy; v2.2 defaults to the
        # dispersive anisotropy (decision 13c), which moves τ_rayl at
        # 1613 nm by -0.4 %. Select the formulation of the golden values.
        params.rayleigh_anisotropy = ConstantAnisotropy()
        model, implicit = model_from_parameters(LinMode(), params;
            compute_aerosol_microphysics_jacobians=false,
            abundance_jacobians=:frozen_cross_section)
        _, explicit = _vg_lin(params, FixedPressureGrid())
        for name in (:τ̇_abs_psurf, :τ̇_rayl_psurf, :τ̇_aer_psurf, :τ̇_abs, :τ̇_aer)
            @test getfield(implicit, name) == getfield(explicit, name)
        end
        @test same(implicit.τ̇_abs_psurf[1], golden["$(kind)_tau_abs_psurf"])
        @test same(implicit.τ̇_rayl_psurf[1], golden["$(kind)_tau_rayl_psurf"])
        @test same(implicit.τ̇_aer_psurf[1], golden["$(kind)_tau_aer_psurf"])
        if kind === :absco
            # 7acab570 had no ABSCO H₂O-broadener term in the q-driven H₂O
            # columns; removing exactly that term must restore its Jacobian.
            profile, grid = model.profile, params.spec_bands[1]
            broadener = zero(implicit.τ̇_abs[1])
            absorber = params.absorption_params.luts[1][1]
            for gas in ("CO2", "CH4")
                CoreRT._accumulate_broadener_tangent!(
                    broadener, 1, absorber, grid, profile.vmr[gas], profile)
            end
            r = profile.vmr_h2o
            CoreRT._accumulate_broadener_tangent!(
                broadener, 1, params.absorption_params.h2o_lut[1], grid, r, profile;
                broadener=r ./ (1 .+ r), dbroadener_dr=inv.((1 .+ r) .^ 2))
            @test any(!iszero, broadener)
            implicit = deepcopy(implicit)
            implicit.τ̇_abs[1] .-= broadener
        end
        result = rt_run(model, implicit, 2, size(implicit.τ̇_abs[1], 1), 1)
        @test Array(result.toa) ≈ golden["$(kind)_toa"] rtol=1e-13
        @test Array(result.toa_jacobian) ≈ golden["$(kind)_toa_jacobian"] rtol=1e-12
    end
end

@testset "Float32 FixedAltitudeGrid tangents" begin
    params = _vg_parameters(Float32, :absco; water=true)
    model, lin = _vg_lin(params, FixedAltitudeGrid())
    @test eltype(lin.τ̇_abs_psurf[1]) === Float32
    h = 0.5f0
    fd = _vg_fd(m -> m.τ_abs[1], params, FixedAltitudeGrid(), h)
    @test _vg_relerr(lin.τ̇_abs_psurf[1], fd) < 2e-3
    @test all(isfinite, lin.τ̇_rayl_psurf[1])
end

@testset "planned constructor forwards the convention" begin
    params = _vg_parameters(Float64, :absco)
    for grid in _VG_GRIDS
        model, planned = model_from_parameters(_VGPressureFlavor(), params;
            vertical_grid=grid)
        _, lin = _vg_lin(params, grid)
        @test planned.τ̇_abs_psurf == lin.τ̇_abs_psurf
        result = rt_run(model, planned; jacobian_basis=:physical)
        fd = _vg_fd(m -> rt_run(m).toa, params, grid, 0.05)
        @test result.toa_jacobian[:, :, :, 1] ≈ fd rtol=5e-5
    end
end

end # outer testset
