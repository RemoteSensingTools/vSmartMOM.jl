using Test
using LinearAlgebra
using vSmartMOM
using vSmartMOM.CoreRT
import AtmosphericAbsorption
import NCDatasets

# Native gas columns are ∂τ/∂VMR = σ·N_dry with cross sections, CIA and
# MT_CKD held fixed, except that the q-driven H₂O columns include the exact
# response of ABSCO H₂O-broadener axes in both abundance modes. These tests
# show, with finite differences of the
# forward absorption, where the model's own forward response is larger than
# the native column, and that such columns are rejected unless explicitly
# accepted.
# Synthetic spectroscopy keeps the file independent of downloaded data.

const _AA = AtmosphericAbsorption
const _GRID = [6199.8, 6199.95, 6200.0, 6200.05, 6200.2]

# A strongly self-broadened water line: γ_self ≫ γ_air.
function _water_line()
    lines = _AA.LineDatabase(
        mol=Int32[1], iso=Int32[1], ν0=[6200.0], S=[1e-25], E_lower=[0.0],
        g_upper=[1.0], γ_air=[0.07], γ_self=[0.30], n_air=[0.7], δ_air=[0.0],
        n_self=[0.5], δ_self=[0.0], Y_LM=[0.0], n_Y_LM=[0.5],
        γ2_air=[0.0], δ2_air=[0.0], νVC=[0.0], η=[0.0], molar_mass=[18.0],
        meta=_AA.SourceMetadata("synthetic abundance test", 296.0, 1013.25))
    partition = _AA.TabulatedPF([100.0, 300.0, 500.0], [1.0, 1.0, 1.0])
    return _AA.LineByLineModel(lines, partition; profile=_AA.Voigt(),
        cpf=_AA.HumlicekWeideman32(), wing_cutoff=10.0, vmr=0.0,
        architecture=_AA.CPU())
end

# A table with a fixed broadener abundance.
function _fixed_broadener_table(scale=1e-23)
    ν, p, T = [6199.8, 6200, 6200.2], [100.0, 600, 1100], [210.0, 290, 330]
    σ = [scale * i * (1 + 0.001pp + 0.004tt) for i in eachindex(ν), pp in p, tt in T]
    return _AA.InterpolationModel(ν, p, T, σ, 0.0, _AA.CPU())
end

# An ABSCO table whose cross section grows with its H₂O-broadener node; with
# `broadener=[0.02]` the single node fixes the broadener abundance. A nonzero
# `curvature` gives the broadener intervals different slopes.
function _absco_table(; broadener=[0.0, 0.05], mol=2, curvature=0.0, FT=Float64)
    ν, p = FT[6199.8, 6200, 6200.2], FT[100.0, 600, 1100]
    T = FT[210.0 230 250; 290 310 330]
    broadener = FT.(broadener)
    σ = zeros(FT, length(ν), length(broadener), size(T, 1), length(p))
    for ip in eachindex(p), it in axes(T, 1), ib in eachindex(broadener), iν in eachindex(ν)
        b = broadener[ib]
        σ[iν, ib, it, ip] = 1e-23 * iν *
            (1 + 0.001p[ip] + 0.004T[it, ip] + 20b + curvature * b^2)
    end
    return _AA.AbscoLUT(mol, -1, ν, p, T, broadener, σ; architecture=_AA.CPU())
end

function _write_cia(path, pair, σ)
    header = rpad(pair, 20) * lpad("6199.000", 10) * lpad("6201.000", 10) *
             lpad("2", 7) * lpad("250.0", 7) * lpad("1.000e-40", 10)
    open(path, "w") do io
        println(io, header)
        println(io, "6199.000 $σ")
        println(io, "6201.000 $σ")
    end
    return path
end

function _write_mtckd(path)
    NCDatasets.NCDataset(path, "c") do ds
        ds.attrib["Title"] = "synthetic MT_CKD"
        ds.attrib["Version_description"] = "abundance contract test"
        NCDatasets.defDim(ds, "nwave", 2)
        NCDatasets.defVar(ds, "wavenumbers", [6000.0, 6400.0], ("nwave",))
        NCDatasets.defVar(ds, "self_absco_ref", [2e-22, 2e-22], ("nwave",))
        NCDatasets.defVar(ds, "for_absco_ref", [1e-23, 1e-23], ("nwave",))
        NCDatasets.defVar(ds, "self_texp", [4.0, 4.0], ("nwave",))
        NCDatasets.defVar(ds, "ref_press", 1013.0, ())
        NCDatasets.defVar(ds, "ref_temp", 296.0, ())
    end
    return path
end

"Two-layer clear-sky scene; `h2o` is the q-driven H₂O absorber or `:disabled`."
function _scene(; h2o=:disabled, q=0.02, fixed=String[], variable=String[],
                luts=Any[], vmr=Dict{String,Any}(), cia_files=String[], mtckd="",
                float_type="Float64")
    params = read_parameters(Dict(
        "radiative_transfer" => Dict(
            "spec_bands" => ["[6199.8, 6199.95, 6200.0, 6200.05, 6200.2]"],
            "surface" => ["LambertianSurfaceScalar(0.2)"],
            "nstreams" => 3, "polarization_type" => "Stokes_I()",
            "truncation" => "NoTruncation()", "depol" => -1,
            "float_type" => float_type, "architecture" => "CPU()"),
        "geometry" => Dict("sza" => 35.0, "vza" => [25.0], "vaz" => [40.0], "obs_alt" => [0]),
        "atmospheric_profile" => Dict("T" => [270.0, 280.0],
            "p" => [100.0, 500.0, 1000.0], "profile_reduction" => -1)))
    params.absorption_params = CoreRT.AbsorptionParameters(
        [fixed], [variable], vmr, _AA.Voigt(), _AA.HumlicekWeideman32(), 10.0,
        [luts], Any[h2o], cia_files, mtckd)
    params.q .= q
    return params
end

_full_solve(model, lin) = rt_run(model, lin, 0, size(lin.τ̇_abs[1], 1), 1)

"Centered difference of one layer's optical depth in the dry-air H₂O ratio r."
function _fd_in_h2o_ratio(τ_of_profile, profile, iz; rel=1e-4)
    function τ_at(r)
        perturbed = deepcopy(profile)
        perturbed.vmr_h2o[iz] = r
        perturbed.vcd_h2o[iz] = r * perturbed.vcd_dry[iz]  # fixed dry column
        return τ_of_profile(perturbed)[:, iz]
    end
    r = profile.vmr_h2o[iz]
    h = rel * r
    return (τ_at(r + h) - τ_at(r - h)) / (2h)
end

_mismatch(native, fd) = norm(native - fd) / norm(fd)

@testset "Self-broadened H₂O abundance tangent" begin
    line = _water_line()
    params = _scene(h2o=line)
    model, lin = model_from_parameters(LinMode(), params)
    profile, Nz, iz = model.profile, 2, 2
    native = lin.τ̇_abs[1][iz, :, iz]   # q-driven H₂O occupies gas slot 1

    function h2o_τ(profile)
        τ = zeros(length(_GRID), Nz)
        CoreRT.compute_h2o_absorption_profile!(τ, line, _GRID, profile)
        return τ
    end
    @test h2o_τ(profile) ≈ model.τ_abs[1]
    fd = _fd_in_h2o_ratio(h2o_τ, profile, iz)
    # The native partial misses ∂σ/∂x_H₂O·x_H₂O; at x≈0.02 with γ_self/γ_air≈4
    # it overestimates the line-centre response by several per cent.
    @test _mismatch(native, fd) > 0.02
    @test native[3] / fd[3] > 1.05
    @test fd ≈ _fd_in_h2o_ratio(h2o_τ, profile, iz; rel=1e-3) rtol=1e-6

    err = try
        _full_solve(model, lin)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("q-driven H₂O", sprint(showerror, err))
    @test occursin("self broadening", sprint(showerror, err))
    @test occursin("abundance_jacobians=:frozen_cross_section", sprint(showerror, err))

    # Explicit acceptance of the native partial; its values are unchanged.
    model_p, lin_p = model_from_parameters(LinMode(), params;
        abundance_jacobians=:frozen_cross_section)
    @test lin_p.τ̇_abs == lin.τ̇_abs
    @test lin_p.abundance_contract.omitted == lin.abundance_contract.omitted
    result = _full_solve(model_p, lin_p)
    @test all(isfinite, result.toa_jacobian)
    @test lin_p.τ̇_abs[1][iz, :, iz] ≈ native

    # The alias forwards the keyword.
    _, lin_alias = model_from_parameters_lin(params;
        abundance_jacobians=:frozen_cross_section)
    @test lin_alias.abundance_contract.abundance_jacobians === :frozen_cross_section

    # Fixed humidity remains a supported workflow without any opt-in.
    model_f, lin_f = model_from_parameters(LinMode(), params;
        compute_h2o_jacobians=false)
    @test model_f.τ_abs ≈ model.τ_abs   # forward helper: a few ulps
    @test all(isfinite, _full_solve(model_f, lin_f).toa_jacobian)

    @test_throws ArgumentError model_from_parameters(LinMode(), params;
        abundance_jacobians=:complete_physical)
end

@testset "Fixed-broadener H₂O table tangent is complete" begin
    table = _fixed_broadener_table(1e-25)
    params = _scene(h2o=table)
    model, lin = model_from_parameters(LinMode(), params)
    @test all(isempty, only(lin.abundance_contract.omitted))
    function h2o_τ(profile)
        τ = zeros(length(_GRID), 2)
        CoreRT.compute_h2o_absorption_profile!(τ, table, _GRID, profile)
        return τ
    end
    fd = _fd_in_h2o_ratio(h2o_τ, model.profile, 2)
    @test lin.τ̇_abs[1][2, :, 2] ≈ fd rtol=1e-8
    @test all(isfinite, _full_solve(model, lin).toa_jacobian)
end

@testset "ABSCO H₂O-broadener axis" begin
    absco = _absco_table()
    params = _scene(fixed=["CO2"], luts=Any[absco], vmr=Dict{String,Any}("CO2" => 4e-4))
    model, lin = model_from_parameters(LinMode(), params)
    function co2_τ(profile)
        τ = zeros(length(_GRID), 2)
        CoreRT.compute_absorption_profile!(τ, absco, _GRID, 4e-4, profile)
        return τ
    end
    @test co2_τ(model.profile) ≈ model.τ_abs[1]
    # Forward CO₂ opacity responds to H₂O through the broadener axis, and the
    # q-driven H₂O column carries that response.
    fd = _fd_in_h2o_ratio(co2_τ, model.profile, 2)
    @test norm(fd) > 0
    @test lin.τ̇_abs[1][2, :, 2] ≈ fd rtol=1e-8
    @test all(isempty, only(lin.abundance_contract.omitted))
    @test all(isfinite, _full_solve(model, lin).toa_jacobian)
    # The exact term is kept under :frozen_cross_section too: same values,
    # nothing omitted.
    _, lin_p = model_from_parameters(LinMode(), params;
        abundance_jacobians=:frozen_cross_section)
    @test lin_p.τ̇_abs == lin.τ̇_abs
    @test lin_p.τ̇_abs_psurf == lin.τ̇_abs_psurf
    @test all(isempty, only(lin_p.abundance_contract.omitted))
    # Fixed humidity: no H₂O columns at all.
    _, lin_f = model_from_parameters(LinMode(), params; compute_h2o_jacobians=false)
    @test all(iszero, lin_f.τ̇_abs[1][1:2, :, :])
    @test all(isfinite, _full_solve(model, lin_f).toa_jacobian)

    # H₂O's own table (broadener x = r/(1+r)) and a variable CO₂ table, both
    # inside the second of three broadener intervals, whose slope differs
    # from the first.
    nodes = [0.0, 0.03, 0.06]
    h2o3 = _absco_table(broadener=nodes, mol=1, curvature=400.0)
    co23 = _absco_table(broadener=nodes, curvature=400.0)
    wet = _scene(h2o=h2o3, variable=["CO2"], luts=Any[co23],
                 vmr=Dict{String,Any}("CO2" => [4e-4, 4e-4]))
    model_w, lin_w = model_from_parameters(LinMode(), wet)
    profile_w = model_w.profile
    @test all(r -> nodes[2] < r / (1 + r) && r < nodes[3], profile_w.vmr_h2o)
    function wet_τ(profile)
        τ = zeros(length(_GRID), 2)
        CoreRT.compute_h2o_absorption_profile!(τ, h2o3, _GRID, profile)
        CoreRT.compute_absorption_profile!(τ, co23, _GRID, 4e-4, profile)
        return τ
    end
    @test wet_τ(profile_w) ≈ model_w.τ_abs[1]
    @test all(isempty, only(lin_w.abundance_contract.omitted))
    for iz in 1:2
        @test lin_w.τ̇_abs[1][iz, :, iz] ≈ _fd_in_h2o_ratio(wet_τ, profile_w, iz) rtol=1e-6
    end
    # The CO₂ column (gas slot 2) remains σ·N_dry.
    co2_σ = vSmartMOM.CoreRT._layer_absorption_cross_section(
        co23, _GRID, profile_w.p_full[2], profile_w.T[2], profile_w.vmr_h2o[2])
    @test lin_w.τ̇_abs[1][2 + 2, :, 2] ≈ co2_σ .* profile_w.vcd_dry[2]
    @test all(isfinite, _full_solve(model_w, lin_w).toa_jacobian)
    _, lin_wp = model_from_parameters(LinMode(), wet;
        abundance_jacobians=:frozen_cross_section)
    @test lin_wp.τ̇_abs == lin_w.τ̇_abs
    @test lin_wp.τ̇_abs_psurf == lin_w.τ̇_abs_psurf
    @test all(isempty, only(lin_wp.abundance_contract.omitted))

    # A single broadener node fixes the broadener abundance for H₂O and CO₂.
    single_h2o = _absco_table(broadener=[0.02], mol=1)
    single_co2 = _absco_table(broadener=[0.02])
    fixed_axis = _scene(h2o=single_h2o, fixed=["CO2"], luts=Any[single_co2],
                        vmr=Dict{String,Any}("CO2" => 4e-4))
    model_s, lin_s = model_from_parameters(LinMode(), fixed_axis)
    @test all(isempty, only(lin_s.abundance_contract.omitted))
    function single_τ(profile)
        τ = zeros(length(_GRID), 2)
        CoreRT.compute_h2o_absorption_profile!(τ, single_h2o, _GRID, profile)
        CoreRT.compute_absorption_profile!(τ, single_co2, _GRID, 4e-4, profile)
        return τ
    end
    @test single_τ(model_s.profile) ≈ model_s.τ_abs[1]
    @test lin_s.τ̇_abs[1][2, :, 2] ≈ _fd_in_h2o_ratio(single_τ, model_s.profile, 2) rtol=1e-8
    @test all(isfinite, _full_solve(model_s, lin_s).toa_jacobian)
end

@testset "ABSCO H₂O-broadener tangent: intervals, clamping, offsets, types" begin
    # The scene's uniform humidity gives one H₂O ratio r in both layers.
    r = first(model_from_parameters(LinMode(), _scene())).profile.vmr_h2o[1]
    x = r / (1 + r)

    # A node between x and r: H₂O's own table (broadener x) sits in the first
    # interval and the CO₂ table (broadener r) in the second, whose slope
    # differs, so a wrong broadener coordinate would select the wrong slope.
    straddle = [0.0, (x + r) / 2, 0.08]
    h2o_s = _absco_table(broadener=straddle, mol=1, curvature=900.0)
    co2_s = _absco_table(broadener=straddle, curvature=900.0)
    params = _scene(h2o=h2o_s, fixed=["CO2"], luts=Any[co2_s],
                    vmr=Dict{String,Any}("CO2" => 4e-4))
    model, lin = model_from_parameters(LinMode(), params)
    function straddle_τ(profile)
        τ = zeros(length(_GRID), 2)
        CoreRT.compute_h2o_absorption_profile!(τ, h2o_s, _GRID, profile)
        CoreRT.compute_absorption_profile!(τ, co2_s, _GRID, 4e-4, profile)
        return τ
    end
    @test straddle_τ(model.profile) ≈ model.τ_abs[1]
    for iz in 1:2
        @test lin.τ̇_abs[1][iz, :, iz] ≈ _fd_in_h2o_ratio(straddle_τ, model.profile, iz) rtol=1e-6
    end

    # Above the last node the forward query clamps, so the response is zero.
    clamped = _absco_table(broadener=[0.0, r / 2])
    params_c = _scene(fixed=["CO2"], luts=Any[clamped], vmr=Dict{String,Any}("CO2" => 4e-4))
    model_c, lin_c = model_from_parameters(LinMode(), params_c)
    function clamped_τ(profile)
        τ = zeros(length(_GRID), 2)
        CoreRT.compute_absorption_profile!(τ, clamped, _GRID, 4e-4, profile)
        return τ
    end
    @test all(iszero, _fd_in_h2o_ratio(clamped_τ, model_c.profile, 2))
    @test all(iszero, lin_c.τ̇_abs[1][1:2, :, :])

    # Fixed O₂, variable CO₂ and H₂O tables in one band: the variable table's
    # LUT offset, and every table's broadener term in the H₂O column.
    nodes = [0.0, 0.03, 0.06]
    o2_m = _absco_table(broadener=nodes, mol=7, curvature=-300.0)
    co2_m = _absco_table(broadener=nodes, curvature=400.0)
    h2o_m = _absco_table(broadener=nodes, mol=1, curvature=200.0)
    mixed = _scene(h2o=h2o_m, fixed=["O2"], variable=["CO2"], luts=Any[o2_m, co2_m],
                   vmr=Dict{String,Any}("O2" => 0.21, "CO2" => [4e-4, 4e-4]))
    model_m, lin_m = model_from_parameters(LinMode(), mixed)
    profile_m = model_m.profile
    function mixed_τ(profile)
        τ = zeros(length(_GRID), 2)
        CoreRT.compute_h2o_absorption_profile!(τ, h2o_m, _GRID, profile)
        CoreRT.compute_absorption_profile!(τ, o2_m, _GRID, 0.21, profile)
        CoreRT.compute_absorption_profile!(τ, co2_m, _GRID, 4e-4, profile)
        return τ
    end
    @test mixed_τ(profile_m) ≈ model_m.τ_abs[1]
    @test all(isempty, only(lin_m.abundance_contract.omitted))
    for iz in 1:2
        @test lin_m.τ̇_abs[1][iz, :, iz] ≈ _fd_in_h2o_ratio(mixed_τ, profile_m, iz) rtol=1e-6
        co2_σ = CoreRT._layer_absorption_cross_section(
            co2_m, _GRID, profile_m.p_full[iz], profile_m.T[iz], profile_m.vmr_h2o[iz])
        @test lin_m.τ̇_abs[1][2 + iz, :, iz] ≈ co2_σ .* profile_m.vcd_dry[iz]
    end
    @test all(isfinite, _full_solve(model_m, lin_m).toa_jacobian)

    # With q ≡ 0 the tables are queried at broadener 0 and the H₂O columns are
    # structural zeros.
    dry = _scene(q=0.0, fixed=["CO2"], luts=Any[co2_m], vmr=Dict{String,Any}("CO2" => 4e-4))
    model_d, lin_d = model_from_parameters(LinMode(), dry)
    for iz in 1:2
        σ0 = CoreRT._layer_absorption_cross_section(
            co2_m, _GRID, model_d.profile.p_full[iz], model_d.profile.T[iz], 0.0)
        @test model_d.τ_abs[1][:, iz] ≈ σ0 .* model_d.profile.vcd_dry[iz] .* 4e-4
    end
    @test all(iszero, lin_d.τ̇_abs[1][1:2, :, :])

    # A Float32 model with Float32 tables reproduces the Float64 tangent at
    # 6200 cm⁻¹, the one grid point (and table node) exact in Float32; the
    # others differ through their rounded spectral weights.
    mixed32 = _scene(h2o=_absco_table(broadener=nodes, mol=1, curvature=200.0, FT=Float32),
                     fixed=["O2"], variable=["CO2"],
                     luts=Any[_absco_table(broadener=nodes, mol=7, curvature=-300.0, FT=Float32),
                              _absco_table(broadener=nodes, curvature=400.0, FT=Float32)],
                     vmr=Dict{String,Any}("O2" => 0.21, "CO2" => [4e-4, 4e-4]),
                     float_type="Float32")
    _, lin_32 = model_from_parameters(LinMode(), mixed32)
    @test eltype(lin_32.τ̇_abs[1]) === Float32
    @test all(isfinite, lin_32.τ̇_abs[1])
    @test lin_32.τ̇_abs[1][:, 3, :] ≈ lin_m.τ̇_abs[1][:, 3, :] rtol=1e-5
end

@testset "ABSCO H₂O-broadener derivative at nodes and outside the table" begin
    # The derivative is the slope of the interval AtmosphericAbsorption
    # evaluates: the right one at the first and interior nodes, the left one at
    # the last node, and zero where the forward query clamps.
    nodes = [0.01, 0.03, 0.06]
    table = _absco_table(broadener=nodes, curvature=400.0)
    p, T = 600.0, 270.0
    σ_at(b) = CoreRT._layer_absorption_cross_section(table, _GRID, p, T, b)
    slope(i) = (σ_at(nodes[i+1]) .- σ_at(nodes[i])) ./ (nodes[i+1] - nodes[i])
    dσ(b) = CoreRT._absorption_broadener_derivative(table, _GRID, p, T, b)
    @test slope(1) != slope(2)
    @test dσ(nodes[1]) ≈ slope(1)
    @test dσ(nodes[2]) ≈ slope(2)
    @test dσ(nodes[3]) ≈ slope(2)
    @test dσ(0.02) ≈ slope(1)
    @test all(iszero, dσ(0.005))
    @test all(iszero, dσ(0.07))
    # The forward query is flat on both clamped sides.
    @test σ_at(0.005) == σ_at(nodes[1])
    @test σ_at(0.07) == σ_at(nodes[3])
end

@testset "MT_CKD continuum abundance response" begin
    mktempdir() do dir
        path = _write_mtckd(joinpath(dir, "mtckd.nc"))
        params = _scene(mtckd=path)
        model, lin = model_from_parameters(LinMode(), params)
        table = vSmartMOM.Absorption.load_mtckd(path)
        function continuum_τ(profile)
            τ = zeros(length(_GRID), 2)
            vSmartMOM.Absorption.compute_τ_h2o_continuum!(
                τ, table, _GRID, profile, profile.vmr_h2o)
            return τ
        end
        @test continuum_τ(model.profile) ≈ model.τ_abs[1]
        fd = _fd_in_h2o_ratio(continuum_τ, model.profile, 2)
        @test norm(fd) > 0
        @test all(iszero, lin.τ̇_abs[1][2, :, 2])
        @test_throws r"MT_CKD H₂O self and foreign continuum" _full_solve(model, lin)
    end
end

@testset "CIA collider abundance response" begin
    mktempdir() do dir
        o2o2 = _write_cia(joinpath(dir, "O2-O2.cia"), "O2-O2", 1e-45)
        n2n2 = _write_cia(joinpath(dir, "N2-N2.cia"), "N2-N2", 1e-45)
        o2 = _fixed_broadener_table(1e-27)
        vmr = Dict{String,Any}("O2" => [0.21, 0.21])

        # Variable O₂ with O₂–O₂ CIA: compare against a forward VMR perturbation.
        params = _scene(q=0.0, variable=["O2"], luts=Any[o2], vmr=vmr, cia_files=[o2o2])
        model, lin = model_from_parameters(LinMode(), params)
        Nz, iz, h = 2, 2, 1e-4
        native = lin.τ̇_abs[1][Nz + iz, :, iz]   # O₂ occupies gas slot 2
        τ_at(v) = begin
            p = deepcopy(params)
            p.absorption_params.vmr["O2"] = [0.21, v]
            model_from_parameters(p).τ_abs[1][:, iz]
        end
        fd = (τ_at(0.21 + h) - τ_at(0.21 - h)) / (2h)
        @test _mismatch(native, fd) > 0.1
        @test fd ≈ (τ_at(0.21 + 10h) - τ_at(0.21 - 10h)) / (20h) rtol=1e-6
        @test_throws r"O2 abundance Jacobians.*O2-O2 collision-induced absorption" _full_solve(model, lin)

        # Every CIA pair responds to H₂O through the dry-air number density.
        humid = _scene(fixed=["O2"], luts=Any[o2], vmr=vmr, cia_files=[n2n2])
        model_h, lin_h = model_from_parameters(LinMode(), humid)
        cia = vSmartMOM.Absorption.load_cia_table(n2n2, _GRID)
        function cia_τ(profile)
            τ = zeros(length(_GRID), Nz)
            vSmartMOM.Absorption.compute_τ_cia!(τ, cia, profile, profile.vmr)
            return τ
        end
        @test norm(_fd_in_h2o_ratio(cia_τ, model_h.profile, iz)) > 0
        @test_throws r"q-driven H₂O.*N2-N2 collision-induced absorption" _full_solve(model_h, lin_h)

        # A zero collider VMR removes the opacity but not its derivative.
        o2n2 = _write_cia(joinpath(dir, "O2-N2.cia"), "O2-N2", 1e-45)
        absent = _scene(q=0.0, variable=["O2"], luts=Any[o2],
                        vmr=Dict{String,Any}("O2" => [0.0, 0.0]), cia_files=[o2n2])
        model_z, lin_z = model_from_parameters(LinMode(), absent)
        τ_zero(v) = begin
            p = deepcopy(absent)
            p.absorption_params.vmr["O2"] = [0.0, v]
            model_from_parameters(p).τ_abs[1][:, iz]
        end
        fd_zero = (τ_zero(h) - τ_zero(0.0)) / h   # τ is linear in v_O₂ here
        @test _mismatch(lin_z.τ̇_abs[1][Nz + iz, :, iz], fd_zero) > 0.1
        @test_throws r"O2-N2 collision-induced absorption" _full_solve(model_z, lin_z)

        # The same O₂ column is complete when the configured CIA pair excludes O₂.
        complete = _scene(q=0.0, variable=["O2"], luts=Any[o2], vmr=vmr, cia_files=[n2n2])
        model_c, lin_c = model_from_parameters(LinMode(), complete)
        @test all(isfinite, _full_solve(model_c, lin_c).toa_jacobian)

        # Fixed O₂ makes no abundance request, and the opt-in accepts the partial.
        fixed = _scene(q=0.0, fixed=["O2"], luts=Any[o2], vmr=vmr, cia_files=[o2o2])
        model_x, lin_x = model_from_parameters(LinMode(), fixed)
        @test all(isfinite, _full_solve(model_x, lin_x).toa_jacobian)
        model_p, lin_p = model_from_parameters(LinMode(), params;
            abundance_jacobians=:frozen_cross_section)
        @test lin_p.τ̇_abs == lin.τ̇_abs
        @test all(isfinite, _full_solve(model_p, lin_p).toa_jacobian)
    end
end

# Selection is enforced where a retrieval plan is compiled, before any RT.
struct _AbundancePlanFlavor <: AbstractJacobianFlavor
    columns::Vector{Int}
end
function CoreRT.jacobian_plan(flavor::_AbundancePlanFlavor, params, model, lin)
    keys = [ParameterKey(:custom, Symbol("column_$c")) for c in flavor.columns]
    push!(keys, ParameterKey(:surface, :albedo; component=1, band=1))
    names = String[string(k.field) for k in keys]
    n = length(flavor.columns)
    layout = ActiveParameterLayout(keys, names, keys, names, flavor.columns, n;
                                   surface_columns=(n + 1):(n + 1))
    return JacobianPlan(flavor, keys, names, [layout])
end

@testset "Plans select only supported abundance columns" begin
    params = _scene(h2o=_water_line(), variable=["CH4"],
                    luts=Any[_fixed_broadener_table()],
                    vmr=Dict{String,Any}("CH4" => [1.8e-6, 1.9e-6]))
    # Native: pressure 1; q-H₂O 2:3; variable CH4 4:5.
    model, planned = model_from_parameters(_AbundancePlanFlavor([1, 4, 5]), params)
    result = rt_run(model, planned)
    @test all(isfinite, result.toa_jacobian)
    @test_throws r"q-driven H₂O abundance Jacobians" model_from_parameters(
        _AbundancePlanFlavor([1, 3]), params)
    # A plan does not unlock a full-layout solve of its base model.
    @test_throws r"self broadening" rt_run(
        model, planned.base, 0, size(planned.τ̇_abs[1], 1), 1)
    model_p, planned_p = model_from_parameters(_AbundancePlanFlavor([1, 3]), params;
        abundance_jacobians=:frozen_cross_section)
    @test all(isfinite, rt_run(model_p, planned_p).toa_jacobian)

    # With q ≡ 0 the H₂O columns are structural zeros, so a plan may not select them.
    dry = deepcopy(params)
    dry.q .= 0
    @test_throws r"q ≡ 0" model_from_parameters(_AbundancePlanFlavor([1, 3]), dry)
    @test all(isfinite, rt_run(model_from_parameters(
        _AbundancePlanFlavor([1, 4, 5]), dry)...).toa_jacobian)

    # Hand-assembled tangents keep both historical constructor forms.
    base = planned.base
    fields = (base.τ̇_abs, base.τ̇_aer, base.lin_aerosol_optics,
              base.τ̇_rayl_psurf, base.τ̇_aer_psurf, base.τ̇_abs_psurf)
    @test RTModelLin(fields...).abundance_contract === nothing
    @test typeof(base)(fields...).abundance_contract === nothing

    # Direct contract checks, independent of model construction.
    contract = CoreRT.AbundanceTangentContract(:complete,
        [["q-driven H₂O", "CO2"]], [[["MT_CKD"], String[]]])
    @test CoreRT._require_complete_abundance_tangents(contract, 1, 3:4, 2) === nothing
    @test_throws ArgumentError CoreRT._require_complete_abundance_tangents(
        contract, 1, [2, 3], 2)
    accepted = CoreRT.AbundanceTangentContract(:frozen_cross_section,
        contract.species, contract.omitted)
    @test CoreRT._require_complete_abundance_tangents(accepted, 1, 1:4, 2) === nothing
end
