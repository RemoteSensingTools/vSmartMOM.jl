# Shared fixture for test_vertical_grid_psurf.jl and its 7acab570 golden-file
# generator (reference/generate_vertical_grid_golden.jl). Keep it free of the
# vertical-grid API except inside functions that only the new code calls.
using LinearAlgebra
using vSmartMOM
using vSmartMOM.CoreRT
import AtmosphericAbsorption

# Synthetic spectroscopy (as in test_absorption_pressure.jl): independent of
# downloaded HITRAN/ABSCO data, with a genuine pressure dependence.
function _vg_line_model(::Type{FT}; water=false) where {FT}
    lines = AtmosphericAbsorption.LineDatabase(
        mol=Int32[water ? 1 : 2], iso=Int32[1], ν0=FT[6200],
        S=FT[water ? 1e-25 : 1e-23], E_lower=FT[0], g_upper=FT[1],
        γ_air=FT[0.07], γ_self=FT[0.30], n_air=FT[0.7],
        δ_air=FT[-0.01], n_self=FT[0.5], δ_self=FT[0.04],
        Y_LM=FT[0], n_Y_LM=FT[0.5], γ2_air=FT[0], δ2_air=FT[0],
        νVC=FT[0], η=FT[0], molar_mass=FT[water ? 18 : 44],
        meta=AtmosphericAbsorption.SourceMetadata("synthetic vertical-grid test", 296.0, 1013.25))
    partition = AtmosphericAbsorption.TabulatedPF(FT[100, 300, 500], FT[1, 1, 1])
    AtmosphericAbsorption.LineByLineModel(lines, partition;
        profile=AtmosphericAbsorption.Voigt(),
        cpf=AtmosphericAbsorption.HumlicekWeideman32(), wing_cutoff=FT(10),
        vmr=zero(FT), architecture=AtmosphericAbsorption.CPU())
end

# Pressure nodes at 100/600/1100 hPa; the test layers' mid-pressures (150,
# 375, 625, 875 hPa) stay well inside node intervals for all FD steps.
function _vg_absco(::Type{FT}) where {FT}
    ν = FT[6199.8, 6200, 6200.2]
    p = FT[100, 600, 1100]
    T = FT[210 230 250; 290 310 330]
    broadener = FT[0, 0.05]
    σ = zeros(FT, length(ν), length(broadener), size(T, 1), length(p))
    for ip in eachindex(p), it in axes(T, 1), ib in eachindex(broadener), iν in eachindex(ν)
        σ[iν, ib, it, ip] = FT(1e-23) * FT(iν) *
            (1 + FT(0.001) * p[ip] + FT(0.000001) * p[ip]^2 +
             FT(0.004) * T[it, ip] + 2 * broadener[ib])
    end
    AtmosphericAbsorption.AbscoLUT(2, -1, ν, p, T, broadener, σ;
        architecture=AtmosphericAbsorption.CPU())
end

function _vg_regular_lut(::Type{FT}) where {FT}
    ν, p, T = FT[6199.8, 6200, 6200.2], FT[100, 600, 1100], FT[210, 290, 330]
    σ = FT[FT(1e-23) * FT(i) *
        (1 + FT(0.001) * pp + FT(0.000001) * pp^2 + FT(0.004) * tt)
        for i in eachindex(ν), pp in p, tt in T]
    AtmosphericAbsorption.InterpolationModel(ν, p, T, σ, zero(FT), AtmosphericAbsorption.CPU())
end

_vg_absorber(FT, kind; water=false) =
    kind === :line ? _vg_line_model(FT; water) :
    kind === :absco ? _vg_absco(FT) : _vg_regular_lut(FT)

function _vg_write_cia(path)
    open(path, "w") do io
        for T in (200.0, 300.0)
            ν = [6199.0, 6200.0, 6201.0]
            σ = T == 200.0 ? [2e-46, 3e-46, 4e-46] : [1e-46, 2e-46, 3e-46]
            header = rpad("O2-O2", 20) * lpad(string(first(ν)), 10) *
                     lpad(string(last(ν)), 10) * lpad(string(length(ν)), 7) *
                     lpad(string(T), 7) * lpad("1.000e-40", 10) *
                     repeat(" ", 33) * lpad("11", 3)
            println(io, header)
            for i in eachindex(ν)
                println(io, "$(ν[i]) $(σ[i])")
            end
        end
    end
    return path
end

"""
Four-layer test scene: mixed fixed/variable trace gases on the given absorber,
q-driven H₂O through the same absorber (`water=true`), optional CIA, and two
aerosols — one altitude-form `LogNormal` and one pressure-form `Normal`.
"""
function _vg_parameters(::Type{FT}, kind; water=false, aerosols=false,
                        cia_path=nothing) where {FT}
    config = Dict(
        "radiative_transfer" => Dict(
            "spec_bands" => ["[6199.8, 6199.95, 6200.0, 6200.05, 6200.2]"],
            "surface" => ["LambertianSurfaceScalar(0.2)"],
            "nstreams" => 3, "polarization_type" => "Stokes_I()",
            "truncation" => "auto", "depol" => -1,
            "float_type" => string(FT), "architecture" => "CPU()"),
        "geometry" => Dict("sza" => 35.0, "vza" => [25.0], "vaz" => [40.0],
                           "obs_alt" => [0]),
        "atmospheric_profile" => Dict(
            "T" => [230.0, 250.0, 270.0, 285.0],
            "p" => [50.0, 250.0, 500.0, 750.0, 1000.0],
            "profile_reduction" => -1))
    if aerosols
        config["scattering"] = Dict(
            "aerosols" => [
                Dict("τ_ref" => 0.08, "μ" => 0.15, "σ" => 1.4, "nᵣ" => 1.42,
                     "nᵢ" => 1e-3, "z₀" => 2.5, "σ₀" => 0.5),
                Dict("τ_ref" => 0.05, "μ" => 0.15, "σ" => 1.4, "nᵣ" => 1.42,
                     "nᵢ" => 1e-3, "p₀" => 600.0, "σp" => 150.0)],
            "r_max" => 3.0, "nquad_radius" => 100, "λ_ref" => 0.77,
            "decomp_type" => "NAI2()")
    end
    params = read_parameters(config)
    absorber = _vg_absorber(FT, kind)
    vmr = Dict{String,Any}("CO2" => FT(4e-4), "CH4" => FT[3e-4, 4e-4, 5e-4, 6e-4],
                           "O2" => FT(0.21))
    h2o = water ? Any[_vg_absorber(FT, kind; water=true)] : Any[:disabled]
    cia = cia_path === nothing ? String[] : [cia_path]
    params.absorption_params = CoreRT.AbsorptionParameters(
        [["CO2"]], [["CH4"]], vmr, AtmosphericAbsorption.Voigt(),
        AtmosphericAbsorption.HumlicekWeideman32(), FT(10), [Any[absorber, absorber]],
        h2o, cia, "")
    params.q .= water ? FT[0.002, 0.004, 0.006, 0.01] : zero(FT)
    return params
end

"Copy `params` with the surface pressure moved to `p_surf` under `grid`."
function _vg_at_psurf(params, grid, p_surf)
    moved = deepcopy(params)
    moved.p = surface_pressure_interfaces(grid, params.p, p_surf)
    return moved
end

# These comparisons check the native gas columns (cross sections held fixed),
# so they opt in to that partial explicitly under the abundance contract.
_vg_lin(params, grid; kwargs...) = model_from_parameters(LinMode(), params;
    vertical_grid=grid, compute_aerosol_microphysics_jacobians=false,
    abundance_jacobians=:frozen_cross_section, kwargs...)

"Central FD of an optical field `f(model)` with respect to p_surf."
function _vg_fd(f, params, grid, h; kwargs...)
    p_surf = params.p[end]
    plus = model_from_parameters(_vg_at_psurf(params, grid, p_surf + h); kwargs...)
    minus = model_from_parameters(_vg_at_psurf(params, grid, p_surf - h); kwargs...)
    return (f(plus) - f(minus)) / (2h)
end

_vg_relerr(a, b) = norm(a - b) / norm(b)
