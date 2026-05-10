#!/usr/bin/env julia

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering: noRS
using LinearAlgebra, Statistics, Printf

const YAML_FAST = "test/test_parameters/JacobianTestFast.yaml"

function rel_stats(analytic, fd; threshold=1e-12)
    mask = abs.(fd) .> threshold
    if !any(mask)
        return (NaN, NaN, 0)
    end
    errs = abs.(analytic[mask] .- fd[mask]) ./ abs.(fd[mask])
    return (maximum(errs), mean(errs), count(mask))
end

function print_stats(label, analytic, fd)
    mx, mn, n = rel_stats(analytic, fd)
    @printf("  %-20s max=%.4e mean=%.4e n=%d\n", label, mx, mn, n)
end

function setup_rs_type(model, iBand::Int)
    pol_type = model.params.polarization_type
    nSpec = size(model.τ_abs[iBand], 1)
    rs = noRS()
    rs.bandSpecLim = [1:nSpec]
    F0 = zeros(Float64, pol_type.n, nSpec)
    F0[1, :] .= 1.0
    rs.F₀ = F0
    return rs
end

function alloc_layers(model, rs_type, nparams::Int, nSpec::Int)
    pol_type = model.params.polarization_type
    quad_points = model.quad_points
    FT = eltype(model.obs_geom.sza)
    arr_type = CoreRT.array_type(model.params.architecture)
    NquadN = quad_points.Nquad * pol_type.n
    dims = (NquadN, NquadN)
    added, added_lin = CoreRT.make_added_layer(LinMode(), rs_type, FT, arr_type, nparams, dims, nSpec)
    I_static = Diagonal(arr_type(Diagonal{FT}(ones(NquadN))))
    return added, added_lin, I_static
end

function make_tau_only_lin(cp::CoreRT.CoreScatteringOpticalProperties, dτ_dir::AbstractVector)
    nSpec = length(cp.τ)
    nμ = size(cp.Z⁺⁺, 1)
    τdot = zeros(Float64, 1, nSpec)
    ωdot = zeros(Float64, 1, nSpec)
    Zppdot = zeros(Float64, 1, nμ, nμ, nSpec)
    Zmpdot = zeros(Float64, 1, nμ, nμ, nSpec)
    τdot[1, :] .= dτ_dir
    return CoreRT.CoreScatteringOpticalPropertiesLin(τdot, ωdot, Zppdot, Zmpdot)
end

function run_elemental_once!(
    model,
    rs_type,
    cp::CoreRT.CoreScatteringOpticalProperties,
    dτ::AbstractVector,
    τ_sum::AbstractVector,
    τ̇_sum::AbstractMatrix,
    ndoubl::Int;
    m::Int=0,
    SFI::Bool=true,
)
    nSpec = length(cp.τ)
    added, added_lin, _ = alloc_layers(model, rs_type, 1, nSpec)
    CoreRT.elemental!(
        model.params.polarization_type,
        SFI,
        τ_sum,
        τ̇_sum,
        dτ,
        rs_type.F₀,
        cp,
        m,
        ndoubl,
        true,
        model.quad_points,
        added,
        added_lin,
        model.params.architecture,
    )
    return added, added_lin
end

function run_tau_analytic_through_doubling!(
    model,
    rs_type,
    cp::CoreRT.CoreScatteringOpticalProperties,
    cplin_tau::CoreRT.CoreScatteringOpticalPropertiesLin,
    dτ::AbstractVector,
    τ_sum::AbstractVector,
    τ̇_sum::AbstractMatrix,
    ndoubl::Int;
    m::Int=0,
    SFI::Bool=true,
)
    nSpec = length(cp.τ)
    added, added_lin, I_static = alloc_layers(model, rs_type, 1, nSpec)

    CoreRT.elemental!(
        model.params.polarization_type,
        SFI,
        τ_sum,
        τ̇_sum,
        dτ,
        rs_type.F₀,
        cp,
        m,
        ndoubl,
        true,
        model.quad_points,
        added,
        added_lin,
        model.params.architecture,
    )

    CoreRT.lin_added_layer_all_params!(
        rs_type,
        model.params.polarization_type,
        SFI,
        model.quad_points,
        cplin_tau,
        added_lin,
        model.params.architecture,
        ndoubl,
    )

    dτ̇_allparams = cplin_tau.τ̇ ./ (2.0^ndoubl)
    expk = CoreRT.array_type(model.params.architecture)(exp.(-dτ / model.quad_points.μ₀))

    CoreRT.doubling_allparams!(
        model.params.polarization_type,
        SFI,
        expk,
        ndoubl,
        added,
        added_lin,
        I_static,
        model.params.architecture,
        dτ̇_allparams,
        model.quad_points.μ₀,
    )

    return added, added_lin
end

function run_forward_through_doubling!(
    model,
    rs_type,
    cp::CoreRT.CoreScatteringOpticalProperties,
    dτ::AbstractVector,
    τ_sum::AbstractVector,
    ndoubl::Int;
    m::Int=0,
    SFI::Bool=true,
)
    nSpec = length(cp.τ)
    added, added_lin, I_static = alloc_layers(model, rs_type, 1, nSpec)
    τ̇_sum_zero = zeros(Float64, 1, nSpec)
    CoreRT.elemental!(
        model.params.polarization_type,
        SFI,
        τ_sum,
        τ̇_sum_zero,
        dτ,
        rs_type.F₀,
        cp,
        m,
        ndoubl,
        true,
        model.quad_points,
        added,
        added_lin,
        model.params.architecture,
    )
    dτ̇_zero = zeros(Float64, 1, nSpec)
    expk = CoreRT.array_type(model.params.architecture)(exp.(-dτ / model.quad_points.μ₀))
    CoreRT.doubling_allparams!(
        model.params.polarization_type,
        SFI,
        expk,
        ndoubl,
        added,
        added_lin,
        I_static,
        model.params.architecture,
        dτ̇_zero,
        model.quad_points.μ₀,
    )
    return added
end

function main()
    params = parameters_from_yaml(YAML_FAST)
    params.architecture = CPU()
    model, lin_model = model_from_parameters(LinMode(), params)
    rs_type = setup_rs_type(model, 1)

    m = 0
    iBand = 1
    layer_opt, layer_opt_lin, _ = CoreRT.constructCoreOpticalProperties(rs_type, iBand, m, model, lin_model)
    _, τ_sum_all, _ = CoreRT.extractEffectiveProps(layer_opt, layer_opt_lin)

    iz = 4
    cp = layer_opt[iz]
    τ_sum = Array(τ_sum_all[:, iz])
    nSpec = length(cp.τ)

    dτ_max = minimum([maximum(cp.τ .* cp.ϖ), 0.001 * minimum(model.quad_points.qp_μ)])
    _, ndoubl = CoreRT.doubling_number(dτ_max, maximum(cp.τ .* cp.ϖ))
    dτ = Array(cp.τ ./ (2.0^ndoubl))

    τ_dir = ones(Float64, nSpec)
    τ̇_sum_zero = zeros(Float64, 1, nSpec)
    cplin_tau = make_tau_only_lin(cp, τ_dir)
    dτ_dir = τ_dir ./ (2.0^ndoubl)

    println("tau stagewise diagnostic (fixed ndoubl)")
    @printf("layer=%d ndoubl=%d nSpec=%d\n\n", iz, ndoubl, nSpec)

    # Stage 1: elemental only, derivative wrt elemental dτ
    added0, added_lin0 = run_elemental_once!(model, rs_type, cp, dτ, τ_sum, τ̇_sum_zero, ndoubl)
    eps_dτ = 1e-7
    added_p, _ = run_elemental_once!(model, rs_type, cp, dτ .+ eps_dτ .* ones(nSpec), τ_sum, τ̇_sum_zero, ndoubl)
    added_m, _ = run_elemental_once!(model, rs_type, cp, dτ .- eps_dτ .* ones(nSpec), τ_sum, τ̇_sum_zero, ndoubl)

    fd_elem_T_dτ = (Array(added_p.t⁺⁺) .- Array(added_m.t⁺⁺)) ./ (2eps_dτ)
    fd_elem_R_dτ = (Array(added_p.r⁻⁺) .- Array(added_m.r⁻⁺)) ./ (2eps_dτ)
    fd_elem_Jp_dτ = (Array(added_p.j₀⁺) .- Array(added_m.j₀⁺)) ./ (2eps_dτ)
    fd_elem_Jm_dτ = (Array(added_p.j₀⁻) .- Array(added_m.j₀⁻)) ./ (2eps_dτ)

    println("elemental core derivative vs FD (w.r.t. dτ)")
    print_stats("T++", Array(@view(added_lin0.ṫ⁺⁺[1,:,:,:])), fd_elem_T_dτ)
    print_stats("R-+", Array(@view(added_lin0.ṙ⁻⁺[1,:,:,:])), fd_elem_R_dτ)
    print_stats("J0+", Array(@view(added_lin0.J̇₀⁺[1,:,:,:])), fd_elem_Jp_dτ)
    print_stats("J0-", Array(@view(added_lin0.J̇₀⁻[1,:,:,:])), fd_elem_Jm_dτ)
    println()

    # Stage 2: elemental mapped to full-τ derivative (before doubling)
    eps_τ = 1e-6
    δdτ = eps_τ .* dτ_dir
    added_pτ, _ = run_elemental_once!(model, rs_type, cp, dτ .+ δdτ, τ_sum, τ̇_sum_zero, ndoubl)
    added_mτ, _ = run_elemental_once!(model, rs_type, cp, dτ .- δdτ, τ_sum, τ̇_sum_zero, ndoubl)
    fd_elem_T_τ = (Array(added_pτ.t⁺⁺) .- Array(added_mτ.t⁺⁺)) ./ (2eps_τ)
    fd_elem_R_τ = (Array(added_pτ.r⁻⁺) .- Array(added_mτ.r⁻⁺)) ./ (2eps_τ)
    fd_elem_Jp_τ = (Array(added_pτ.j₀⁺) .- Array(added_mτ.j₀⁺)) ./ (2eps_τ)
    fd_elem_Jm_τ = (Array(added_pτ.j₀⁻) .- Array(added_mτ.j₀⁻)) ./ (2eps_τ)

    an_elem_T_τ = Array(@view(added_lin0.ṫ⁺⁺[1,:,:,:])) .* reshape(dτ_dir, 1, 1, nSpec)
    an_elem_R_τ = Array(@view(added_lin0.ṙ⁻⁺[1,:,:,:])) .* reshape(dτ_dir, 1, 1, nSpec)
    an_elem_Jp_τ = Array(@view(added_lin0.J̇₀⁺[1,:,:,:])) .* reshape(dτ_dir, 1, 1, nSpec)
    an_elem_Jm_τ = Array(@view(added_lin0.J̇₀⁻[1,:,:,:])) .* reshape(dτ_dir, 1, 1, nSpec)

    println("elemental mapped to full-τ vs FD (before chain/doubling)")
    print_stats("T++", an_elem_T_τ, fd_elem_T_τ)
    print_stats("R-+", an_elem_R_τ, fd_elem_R_τ)
    print_stats("J0+", an_elem_Jp_τ, fd_elem_Jp_τ)
    print_stats("J0-", an_elem_Jm_τ, fd_elem_Jm_τ)
    println()

    # Stage 3: analytic through chain + doubling, compare to FD through forward + doubling
    added_an, added_lin_an = run_tau_analytic_through_doubling!(
        model, rs_type, cp, cplin_tau, dτ, τ_sum, τ̇_sum_zero, ndoubl
    )
    added_fd_p = run_forward_through_doubling!(model, rs_type, cp, dτ .+ δdτ, τ_sum, ndoubl)
    added_fd_m = run_forward_through_doubling!(model, rs_type, cp, dτ .- δdτ, τ_sum, ndoubl)

    fd_dbl_T_τ = (Array(added_fd_p.t⁺⁺) .- Array(added_fd_m.t⁺⁺)) ./ (2eps_τ)
    fd_dbl_R_τ = (Array(added_fd_p.r⁻⁺) .- Array(added_fd_m.r⁻⁺)) ./ (2eps_τ)
    fd_dbl_Jp_τ = (Array(added_fd_p.j₀⁺) .- Array(added_fd_m.j₀⁺)) ./ (2eps_τ)
    fd_dbl_Jm_τ = (Array(added_fd_p.j₀⁻) .- Array(added_fd_m.j₀⁻)) ./ (2eps_τ)

    println("after chain + doubling (ap_ derivatives) vs FD")
    print_stats("T++", Array(@view(added_lin_an.ap_ṫ⁺⁺[1,:,:,:])), fd_dbl_T_τ)
    print_stats("R-+", Array(@view(added_lin_an.ap_ṙ⁻⁺[1,:,:,:])), fd_dbl_R_τ)
    print_stats("J0+", Array(@view(added_lin_an.ap_J̇₀⁺[1,:,:,:])), fd_dbl_Jp_τ)
    print_stats("J0-", Array(@view(added_lin_an.ap_J̇₀⁻[1,:,:,:])), fd_dbl_Jm_τ)

    return nothing
end

main()
