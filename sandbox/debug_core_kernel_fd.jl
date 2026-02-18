#!/usr/bin/env julia
# Isolate core RT derivative accuracy from outer parameter propagation.
# This script compares analytic derivatives from rt_kernel_lin core pathways
# against central finite differences for directional perturbations in:
#   1) tau
#   2) omega
#   3) Z (both Z++ and Z-+)
#   4) tau_sum (SFI beam attenuation pathway)
#
# It uses a synthetic "single-layer as TOA" setup so interaction terms are not mixed in.

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

function zero_lin_like(cp::CoreRT.CoreScatteringOpticalProperties, nparams::Int)
    nSpec = length(cp.τ)
    nμ = size(cp.Z⁺⁺, 1)
    τdot = zeros(eltype(cp.τ), nparams, nSpec)
    ωdot = zeros(eltype(cp.τ), nparams, nSpec)
    Zppdot = zeros(eltype(cp.τ), nparams, nμ, nμ, nSpec)
    Zmpdot = zeros(eltype(cp.τ), nparams, nμ, nμ, nSpec)
    return CoreRT.CoreScatteringOpticalPropertiesLin(τdot, ωdot, Zppdot, Zmpdot)
end

function run_single_layer_kernel!(
    model,
    rs_type,
    cp::CoreRT.CoreScatteringOpticalProperties,
    cplin::CoreRT.CoreScatteringOpticalPropertiesLin,
    τ_sum::AbstractVector,
    τ̇_sum::AbstractMatrix;
    m::Int=0,
    SFI::Bool=true,
)
    pol_type = model.params.polarization_type
    quad_points = model.quad_points
    qp_μN = quad_points.qp_μN
    FT = eltype(model.obs_geom.sza)
    arr_type = CoreRT.array_type(model.params.architecture)
    NquadN = quad_points.Nquad * pol_type.n
    nSpec = length(cp.τ)
    dims = (NquadN, NquadN)
    nparams = size(cplin.τ̇, 1)

    added, added_lin = CoreRT.make_added_layer(LinMode(), rs_type, FT, arr_type, nparams, dims, nSpec)
    comp, comp_lin = CoreRT.make_composite_layer(LinMode(), rs_type, FT, arr_type, nparams, dims, nSpec)
    I_static = Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))))

    # For iz==1, scattering_interface is not used in interaction dispatch.
    # We still pass a valid type.
    scattering_interface = CoreRT.ScatteringInterface_01()

    CoreRT.rt_kernel!(
        rs_type,
        pol_type,
        SFI,
        added,
        added_lin,
        comp,
        comp_lin,
        cp,
        cplin,
        scattering_interface,
        τ_sum,
        τ̇_sum,
        m,
        quad_points,
        I_static,
        model.params.architecture,
        qp_μN,
        1,  # synthetic TOA layer
    )

    return comp, comp_lin
end

function run_case(
    model,
    rs_type,
    cp_base::CoreRT.CoreScatteringOpticalProperties,
    τ_sum_base::AbstractVector,
    case_name::String,
    dτ::AbstractVector,
    dω::AbstractVector,
    dZpp,
    dZmp,
    dτsum::AbstractVector,
    eps::Float64,
)
    nSpec = length(cp_base.τ)
    nμ = size(cp_base.Z⁺⁺, 1)
    nparams = 1

    τdot = zeros(Float64, nparams, nSpec)
    ωdot = zeros(Float64, nparams, nSpec)
    Zppdot = zeros(Float64, nparams, nμ, nμ, nSpec)
    Zmpdot = zeros(Float64, nparams, nμ, nμ, nSpec)
    τdot[1, :] .= dτ
    ωdot[1, :] .= dω
    Zppdot[1, :, :, :] .= dZpp
    Zmpdot[1, :, :, :] .= dZmp
    cplin_dir = CoreRT.CoreScatteringOpticalPropertiesLin(τdot, ωdot, Zppdot, Zmpdot)
    τ̇_sum = reshape(dτsum, 1, nSpec)

    # Analytic directional derivative
    comp, comp_lin = run_single_layer_kernel!(model, rs_type, cp_base, cplin_dir, τ_sum_base, τ̇_sum)
    an_T = Array(@view(comp_lin.Ṫ⁺⁺[1, :, :, :]))
    an_R = Array(@view(comp_lin.Ṙ⁻⁺[1, :, :, :]))
    an_Jp = Array(@view(comp_lin.J̇₀⁺[1, :, :, :]))
    an_Jm = Array(@view(comp_lin.J̇₀⁻[1, :, :, :]))

    # Finite-difference directional derivative
    cp_plus = CoreRT.CoreScatteringOpticalProperties(
        cp_base.τ .+ eps .* dτ,
        cp_base.ϖ .+ eps .* dω,
        cp_base.Z⁺⁺ .+ eps .* dZpp,
        cp_base.Z⁻⁺ .+ eps .* dZmp,
    )
    cp_minus = CoreRT.CoreScatteringOpticalProperties(
        cp_base.τ .- eps .* dτ,
        cp_base.ϖ .- eps .* dω,
        cp_base.Z⁺⁺ .- eps .* dZpp,
        cp_base.Z⁻⁺ .- eps .* dZmp,
    )
    zero_lin = zero_lin_like(cp_base, 1)
    τ̇_zero = zeros(Float64, 1, nSpec)

    comp_p, _ = run_single_layer_kernel!(model, rs_type, cp_plus, zero_lin, τ_sum_base .+ eps .* dτsum, τ̇_zero)
    comp_m, _ = run_single_layer_kernel!(model, rs_type, cp_minus, zero_lin, τ_sum_base .- eps .* dτsum, τ̇_zero)

    fd_T = (Array(comp_p.T⁺⁺) .- Array(comp_m.T⁺⁺)) ./ (2eps)
    fd_R = (Array(comp_p.R⁻⁺) .- Array(comp_m.R⁻⁺)) ./ (2eps)
    fd_Jp = (Array(comp_p.J₀⁺) .- Array(comp_m.J₀⁺)) ./ (2eps)
    fd_Jm = (Array(comp_p.J₀⁻) .- Array(comp_m.J₀⁻)) ./ (2eps)

    tmax, tmean, tn = rel_stats(an_T, fd_T)
    rmax, rmean, rn = rel_stats(an_R, fd_R)
    jpmax, jpmean, jpn = rel_stats(an_Jp, fd_Jp)
    jmmax, jmmean, jmn = rel_stats(an_Jm, fd_Jm)

    println("case=$case_name (eps=$eps)")
    @printf("  T++: max=%.4e mean=%.4e n=%d\n", tmax, tmean, tn)
    @printf("  R-+: max=%.4e mean=%.4e n=%d\n", rmax, rmean, rn)
    @printf("  J0+: max=%.4e mean=%.4e n=%d\n", jpmax, jpmean, jpn)
    @printf("  J0-: max=%.4e mean=%.4e n=%d\n", jmmax, jmmean, jmn)
    println()
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

    # Choose a strongly scattering layer, but run it in synthetic TOA mode.
    iz = 4
    cp = layer_opt[iz]
    τ_sum_base = τ_sum_all[:, iz]

    nSpec = length(cp.τ)
    nμ = size(cp.Z⁺⁺, 1)

    dτ = ones(Float64, nSpec)
    dω = ones(Float64, nSpec)
    dZpp = Array(cp.Z⁺⁺)
    dZmp = Array(cp.Z⁻⁺)
    dτsum = ones(Float64, nSpec)
    zeros_τ = zeros(Float64, nSpec)
    zeros_Z = zeros(Float64, nμ, nμ, nSpec)

    println("core-only rt_kernel directional derivative check")
    println("layer=$iz, nSpec=$nSpec, nμ=$nμ, m=$m")
    println("(single-layer TOA harness; outer propagation is bypassed)\n")

    run_case(model, rs_type, cp, τ_sum_base, "tau", dτ, zeros_τ, zeros_Z, zeros_Z, zeros_τ, 1e-6)
    run_case(model, rs_type, cp, τ_sum_base, "omega", zeros_τ, dω, zeros_Z, zeros_Z, zeros_τ, 1e-6)
    run_case(model, rs_type, cp, τ_sum_base, "Z", zeros_τ, zeros_τ, dZpp, dZmp, zeros_τ, 1e-6)
    run_case(model, rs_type, cp, τ_sum_base, "tau_sum", zeros_τ, zeros_τ, zeros_Z, zeros_Z, dτsum, 1e-6)
end

main()
