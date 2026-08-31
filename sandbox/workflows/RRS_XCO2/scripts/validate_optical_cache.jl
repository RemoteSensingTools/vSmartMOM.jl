#!/usr/bin/env julia

include(joinpath(@__DIR__, "common.jl"))
using .RRSXCO2Common
using vSmartMOM
using vSmartMOM.CoreRT
using Test

params = load_parameters()
params.architecture = CPU()
# Retain both band extrema and interior samples so the aerosol spectral
# interpolation and phase derivatives remain genuinely wavelength-dependent.
params.spec_bands = [band[unique(round.(Int, range(1, length(band); length=9)))]
                     for band in params.spec_bands]
model, lin_model = model_from_parameters(LinMode(), params; external_solar=true)
rs = vSmartMOM.InelasticScattering.noRS{CoreRT.float_type(model)}()

function samefield(a, b; rtol=2e-13, atol=2e-15)
    a === nothing && return b === nothing
    b === nothing && return false
    return isapprox(Array(a), Array(b); rtol, atol)
end

@testset "public linearized RT cache path" begin
    n_aerosol = length(params.scattering_params.rt_aerosols)
    n_gas = size(lin_model.τ̇_abs[1], 1)
    n_surface = length(params.brdf[1].legendre_coeff)
    result = rt_run(model, lin_model, n_aerosol, n_gas, n_surface; i_band=1)
    @test all(isfinite, Array(result.toa))
    @test all(isfinite, Array(result.toa_jacobian))
end

@testset "m-invariant optical cache" begin
    for ib in 1:3
        bands = (ib,)
        cache = CoreRT.build_m_invariant_cache_lin(bands, model, lin_model)
        for m in (0, 1, 3)
            old, oldlin, _ = CoreRT.constructCoreOpticalProperties(
                rs, bands, m, model, lin_model)
            new, newlin, _ = CoreRT.constructCoreOpticalProperties(
                rs, bands, m, model, lin_model, cache)
            for iz in eachindex(old)
                for name in (:τ, :ϖ, :Z⁺⁺, :Z⁻⁺, :Z₀⁺, :Z₀⁻)
                    @test samefield(getfield(old[iz], name), getfield(new[iz], name))
                end
                for name in (:τ̇, :ϖ̇, :Ż⁺⁺, :Ż⁻⁺, :Ż₀⁺, :Ż₀⁻)
                    @test samefield(getfield(oldlin[iz], name), getfield(newlin[iz], name))
                end
            end
        end
    end

    # The public linearized entry point supplies a scalar i_band. Ensure the
    # cache path normalizes it identically to the internal one-band tuple.
    scalar_cache = CoreRT.build_m_invariant_cache_lin(1, model, lin_model)
    scalar = CoreRT.constructCoreOpticalProperties(rs, 1, 0, model, lin_model,
                                                    scalar_cache)
    tuple_cache = CoreRT.build_m_invariant_cache_lin((1,), model, lin_model)
    tuple = CoreRT.constructCoreOpticalProperties(rs, (1,), 0, model, lin_model,
                                                   tuple_cache)
    @test samefield(scalar[1][1].τ, tuple[1][1].τ)
    @test samefield(scalar[2][1].τ̇, tuple[2][1].τ̇)
end
