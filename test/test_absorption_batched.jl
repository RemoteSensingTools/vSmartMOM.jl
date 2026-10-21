# Batched vs per-layer absorption accumulation: bit-exact, on real tables.
#
# `compute_absorption_profile!` gained a single-kernel fast path for tabulated
# models (AtmosphericAbsorption.compute_cross_section_profile). Its claim is
# bit-exactness against the per-layer loop on the CPU backend — same bracket
# and clamp semantics, same blend order, same per-element vcd*vmr weighting.
# This test makes both paths fill tau_abs on a REAL ABSCO table (skipping
# loudly when the table is absent — campaign-machine test) and on a synthetic
# InterpolationModel, and compares with `==`.
#
#   julia -t 4 --project=. test/test_absorption_batched.jl
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using Test
using vSmartMOM
using vSmartMOM.CoreRT: compute_absorption_profile!, AtmosphericProfile
import AtmosphericAbsorption
using AtmosphericAbsorption: Architectures

if !isdefined(AtmosphericAbsorption, :compute_cross_section_profile)
    @warn "installed AtmosphericAbsorption lacks compute_cross_section_profile — SKIPPING"
    exit(0)
end

# A small synthetic profile in the shape compute_absorption_profile! consumes.
function fake_profile(FT, n)
    p_full = collect(FT, range(30.0, 1000.0; length = n))
    T = collect(FT, range(210.0, 295.0; length = n))
    vmr_h2o = collect(FT, range(0.0, 0.02; length = n))
    vcd_dry = collect(FT, range(2.0e21, 4.0e21; length = n))
    q = vmr_h2o ./ 2
    p_half = collect(FT, range(25.0, 1010.0; length = n + 1))
    vcd_h2o = vcd_dry .* vmr_h2o
    Δz = fill(FT(500), n)
    return AtmosphericProfile(T, p_full, q, p_half, vmr_h2o, vcd_dry, vcd_h2o,
                              Dict{String, FT}("CO2" => FT(4.0e-4)), Δz)
end

const ABSCO = expanduser("~/data/ABSCO/v5.2_final/co2_v52.hdf")

@testset "batched absorption accumulation" begin
    FT = Float32
    n = 24
    prof = fake_profile(FT, n)
    grid = collect(FT, range(6180.0, 6260.0; length = 2001))
    vmr_gas = fill(FT(4.0e-4), n)

    models = Any[]
    if isfile(ABSCO)
        push!(models, ("ABSCO co2_v52",
                       AtmosphericAbsorption.read_oco2_absco(ABSCO, :wco2; FT,
                           architecture = Architectures.CPU(), broadener_vmr = :all)))
    else
        @warn "no ABSCO table at $ABSCO — nothing to test on this machine"
        exit(0)
    end
    # NOTE deliberately no AA InterpolationModel arm here: vSmartMOM's per-layer
    # loop never had an `absorption_cross_section` method for that type (JLD2
    # tables load as vSmartMOM's own legacy model), so the only tabulated model
    # both paths serve is AbscoLUT. AA-side bitwise coverage for
    # InterpolationModel lives in AtmosphericAbsorption/test_profile_batch.jl.

    for (name, model) in models
        @testset "$name" begin
            τ_loop = zeros(FT, length(grid), n)
            τ_fast = zeros(FT, length(grid), n)
            compute_absorption_profile!(τ_loop, model, grid, vmr_gas, prof;
                                        batched = false)
            compute_absorption_profile!(τ_fast, model, grid, vmr_gas, prof;
                                        batched = true)
            @test τ_fast == τ_loop                      # bit-exact, not approx
            @test any(!iszero, τ_fast)                  # and not trivially zero
            # scalar gas VMR and an explicit broadener vector go the same way
            τ_loop2 = zeros(FT, length(grid), n); τ_fast2 = zeros(FT, length(grid), n)
            compute_absorption_profile!(τ_loop2, model, grid, FT(4.0e-4), prof;
                                        self_broadener_vmr = prof.vmr_h2o,
                                        batched = false)
            compute_absorption_profile!(τ_fast2, model, grid, FT(4.0e-4), prof;
                                        self_broadener_vmr = prof.vmr_h2o,
                                        batched = true)
            @test τ_fast2 == τ_loop2
        end
    end
end
println("batched absorption tests passed")
