# From test/: julia --project=. ../docs/dev_notes/pressure_source_followup_2026-09-08/batch_update_failure_probe.jl
# Historical defect reproduction, retained as audit evidence. Its assertions
# intentionally expect the old broken behavior and fail on the corrected code.
# For current regression coverage use test/test_update_model_transaction.jl.
# Uses synthetic data only.
using vSmartMOM
using vSmartMOM.CoreRT
using Interpolations
import AtmosphericAbsorption
using LinearAlgebra

for reduction in (-1, 1)
println("PROFILE_REDUCTION = ", reduction)
params = read_parameters(Dict(
    "radiative_transfer" => Dict(
        "spec_bands" => ["[6199.8, 6200.0, 6200.2]"],
        "surface" => ["LambertianSurfaceScalar(0.2)"],
        "nstreams" => 3, "polarization_type" => "Stokes_I()",
        "truncation" => "NoTruncation()", "depol" => -1,
        "float_type" => "Float64", "architecture" => "CPU()"),
    "geometry" => Dict("sza" => 35.0, "vza" => [25.0], "vaz" => [40.0], "obs_alt" => [0]),
    "atmospheric_profile" => Dict("T" => [270.0, 280.0],
        "p" => [100.0, 500.0, 1000.0], "profile_reduction" => reduction)))
ν = range(6199.0, 6201.0; length=3)
p = range(100.0, 1100.0; length=3)
T = range(200.0, 300.0; length=3)
table = [1e-24 * (1 + 0.001pp + 0.002tt + i) for i in eachindex(ν), pp in p, tt in T]
itp = interpolate(table, BSpline(Linear()))
absorber = vSmartMOM.Absorption.InterpolationModel(itp, 2, 1, ν, p, T)
params.absorption_params = CoreRT.AbsorptionParameters(
    [["CO2"]], [String[]], Dict{String,Any}("CO2"=>4e-4),
    AtmosphericAbsorption.Voigt(), AtmosphericAbsorption.HumlicekWeideman32(),
    10.0, [Any[absorber]], Any[:disabled], String[], "")
ctx = BatchContext(params)
println("current_T_aliases_live_profile = ", ctx.current_T === ctx.model.profile.T)
println("params_T_aliases_live_profile = ", params.T === ctx.model.profile.T)
old_model_T = copy(ctx.model.profile.T)
old_current_T = copy(ctx.current_T)
old_τ = copy(ctx.model.τ_abs[1])
old_radiance = rt_run(ctx.model).toa
failure = try
    update_model!(ctx; T=[350.0, 360.0])
    nothing
catch err
    err
end
println("exception_type = ", typeof(failure))
println("exception = ", sprint(showerror, failure))
println("old_model_T = ", old_model_T)
println("model_T_after_failed_update = ", ctx.model.profile.T)
println("remembered_T_after_failed_update = ", ctx.current_T)
println("remembered_T_unchanged = ", ctx.current_T == old_current_T)
println("live_model_T_changed = ", ctx.model.profile.T != old_model_T)
println("old_absorption_nonzero = ", any(!iszero, old_τ))
println("absorption_now_all_zero = ", all(iszero, ctx.model.τ_abs[1]))
new_radiance = rt_run(ctx.model).toa
println("rt_run_accepts_failed_context_model = true")
println("relative_radiance_change = ", norm(new_radiance-old_radiance)/norm(old_radiance))
@assert failure !== nothing
@assert (ctx.current_T == old_current_T) == (reduction != -1)
@assert ctx.model.profile.T != old_model_T
@assert any(!iszero, old_τ)
@assert all(iszero, ctx.model.τ_abs[1])
println("REPRODUCTION_CONFIRMED")

end
