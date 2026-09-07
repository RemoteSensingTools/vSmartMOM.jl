# Isolated full-OE comparison: direct physical/matrix versus local/source + shared LUTs.
# Run from test/ with the study data environment in ../suniti_replay/README.md.
using vSmartMOM, CUDA, NCDatasets, JLD2, Statistics, LinearAlgebra, Logging, TOML, SHA
push!(LOAD_PATH,pkgdir(vSmartMOM))
CUDA.device!(0)
CUDA.allowscalar(false)
const study = ENV["STUDY_ROOT"]
include(joinpath(study,"inversion/OptimalEstimation.jl"))
include(joinpath(study,"inversion/VSmartMOMForward.jl"))
using .VSmartMOMForward

quiet(f) = redirect_stdout(devnull) do
    with_logger(NullLogger()) do
        f()
    end
end

function template_state(params)
    (; p=copy(params.p),T=copy(params.T),q=copy(params.q),
       bands=deepcopy(params.spec_bands),vmr=deepcopy(params.absorption_params.vmr),
       aerosol=[(a.τ_ref,a.profile.μ,a.profile.σ) for a in params.scattering_params.rt_aerosols],
       surface=[copy(s.legendre_coeff) for s in params.brdf])
end

function table_hashes(params)
    tables = Any[lut for band in params.absorption_params.luts for lut in band]
    append!(tables,[lut for lut in params.absorption_params.h2o_lut
                    if lut !== nothing && lut !== :disabled])
    map(tables) do lut
        coefficients = lut.itp.coefs
        while parent(coefficients) !== coefficients
            coefficients = parent(coefficients)
        end
        bytes2hex(sha256(reinterpret(UInt8,vec(coefficients))))
    end
end


using .OptimalEstimation
const output_dir = ENV["REPLAY_OUTPUT"]
mkpath(output_dir)
source = read(joinpath(study,"inversion/VSmartMOMForward.jl"),String)
begin_at = findfirst("function evaluate_oco_forward(",source).start
stop_at = findfirst("\n(evaluator::OCOForwardEvaluator)",source).start
original = source[begin_at:prevind(source,stop_at)]
@assert occursin("enumerate(BAND_SPECS)",original)
original = replace(original,"enumerate(BAND_SPECS)"=>"enumerate(BAND_SPECS[1:1])",
    "undef, 3)"=>"undef, 1)")
solver = "rt_run_lin(model, lin_model; i_band=iband, sources)"
@assert occursin(solver,original)
for mode in (:reference,:optimized)
    method = replace(original,"function evaluate_oco_forward("=>"function evaluate_$mode(",
        solver => mode == :reference ?
        "rt_run_lin(model, lin_model; i_band=iband, sources, jacobian_basis=:physical, jacobian_adding=:matrix)" :
        "rt_run_lin(model, lin_model; i_band=iband, sources, jacobian_basis=:local, jacobian_adding=:source)")
    mode == :optimized && (method=replace(method,"deepcopy(evaluator.base_parameters)"=>
        "copy_parameters(evaluator.base_parameters; share_luts=true)"))
    Base.include_string(VSmartMOMForward,method,"isolated_$mode.jl")
end
evaluator = quiet(()->OCOForwardEvaluator(;architecture=:GPU,float_type=Float64,nstreams=9))
VSmartMOMForward.set_fixed_upper_co2_ppm!(evaluator,400.0)

for mode in (:reference,:optimized)
    state=JLD2.load(joinpath(output_dir,"state035_corrected_siffalse-$mode.jld2"),"state")
    println("FLOAT64 O2 state=$mode");flush(stdout)
    result=quiet(()->VSmartMOMForward.evaluate_optimized(evaluator,state))
    JLD2.jldsave(joinpath(output_dir,"precision64-$mode.jld2");state,y=result.measurement,K=result.jacobian)
end
println("FLOAT64 COMPLETE")
