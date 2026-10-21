# Isolated adapter comparison: only the parameter-copy policy changes.
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

function copy_probe(params,share)
    t = @timed copy_parameters(params;share_luts=share)
    original = params.absorption_params.luts[1][1].itp.coefs
    copied = t.value.absorption_params.luts[1][1].itp.coefs
    @assert (original === copied) == share
    (;seconds=t.time,bytes=t.bytes)
end

const saved = ENV["REPLAY_STATE"]
const output_dir = ENV["REPLAY_OUTPUT"]
mkpath(output_dir)
state, noise = NCDataset(saved) do ds
    (Array(ds["final_state"][:]),Array(ds["noise_standard_deviation"][:]))
end
println("PACKAGE ",pathof(vSmartMOM));flush(stdout)
evaluator = quiet(()->OCOForwardEvaluator(;architecture=:GPU,float_type=Float32,nstreams=9))
template_before = template_state(evaluator.base_parameters)
tables_before = table_hashes(evaluator.base_parameters)
println("PREPARED grids=",length.(evaluator.base_parameters.spec_bands));flush(stdout)

source = read(joinpath(study,"inversion/VSmartMOMForward.jl"),String)
begin_at = findfirst("function evaluate_oco_forward(",source).start
stop_at = findfirst("\n(evaluator::OCOForwardEvaluator)",source).start
original_method = source[begin_at:prevind(source,stop_at)]
solver_call = "rt_run_lin(model, lin_model; i_band=iband, sources)"
copy_call = "deepcopy(evaluator.base_parameters)"
@assert occursin(solver_call,original_method) && occursin(copy_call,original_method)
for mode in (:copied,:shared)
    method = replace(original_method,"function evaluate_oco_forward(" =>
        "function evaluate_oco_forward_$mode(", solver_call =>
        "rt_run_lin(model, lin_model; i_band=iband, sources, jacobian_basis=:local, jacobian_adding=:source)")
    mode == :shared && (method=replace(method,copy_call =>
        "copy_parameters(evaluator.base_parameters; share_luts=true)"))
    Base.include_string(VSmartMOMForward,method,"isolated_$(mode)_lut_replay.jl")
end

copy_records = Dict()
for share in (false,true)
    copy_probe(evaluator.base_parameters,share) # warm the copy helper
    samples = map(1:3) do _
        GC.gc()
        copy_probe(evaluator.base_parameters,share)
    end
    copy_records[string(share)] = Dict("seconds"=>median(s.seconds for s in samples),
        "bytes"=>median(s.bytes for s in samples))
end
println("COPY_PROBES ",copy_records);flush(stdout)

records = []
outputs = Dict()
for mode in (:copied,:shared)
    fn = getproperty(VSmartMOMForward,Symbol("evaluate_oco_forward_$mode"))
    println("WARM mode=$mode");flush(stdout)
    quiet(()->fn(evaluator,state));CUDA.synchronize()
    samples = []
    for trial in 1:3
        GC.gc();CUDA.reclaim();CUDA.synchronize()
        timed = @timed begin
            result = quiet(()->fn(evaluator,state))
            CUDA.synchronize()
            result
        end
        result = timed.value
        record = Dict("trial"=>trial,"seconds"=>timed.time,"host_bytes"=>timed.bytes,
            "timing"=>Dict(string(k)=>v for (k,v) in pairs(result.timing)))
        push!(samples,record)
        println("SAMPLE mode=$mode ",record);flush(stdout)
        outputs[mode] = result
    end
    result = outputs[mode]
    JLD2.jldsave(joinpath(output_dir,"$mode-output.jld2");
        y=result.measurement,K=result.jacobian,state)
    push!(records,Dict("mode"=>string(mode),"samples"=>samples))
    open(io->TOML.print(io,Dict("package"=>pathof(vSmartMOM),"saved_state"=>saved,
        "records"=>records,"copy_probes"=>copy_records)),joinpath(output_dir,"timings.toml"),"w")
end
@assert isequal(outputs[:copied].measurement,outputs[:shared].measurement)
@assert isequal(outputs[:copied].jacobian,outputs[:shared].jacobian)

# A different trial exercises all state mutation families before replaying A.
perturbed = copy(state)
perturbed[1] += 0.1
perturbed[VSmartMOMForward.CO2_RANGE] .*= 1.001
perturbed[VSmartMOMForward.LOG_AOD_RANGE] .+= 0.001
perturbed[VSmartMOMForward.LOG_HEIGHT_RANGE] .+= 0.001
perturbed[first(VSmartMOMForward.SURFACE_RANGE)] += 0.001
perturbed[first(VSmartMOMForward.SIF_RANGE)] += 1e-6
println("ISOLATION A -> B -> A");flush(stdout)
copied_b = quiet(()->VSmartMOMForward.evaluate_oco_forward_copied(evaluator,perturbed))
shared_b = quiet(()->VSmartMOMForward.evaluate_oco_forward_shared(evaluator,perturbed))
shared_a = quiet(()->VSmartMOMForward.evaluate_oco_forward_shared(evaluator,state))
@assert isequal(copied_b.measurement,shared_b.measurement)
@assert isequal(copied_b.jacobian,shared_b.jacobian)
@assert isequal(outputs[:shared].measurement,shared_a.measurement)
@assert isequal(outputs[:shared].jacobian,shared_a.jacobian)
@assert !isequal(shared_b.measurement,shared_a.measurement)
@assert isequal(template_before,template_state(evaluator.base_parameters))
@assert tables_before == table_hashes(evaluator.base_parameters)
for (name,result,trial_state) in (("copied-b",copied_b,perturbed),
                                 ("shared-b",shared_b,perturbed),("shared-a",shared_a,state))
    JLD2.jldsave(joinpath(output_dir,"$name-output.jld2");
        y=result.measurement,K=result.jacobian,state=trial_state)
end
open(io->TOML.print(io,Dict("exact_copy_policy_parity"=>true,
    "exact_perturbed_state_parity"=>true,"exact_repeat_state_parity"=>true,
    "template_unchanged"=>true,"table_coefficients_unchanged"=>true,
    "table_sha256"=>tables_before)),joinpath(output_dir,"isolation.toml"),"w")
println("PASS: exact outputs for both states; template and LUT coefficients unchanged.")
