# Isolated single-state replay of the inspected RRS_XCO2 adapter.
# See README.md for required data/environment and the exact source hashes.
# Only outputs under REPLAY_OUTPUT are written; no campaign runner is included.
using vSmartMOM, CUDA, NCDatasets, JLD2, Statistics, LinearAlgebra, Logging, TOML
push!(LOAD_PATH,pkgdir(vSmartMOM)) # expose the study's direct workflow dependencies
CUDA.device!(0)
CUDA.allowscalar(false)
const study = ENV["STUDY_ROOT"]
include(joinpath(study,"inversion/OptimalEstimation.jl"))
include(joinpath(study,"inversion/VSmartMOMForward.jl"))
using .VSmartMOMForward
function quiet(f)
    redirect_stdout(devnull) do
        with_logger(NullLogger()) do
            f()
        end
    end
end

const saved = ENV["REPLAY_STATE"]
mkpath(ENV["REPLAY_OUTPUT"])
state, expected_y, expected_K, noise = NCDataset(saved) do ds
    (Array(ds["final_state"][:]),Array(ds["final_forward_model"][:]),
     Array(ds["final_jacobian"][:,:]),Array(ds["noise_standard_deviation"][:]))
end
println("PACKAGE ",pathof(vSmartMOM));flush(stdout)
evaluator = quiet(()->OCOForwardEvaluator(;architecture=:GPU,float_type=Float32,nstreams=9))
println("PREPARED grids=",length.(evaluator.base_parameters.spec_bands));flush(stdout)
GC.gc()
copy_test = @timed deepcopy(evaluator.base_parameters)
original = evaluator.base_parameters.absorption_params.luts[1][1]
copied = copy_test.value.absorption_params.luts[1][1]
println("PARAMETER_COPY seconds=$(copy_test.time) bytes=$(copy_test.bytes) lut_type=$(typeof(original)) same_lut=$(original === copied) same_interpolator=$(original.itp === copied.itp)")
copy_test=nothing; copied=nothing; GC.gc(); CUDA.reclaim()

# The study adapter and state mapping are unchanged. Only the solver keywords
# vary for the optimized comparison, and all outputs remain under /tmp.
modes = get(ENV,"REPLAY_FAST","false")=="true" ? (:matrix,:source) : (:matrix,)
if :source in modes
    source = read(joinpath(study,"inversion/VSmartMOMForward.jl"),String)
    begin_at=findfirst("function evaluate_oco_forward(",source).start
    stop_at=findfirst("\n(evaluator::OCOForwardEvaluator)",source).start
    method=source[begin_at:prevind(source,stop_at)]
    method=replace(method,"function evaluate_oco_forward("=>"function evaluate_oco_forward_fast(")
    @assert occursin("rt_run_lin(model, lin_model; i_band=iband, sources)",method) "Study adapter solver call changed; review the replay substitution."
    method=replace(method,"rt_run_lin(model, lin_model; i_band=iband, sources)"=>
        "rt_run_lin(model, lin_model; i_band=iband, sources, jacobian_basis=:local, jacobian_adding=:source)")
    Base.include_string(VSmartMOMForward,method,"isolated_fast_replay.jl")
end
records=[]
for mode in modes
    f=mode==:matrix ? ()->evaluator(state) :
        ()->VSmartMOMForward.evaluate_oco_forward_fast(evaluator,state)
    println("WARM mode=$mode"); flush(stdout)
    warm=quiet(f); CUDA.synchronize()
    println("WARM_DONE mode=$mode timing=$(warm.timing) archived_y_error=$(maximum(abs,warm.measurement-expected_y)) archived_K_error=$(maximum(abs,warm.jacobian-expected_K)) noise_scaled_error=$(maximum(abs,(warm.measurement-expected_y)./noise))"); flush(stdout)
    samples=[]
    for trial in 1:3
        GC.gc(); CUDA.reclaim(); CUDA.synchronize()
        t=@timed begin
            evaluation=quiet(f)
            CUDA.synchronize()
            evaluation
        end
        result=t.value
        r=Dict("trial"=>trial,"seconds"=>t.time,"host_bytes"=>t.bytes,
            "gctime"=>t.gctime,"timing"=>Dict(string(k)=>v for (k,v) in pairs(result.timing)),
            "max_y_error"=>maximum(abs,result.measurement-expected_y),
            "max_K_error"=>maximum(abs,result.jacobian-expected_K),
            "max_noise_scaled_error"=>maximum(abs,(result.measurement-expected_y)./noise))
        push!(samples,r)
        println("SAMPLE mode=$mode ",r);flush(stdout)
        if trial==3
            JLD2.jldsave(joinpath(ENV["REPLAY_OUTPUT"],"$mode-output.jld2");
                y=result.measurement,K=result.jacobian,state)
        end
    end
    push!(records,Dict("mode"=>string(mode),"samples"=>samples))
    open(io->TOML.print(io,Dict("package"=>pathof(vSmartMOM),"saved_state"=>saved,
        "records"=>records)),joinpath(ENV["REPLAY_OUTPUT"],"timings.toml"),"w")
end

# Compare the two current implementations after timing. The archived result
# remains a reported reference: the precision repair intentionally changes its
# former Float32/Float64 mixture. Coordinate units differ across K columns, so
# use a norm for each column rather than one norm dominated by large columns.
if :source in modes
    matrix = JLD2.load(joinpath(ENV["REPLAY_OUTPUT"],"matrix-output.jld2"))
    source = JLD2.load(joinpath(ENV["REPLAY_OUTPUT"],"source-output.jld2"))
    noise_error = maximum(abs,(source["y"]-matrix["y"])./noise)
    column_errors = [norm(source["K"][:,p]-matrix["K"][:,p]) /
                     max(norm(matrix["K"][:,p]),eps(Float64))
                     for p in axes(matrix["K"],2)]
    @assert noise_error < 1e-3 "Source/matrix radiance disagreement exceeds 0.001 noise units"
    @assert maximum(column_errors) < 3e-4 "Source/matrix Jacobian column disagreement"
    println("PARITY noise_units=$noise_error max_column_relative_l2=$(maximum(column_errors))")
end
