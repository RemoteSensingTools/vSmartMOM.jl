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
evaluator = quiet(()->OCOForwardEvaluator(;architecture=:GPU,float_type=Float32,nstreams=9))
VSmartMOMForward.set_fixed_upper_co2_ppm!(evaluator,400.0)
template_before = template_state(evaluator.base_parameters)
hashes_before = table_hashes(evaluator.base_parameters)

# A focused same-state comparison separates propagation error from the
# difference between the terminal states chosen by two OE trajectories.
struct MomentLogger <: AbstractLogger
    records::Vector{Dict{String,Any}}
end
Logging.min_enabled_level(::MomentLogger) = Logging.Info
Logging.shouldlog(::MomentLogger, args...) = true
Logging.catch_exceptions(::MomentLogger) = false
function Logging.handle_message(log::MomentLogger, level, message, _module, group, id, file, line; kwargs...)
    occursin("Fourier series converged",string(message)) || return
    push!(log.records,Dict(string(k)=>v for (k,v) in kwargs))
end
const case = "state035_corrected_siffalse"
summary = Dict[]
for state_mode in (:reference,:optimized)
    state = JLD2.load(joinpath(output_dir,"$case-$state_mode.jld2"),"state")
    for mode in (:reference,:optimized)
        fn = getproperty(VSmartMOMForward,Symbol("evaluate_$mode"))
        log = MomentLogger(Dict{String,Any}[])
        println("DIAGNOSTIC state=$state_mode solver=$mode");flush(stdout)
        result = redirect_stdout(devnull) do
            with_logger(log) do
                fn(evaluator,state)
            end
        end
        JLD2.jldsave(joinpath(output_dir,"diagnostic-$state_mode-$mode.jld2");
            state,y=result.measurement,K=result.jacobian)
        push!(summary,Dict("state"=>string(state_mode),"solver"=>string(mode),"moments"=>log.records))
        open(io->TOML.print(io,Dict("records"=>summary)),joinpath(output_dir,"diagnostic.toml"),"w")
        println("MOMENTS ",summary[end]);flush(stdout)
    end
end

# Hold the Fourier order policy fixed to test its role in any discrepancy
# between the nearby terminal states. This does not change the OE replay.
numerics = evaluator.base_parameters.numerics
options = (; (name=>getfield(numerics,name) for name in fieldnames(typeof(numerics)))...)
evaluator.base_parameters.numerics = typeof(numerics)(;
    merge(options,(fourier_convergence=vSmartMOM.CoreRT.AllFourierMoments(),))...)
for state_mode in (:reference,:optimized)
    state = JLD2.load(joinpath(output_dir,"$case-$state_mode.jld2"),"state")
    println("FIXED FOURIER state=$state_mode");flush(stdout)
    result = quiet(()->VSmartMOMForward.evaluate_optimized(evaluator,state))
    JLD2.jldsave(joinpath(output_dir,"diagnostic-$state_mode-allmoments.jld2");
        state,y=result.measurement,K=result.jacobian)
end
println("DIAGNOSTICS COMPLETE")
@assert isequal(template_before,template_state(evaluator.base_parameters))
@assert hashes_before == table_hashes(evaluator.base_parameters)
open(io->TOML.print(io,Dict("template_unchanged"=>true,
    "table_coefficients_unchanged"=>true,"table_sha256"=>hashes_before)),
    joinpath(output_dir,"diagnostic-isolation.toml"),"w")

# Inspect the other discrete RT choice without constructing phase matrices.
discretization = Dict[]
for state_mode in (:reference,:optimized)
    state = JLD2.load(joinpath(output_dir,"$case-$state_mode.jld2"),"state")
    params = copy_parameters(evaluator.base_parameters;share_luts=true)
    VSmartMOMForward.apply_retrieval_state!(params,state,evaluator.tau_ref_scale;
        fixed_upper_co2_vmr=evaluator.fixed_upper_co2_vmr)
    model, planned = quiet(()->model_from_parameters(OCO_RRS_synth(),params;external_solar=true))
    for ib in 1:3
        cr = vSmartMOM.CoreRT
        cache = cr.build_m_invariant_cache_lin(ib,model,planned.base;
            active_layout=band_layout(planned.plan,ib))
        local_cache = cr.build_local_jacobian_cache(ib,model,planned.base,cache)
        counts = Int[]
        scatter = Float64[]
        for layer in local_cache.layers
            optics = cr.CoreScatteringOpticalProperties(layer.τ,layer.ϖ,nothing,nothing)
            _,nd = cr.get_dtau_ndoubl(optics,model.quad_points;
                dτ_max_threshold=model.numerics.dτ_max_threshold,
                dτ_min_floor=model.numerics.dτ_min_floor)
            push!(counts,nd)
            push!(scatter,maximum(layer.τ .* layer.ϖ))
        end
        push!(discretization,Dict("state"=>string(state_mode),"band"=>ib,
            "ndoubl"=>counts,"max_scattering_depth"=>scatter))
    end
end
open(io->TOML.print(io,Dict("records"=>discretization)),
    joinpath(output_dir,"diagnostic-doubling.toml"),"w")
println("DISCRETIZATION COMPLETE")
