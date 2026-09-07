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
root = joinpath(study,"bottom_layer_XCO2_retrievals/retrievals_acos_mapped_tapered_vertical_correlation_nosif")
results = Dict[]
for (scene,kind,sif) in ((1,"corrected",false),(1,"uncorrected",false),
                          (35,"corrected",false),(35,"uncorrected",false),
                          (35,"corrected",true),(35,"uncorrected",true))
    id = "state$(lpad(scene,3,'0'))_$(kind)_sif$(sif)"
    path = joinpath(root,kind,"retrieval_state$(lpad(scene,3,'0'))_perturbation10.nc")
    data = NCDataset(path) do ds
        settings = OESettings(; (k=>ds.attrib[string(k)] for k in
            (:convergence_threshold,:maximum_iterations,:maximum_divergences,
             :maximum_band_chi_squared,:initial_gamma))...)
        (; y=Array(ds["measurement_perturbed"][:]),variance=Array(ds["Se_diagonal"][:]),
           xa=Array(ds["a_priori_state"][:]),Sa=Array(ds["a_priori_covariance"][:,:]),
           initial=Array(ds["state_at_trial"][:,1]),archived=Array(ds["final_state"][:]),
           settings)
    end
    y = copy(data.y)
    # Controlled SIF extension of the archived observation, using the physical
    # reference at the archived final atmospheric state. This is synthetic,
    # not a replay of a released SIF truth campaign.
    if sif
        a = copy(data.archived); a[29:30] .= 0
        b = copy(a)
        shape = VSmartMOMForward.RRSXCO2Common.campaign_sif_state()
        b[29:30] .= (shape.SIF760,shape.mSIF)
        off = quiet(()->VSmartMOMForward.evaluate_reference(evaluator,a))
        on = quiet(()->VSmartMOMForward.evaluate_reference(evaluator,b))
        y .+= on.measurement .- off.measurement
    end
    pair = Dict()
    for mode in (:reference,:optimized)
        fn = getproperty(VSmartMOMForward,Symbol("evaluate_$mode"))
        println("START $id $mode");flush(stdout)
        callback = r -> begin
            println("TRIAL $id $mode n=$(r.trial) accepted=$(r.accepted) cost=$(r.total_cost) step=$(r.d_sigma_sq_scaled)");flush(stdout)
        end
        GC.gc();CUDA.reclaim()
        timed = @timed solve_optimal_estimation(x->quiet(()->fn(evaluator,x)),y,
            data.variance,data.xa,data.Sa;initial_state=data.initial,
            settings=data.settings,record_callback=callback)
        r = timed.value
        pair[mode] = r
        JLD2.jldsave(joinpath(output_dir,"$id-$mode.jld2");
            state=r.final_state,y=r.final_measurement,K=r.final_jacobian,
            posterior=r.posterior_covariance,averaging_kernel=r.averaging_kernel,
            states=hcat((v.state for v in r.records)...),
            costs=[v.total_cost for v in r.records],accepted=[v.accepted for v in r.records],
            observation=y,variance=data.variance,xa=data.xa,Sa=data.Sa)
        residual = r.final_measurement-y
        delta = r.final_state-data.xa
        record = Dict("case"=>id,"mode"=>string(mode),"input"=>path,
            "input_sha256"=>bytes2hex(sha256(read(path))),"synthetic_sif"=>sif,
            "seconds"=>timed.time,"outcome"=>r.outcome,"converged"=>r.converged,
            "trials"=>length(r.records),"rejected"=>count(v->!v.accepted,r.records),
            "cost"=>sum(abs2,residual./sqrt.(data.variance))+dot(delta,data.Sa\delta),
            "xco2_ppm"=>VSmartMOMForward.column_averaged_co2_ppm(evaluator,r.final_state),
            "chi_squared"=>r.final_band_chi_squared)
        push!(results,record)
        open(io->TOML.print(io,Dict("records"=>results)),joinpath(output_dir,"results.toml"),"w")
        println("DONE ",record);flush(stdout)
    end
    a,b = pair[:reference],pair[:optimized]
    println("COMPARE $id max_state_sigma=",maximum(abs.(a.final_state-b.final_state)./sqrt.(diag(data.Sa))),
        " max_noise=",maximum(abs.(a.final_measurement-b.final_measurement)./sqrt.(data.variance)))
    @assert isequal(template_before,template_state(evaluator.base_parameters))
    @assert hashes_before == table_hashes(evaluator.base_parameters)
end
println("COMPLETE template and LUTs unchanged")
