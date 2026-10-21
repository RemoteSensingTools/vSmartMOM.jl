# Full three-band retrievals at Float32 and Float64, with matched elemental
# controls. A second Float64 solve uses exact promoted Float32 spectral grids.
# All runs retain the archived observation, prior, initial state, and stopping rule.
source = read(joinpath(@__DIR__,"replay.jl"),String)
stop = findfirst("evaluator = quiet(()->OCOForwardEvaluator",source).start
Base.include_string(Main,source[1:prevind(source,stop)],"retrieval_precision_helpers.jl")

function settings_matched!(e)
    n=e.base_parameters.numerics
    FT=e.base_parameters.float_type
    options=(;(f=>getfield(n,f) for f in fieldnames(typeof(n)))...)
    e.base_parameters.numerics=typeof(n)(;merge(options,
        (dτ_max_threshold=FT(Float32(0.001)),dτ_min_floor=FT(1024eps(Float32))))...)
    VSmartMOMForward.set_fixed_upper_co2_ppm!(e,400.0)
end
function xco2_gradient(e,x)
    gradient=zeros(length(x))
    for j in 1:13
        h=j==1 ? 1e-3 : 1e-6
        a,b=copy(x),copy(x)
        a[j]-=h; b[j]+=h
        gradient[j]=(column_averaged_co2_ppm(e,b)-column_averaged_co2_ppm(e,a))/(2h)
    end
    gradient
end
function projection(e,baseline,delta)
    K=baseline["K"]; variance=baseline["variance"]; Sa=baseline["Sa"]
    scale=sqrt.(diag(Sa)); B=K.*reshape(scale,1,:)./sqrt.(variance)
    correlation=Sa./(scale*scale')
    dz=-(B'B+inv(Symmetric(correlation)))\(B'*(delta./sqrt.(variance)))
    dx=scale.*dz
    gradient=xco2_gradient(e,baseline["state"])
    Dict("xco2_shift_ppm"=>dot(gradient,dx),"state_shift"=>dx,
        "max_state_prior_sigma"=>maximum(abs,dz),
        "psurf_shift_hpa"=>dx[1],"active_layer_co2_shift_ppm"=>dx[2:13].*1e6)
end
function main_precision_retrieval()
    path=joinpath(study,"bottom_layer_XCO2_retrievals/retrievals_acos_mapped_tapered_vertical_correlation_nosif/corrected/retrieval_state035_perturbation10.nc")
    data=NCDataset(path) do ds
        settings=OESettings(; (k=>ds.attrib[string(k)] for k in
            (:convergence_threshold,:maximum_iterations,:maximum_divergences,
             :maximum_band_chi_squared,:initial_gamma))...)
        (;y=Array(ds["measurement_perturbed"][:]),variance=Array(ds["Se_diagonal"][:]),
          xa=Array(ds["a_priori_state"][:]),Sa=Array(ds["a_priori_covariance"][:,:]),
          initial=Array(ds["state_at_trial"][:,1]),settings)
    end
    baseline=JLD2.load(joinpath(output_dir,"state035_corrected_siffalse-optimized.jld2"))
    @assert baseline["observation"]==data.y && baseline["variance"]==data.variance
    @assert baseline["xa"]==data.xa && baseline["Sa"]==data.Sa
    println("PREPARE precision-retrieval evaluators");flush(stdout)
    lo=quiet(()->OCOForwardEvaluator(;architecture=:GPU,float_type=Float32,nstreams=9))
    hi=quiet(()->OCOForwardEvaluator(;architecture=:GPU,float_type=Float64,nstreams=9))
    settings_matched!(lo);settings_matched!(hi)
    # Local projections of the isolated O2 differences, keeping CO2-band delta=0.
    projections=Dict[]
    reference=JLD2.load(joinpath(output_dir,"frozen-optimized-prep32-rt32.jld2"),"y")
    for (label,file) in (("RT_only_O2","frozen-optimized-prep32-rt64"),
                         ("full_O2_matched_grid","precision64-grid32-optimized"),
                         ("full_O2_native_grid","frozen-optimized-prep64-rt64"))
        y=JLD2.load(joinpath(output_dir,file*".jld2"),"y")
        delta=zeros(length(data.y));delta[1:length(y)].=y-reference
        p=projection(hi,baseline,delta);p["configuration"]=label
        push!(projections,p)
        println("PROJECTION $label xco2_ppm=",p["xco2_shift_ppm"]);flush(stdout)
    end
    records=Dict[]
    for (label,e) in (("Float32",lo),("Float64_native_grid",hi),("Float64_matched_grid",hi))
        if label=="Float64_matched_grid"
            for band in eachindex(hi.base_parameters.spec_bands)
                hi.base_parameters.spec_bands[band]=Float64.(lo.base_parameters.spec_bands[band])
            end
        end
        println("FIXED STATE $label");flush(stdout)
        fixed=quiet(()->VSmartMOMForward.evaluate_optimized(e,baseline["state"]))
        if label=="Float32"
            @assert fixed.measurement==baseline["y"] && fixed.jacobian==baseline["K"]
        end
        p=projection(hi,baseline,fixed.measurement-baseline["y"])
        p["configuration"]=label*"_all_bands";push!(projections,p)
        println("START PRECISION RETRIEVAL $label");flush(stdout)
        callback=r->begin
            println("TRIAL $label n=$(r.trial) accepted=$(r.accepted) cost=$(r.total_cost) step=$(r.d_sigma_sq_scaled)");flush(stdout)
        end
        GC.gc();CUDA.reclaim()
        timed=@timed solve_optimal_estimation(x->quiet(()->VSmartMOMForward.evaluate_optimized(e,x)),
            data.y,data.variance,data.xa,data.Sa;initial_state=data.initial,settings=data.settings,
            record_callback=callback)
        r=timed.value
        label=="Float32" && @assert r.final_state==baseline["state"]
        gradient=xco2_gradient(hi,r.final_state)
        dx=r.final_state-baseline["state"]
        file="precision-retrieval-$label.jld2"
        JLD2.jldsave(joinpath(output_dir,file);state=r.final_state,y=r.final_measurement,
            K=r.final_jacobian,posterior=r.posterior_covariance,averaging_kernel=r.averaging_kernel,
            states=hcat((v.state for v in r.records)...),
            costs=[v.total_cost for v in r.records],accepted=[v.accepted for v in r.records],
            fixed_y=fixed.measurement,fixed_K=fixed.jacobian,
            observation=data.y,variance=data.variance,xa=data.xa,Sa=data.Sa)
        record=Dict("configuration"=>label,"converged"=>r.converged,"outcome"=>r.outcome,
            "trials"=>length(r.records),"accepted"=>[v.accepted for v in r.records],
            "cost"=>sum(abs2,(r.final_measurement-data.y)./sqrt.(data.variance))+
                dot(r.final_state-data.xa,data.Sa\(r.final_state-data.xa)),
            "xco2_native_ppm"=>column_averaged_co2_ppm(e,r.final_state),
            "xco2_common64_ppm"=>column_averaged_co2_ppm(hi,r.final_state),
            "xco2_shift_common64_ppm"=>column_averaged_co2_ppm(hi,r.final_state)-
                column_averaged_co2_ppm(hi,baseline["state"]),
            "posterior_xco2_sigma_ppm"=>sqrt(dot(gradient,r.posterior_covariance*gradient)),
            "state_shift"=>dx,"max_state_prior_sigma"=>maximum(abs,dx./sqrt.(diag(data.Sa))),
            "psurf_shift_hpa"=>dx[1],"active_layer_co2_shift_ppm"=>dx[2:13].*1e6,
            "chi_squared"=>r.final_band_chi_squared,"seconds"=>timed.time,
            "file"=>file,"sha256"=>bytes2hex(sha256(read(joinpath(output_dir,file)))))
        push!(records,record)
        settings=Dict(string(f)=>getfield(data.settings,f) for f in fieldnames(typeof(data.settings)))
        open(io->TOML.print(io,Dict("records"=>records,"linearized_projections"=>projections,
            "settings"=>settings,"input_sha256"=>bytes2hex(sha256(read(path))))),
            joinpath(output_dir,"precision-retrieval.toml"),"w")
        println("DONE $label xco2_shift=",record["xco2_shift_common64_ppm"]," converged=",r.converged);flush(stdout)
    end
    println("PRECISION RETRIEVAL COMPLETE")
end
main_precision_retrieval()
