# Reuse the Float64 probe, retaining the Float32 elemental floor and threshold.
source = read(joinpath(@__DIR__,"precision.jl"),String)
marker = "for mode in (:reference,:optimized)\n    state=JLD2.load"
@assert occursin(marker,source)
setup = """
numerics = evaluator.base_parameters.numerics
options = (; (name=>getfield(numerics,name) for name in fieldnames(typeof(numerics)))...)
evaluator.base_parameters.numerics = typeof(numerics)(;
    merge(options,(dτ_max_threshold=Float64(Float32(0.001)),
                   dτ_min_floor=Float64(1024eps(Float32)),))...)
"""
source = replace(source,marker=>setup*"\n"*marker,
                 "precision64-"=>"precision64-matched-")
Base.include_string(Main,source,"precision_matched_probe.jl")

# Confirm the matched-floor O2 doubling counts against the Float32 trace.
records=Dict[]
for mode in (:reference,:optimized)
    state=JLD2.load(joinpath(output_dir,"state035_corrected_siffalse-$mode.jld2"),"state")
    params=copy_parameters(evaluator.base_parameters;share_luts=true)
    VSmartMOMForward.apply_retrieval_state!(params,state,evaluator.tau_ref_scale;
        fixed_upper_co2_vmr=evaluator.fixed_upper_co2_vmr)
    model,planned=quiet(()->model_from_parameters(OCO_RRS_synth(),params;external_solar=true))
    cr=vSmartMOM.CoreRT
    cache=cr.build_m_invariant_cache_lin(1,model,planned.base;
        active_layout=band_layout(planned.plan,1))
    lc=cr.build_local_jacobian_cache(1,model,planned.base,cache)
    counts=map(lc.layers) do layer
        optics=cr.CoreScatteringOpticalProperties(layer.τ,layer.ϖ,nothing,nothing)
        last(cr.get_dtau_ndoubl(optics,model.quad_points;
            dτ_max_threshold=model.numerics.dτ_max_threshold,
            dτ_min_floor=model.numerics.dτ_min_floor))
    end
    push!(records,Dict("state"=>string(mode),"ndoubl"=>counts,
        "d_tau_min_floor"=>model.numerics.dτ_min_floor,
        "d_tau_max_threshold"=>model.numerics.dτ_max_threshold))
end
open(io->TOML.print(io,Dict("records"=>records)),
    joinpath(output_dir,"precision-matched-doubling.toml"),"w")
println("MATCHED-FLOOR COMPLETE")
