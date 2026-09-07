# Float64 preparation + RT evaluated at the exact Float32 solve-grid nodes.
# Changing a grid label AFTER RT is only an instrument test; this experiment
# rebuilds all spectral optics and sources at the supplied coordinates.
source = read(joinpath(@__DIR__,"precision_frozen.jl"),String)
@assert endswith(source,"main()\n")
Base.include_string(Main,chop(source;tail=length("main()\n")),"frozen_helpers.jl")

function grid_main()
    println("PREPARE grid-control evaluators"); flush(stdout)
    lo = quiet(() -> V.OCOForwardEvaluator(;architecture=:GPU,float_type=Float32,nstreams=9))
    hi = quiet(() -> V.OCOForwardEvaluator(;architecture=:GPU,float_type=Float64,nstreams=9))
    V.set_fixed_upper_co2_ppm!(hi,400.0)
    n = hi.base_parameters.numerics
    options = (; (f=>getfield(n,f) for f in fieldnames(typeof(n)))...)
    hi.base_parameters.numerics = typeof(n)(;merge(options,
        (dτ_max_threshold=Float64(Float32(0.001)),dτ_min_floor=Float64(1024eps(Float32))))...)
    original_grid = copy(hi.base_parameters.spec_bands[1])
    for band in eachindex(hi.base_parameters.spec_bands)
        hi.base_parameters.spec_bands[band] = Float64.(lo.base_parameters.spec_bands[band])
    end
    records = Dict[]
    for mode in (:reference,:optimized)
        state = JLD2.load(joinpath(output_dir,"state035_corrected_siffalse-$mode.jld2"),"state")
        println("GRID CONTROL state=$mode"); flush(stdout)
        setup = quiet(() -> prepare(hi,state))
        @assert setup.params.spec_bands[1] == Float64.(lo.base_parameters.spec_bands[1])
        result = quiet(() -> rt_run_lin(setup.model,setup.planned;
            i_band=1,sources=setup.sources,jacobian_basis=:local,jacobian_adding=:source))
        observed = observe(result,setup,hi)
        file = "precision64-grid32-$mode.jld2"
        JLD2.jldsave(joinpath(output_dir,file);state,observed...)
        cache = C.build_m_invariant_cache_lin(1,setup.model,setup.planned.base;
            active_layout=band_layout(setup.planned.plan,1))
        local_cache = C.build_local_jacobian_cache(1,setup.model,setup.planned.base,cache)
        counts = [last(C.get_dtau_ndoubl(C.CoreScatteringOpticalProperties(
            l.τ,l.ϖ,nothing,nothing),setup.model.quad_points;
            dτ_max_threshold=setup.model.numerics.dτ_max_threshold,
            dτ_min_floor=setup.model.numerics.dτ_min_floor)) for l in local_cache.layers]
        push!(records,Dict("state"=>string(mode),"file"=>file,
            "m_used"=>C._LAST_FOURIER_M_USED[],"ndoubl"=>counts,
            "sha256"=>bytes2hex(sha256(read(joinpath(output_dir,file))))))
    end
    grid = hi.base_parameters.spec_bands[1]
    open(io -> TOML.print(io,Dict("records"=>records,
        "grid_max_difference_cm_inverse"=>maximum(abs,grid-original_grid),
        "grid_max_difference_nm"=>maximum(abs,1e7./grid-1e7./original_grid),
        "n_spectral_points"=>length(grid))),joinpath(output_dir,"precision-grid.toml"),"w")
    println("GRID CONTROL COMPLETE")
end
grid_main()
