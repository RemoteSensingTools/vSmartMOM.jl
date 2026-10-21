# Diagnostic only: freeze the production optical/source boundary, then change RT FT.
# Run from test/ with the same study environment as precision_matched.jl.
# No study files or production methods are modified. The cloned driver uses the
# production kernels; an identity run must reproduce the public API exactly.
using vSmartMOM, CUDA, NCDatasets, JLD2, LinearAlgebra, Logging, TOML, SHA
push!(LOAD_PATH, pkgdir(vSmartMOM))
CUDA.device!(0)
CUDA.allowscalar(false)
BLAS.set_num_threads(1)
const C = vSmartMOM.CoreRT
const study = ENV["STUDY_ROOT"]
include(joinpath(study, "inversion/OptimalEstimation.jl"))
include(joinpath(study, "inversion/VSmartMOMForward.jl"))
const V = VSmartMOMForward
quiet(f) = redirect_stdout(() -> with_logger(f, NullLogger()), devnull)
const output_dir = ENV["REPLAY_OUTPUT"]
const frozen = Dict{Any,Any}()
const target_precision = Ref{DataType}(Float32)

# Explicit supported payloads: fail on an unrecognized new boundary type.
# Memoization preserves shared phase bases and avoids spectral copies per layer.
function precision_cast(x, FT, memo=IdDict())
    x isa AbstractFloat && return FT(x)
    x isa Union{Integer,Nothing,Symbol} && return x
    haskey(memo, x) && return memo[x]
    result = if x isa AbstractArray{<:AbstractFloat}
        eltype(x) === FT ? x : FT.(x)
    elseif x isa AbstractArray || x isa Tuple
        map(v -> precision_cast(v, FT, memo), x)
    elseif x isa C.PreparedSolarBeam
        f = precision_cast(x.F₀, FT, memo)
        C.PreparedSolarBeam{FT,typeof(f)}(f)
    elseif x isa C.PreparedSurfaceSIF
        f = precision_cast(x.SIF₀, FT, memo)
        df = precision_cast(x.SIḞ₀, FT, memo)
        C.PreparedSurfaceSIF{FT,typeof(f),typeof(df)}(f, df, x.n_parameters)
    elseif x isa C.NoSource
        x
    elseif x isa Union{C.QuadPoints,C.LambertianSurfaceLegendre,
                       C.CoreScatteringOpticalProperties,
                       C.CoreScatteringOpticalPropertiesLin,
                       C.LocalOpticalJacobian,C.SourceSet}
        Base.typename(typeof(x)).wrapper(
            (precision_cast(getfield(x, n), FT, memo) for n in fieldnames(typeof(x)))...)
    else
        error("Unsupported frozen-boundary payload: $(typeof(x))")
    end
    memo[x] = result
    result
end
frozen_input(f, key) = precision_cast(get!(f, frozen, key), target_precision[])

# Extract only the final full driver. Each replacement must match exactly once.
driver_path = joinpath(pkgdir(vSmartMOM), "src/CoreRT/rt_run_lin.jl")
driver = read(driver_path, String)
driver = driver[findfirst("# Full multiple scattering\n", driver).start:end]
function replace_once(source, old, new)
    @assert length(findall(old, source)) == 1 old
    replace(source, old => new)
end
driver = replace_once(driver, "function rt_run(", "function precision_frozen_run(")
driver = replace_once(driver, "(; obs_alt, sza, vza, vaz) = model.obs_geom",
    "(; obs_alt, sza, vza, vaz) = model.obs_geom\n" *
    "    sza, vza, vaz = Main.frozen_input(() -> (sza,vza,vaz), :angles)\n" *
    "    frozen_quad = Main.frozen_input(() -> model.quad_points, :quad)")
driver = replace(driver, "= model.quad_points #" => "= frozen_quad #",
    "quad_points = model.quad_points\n" => "quad_points = frozen_quad\n")
driver = replace_once(driver, "brdf = get_surface(model, iBand)",
    "brdf = Main.frozen_input(() -> get_surface(model, iBand), :surface)")
driver = replace_once(driver, "FT = eltype(sza)", "FT = Main.target_precision[]")
driver = replace_once(driver,
    "prepare_sources(effective_sources, FT, pol_type.n, nSpec, arr_type)",
    "Main.frozen_input(() -> prepare_sources(effective_sources, float_type(model), pol_type.n, nSpec, arr_type), :sources)")
driver = replace_once(driver, "weight = m == 0 ? FT(0.5/π) : FT(1.0/π)",
    "weight = Main.frozen_input(() -> (m == 0 ? float_type(model)(0.5/π) : float_type(model)(1.0/π)), (:weight,m))")
driver = replace_once(driver,
    "(local_basis ? construct_local_optical_jacobians : constructCoreOpticalProperties)(\n                RS_type, iBand, m, model, lin_model, m_invariant_cache)",
    "Main.frozen_input(() -> (local_basis ? construct_local_optical_jacobians : constructCoreOpticalProperties)(\n                RS_type, iBand, m, model, lin_model, m_invariant_cache), (:optics,m))")
Base.include_string(C, driver, "precision_frozen_driver.jl")

function prepare(evaluator, state)
    params = copy_parameters(evaluator.base_parameters; share_luts=true)
    physical = V.apply_retrieval_state!(params, state, evaluator.tau_ref_scale;
        fixed_upper_co2_vmr=evaluator.fixed_upper_co2_vmr)
    model, planned = model_from_parameters(OCO_RRS_synth(), params; external_solar=true)
    sources = V.RRSXCO2Common.sources_for_band(params, 1;
        SIF760=physical.SIF760, mSIF=physical.mSIF, solar_T=evaluator.solar_transmission)
    (; params, physical, model, planned, sources)
end

function observe(result, setup, evaluator)
    stokes, local_K = V._canonical_toa(result)
    global_K = Array(globalize_jacobian(local_K, result.layout))
    jac = Array(@view global_K[1,:,:,:])
    spec = V.BAND_SPECS[1]
    wavelength = 1e7 ./ Float64.(setup.params.spec_bands[1])
    y = V.process_stokes_spectrum(wavelength, stokes, evaluator.coefficients[spec.name], spec)
    K = V.process_stokes_jacobian(wavelength, jac, evaluator.coefficients[spec.name], spec)
    columns = first(V.SURFACE_RANGE):first(V.SURFACE_RANGE)+2
    K[:,columns] = K[:,columns] * Float64.(V.surface_coefficient_transform(setup.params,1))
    K[:,V.LOG_AOD_RANGE] .*= reshape(Float64.(setup.physical.tau_ref),1,:)
    K[:,V.LOG_HEIGHT_RANGE] .*= reshape(Float64.(setup.physical.aerosol_height),1,:)
    @assert all(isfinite,y) && all(isfinite,K)
    (; y, K, stokes)
end

function frozen_solve(setup, FT)
    target_precision[] = FT
    (;model,planned,sources) = setup
    C.precision_frozen_run(vSmartMOM.InelasticScattering.noRS{FT}(),
        model, planned.base, C.n_aerosols(model), size(planned.base.τ̇_abs[1],1),
        C.surface_parameter_count(C.get_surface(model,1)), 1;
        sources, active_layout=band_layout(planned.plan,1),
        jacobian_basis=:local, jacobian_adding=:source)
end

function main()
    records = Dict[]
    evaluators = Dict()
    for FT in (Float32,Float64)
        println("PREPARE evaluator $FT"); flush(stdout)
        evaluator = quiet(() -> V.OCOForwardEvaluator(;architecture=:GPU,float_type=FT,nstreams=9))
        V.set_fixed_upper_co2_ppm!(evaluator,400.0)
        n = evaluator.base_parameters.numerics
        options = (; (f=>getfield(n,f) for f in fieldnames(typeof(n)))...)
        evaluator.base_parameters.numerics = typeof(n)(;merge(options,
            (dτ_max_threshold=FT(Float32(0.001)),dτ_min_floor=FT(1024eps(Float32))))...)
        evaluators[FT] = evaluator
    end
    for state_mode in (:reference,:optimized)
        state = JLD2.load(joinpath(output_dir,"state035_corrected_siffalse-$state_mode.jld2"),"state")
        for prepFT in (Float32,Float64)
            println("BUILD state=$state_mode preparation=$prepFT"); flush(stdout)
            setup = quiet(() -> prepare(evaluators[prepFT],state))
            empty!(frozen)
            native = quiet(() -> rt_run_lin(setup.model,setup.planned;
                i_band=1,sources=setup.sources,jacobian_basis=:local,jacobian_adding=:source))
            for solveFT in (prepFT, prepFT === Float32 ? Float64 : Float32)
                println("SOLVE state=$state_mode preparation=$prepFT RT=$solveFT"); flush(stdout)
                result = quiet(() -> frozen_solve(setup,solveFT))
                CUDA.synchronize()
                identity = prepFT === solveFT
                if identity
                    @assert Array(result.toa) == Array(native.toa)
                    @assert Array(result.toa_jacobian) == Array(native.toa_jacobian)
                end
                observed = observe(result,setup,evaluators[prepFT])
                filename = "frozen-$state_mode-prep$(sizeof(prepFT)*8)-rt$(sizeof(solveFT)*8).jld2"
                JLD2.jldsave(joinpath(output_dir,filename);state,observed...)
                # Compare instrument-coordinate preparation using the SAME Stokes array.
                instrument_y = Dict{String,Any}()
                for instFT in (Float32,Float64)
                    e = evaluators[instFT]
                    wavelength = 1e7 ./ Float64.(e.base_parameters.spec_bands[1])
                    spec = V.BAND_SPECS[1]
                    instrument_y[string(instFT)] = V.process_stokes_spectrum(
                        wavelength,observed.stokes,e.coefficients[spec.name],spec)
                end
                JLD2.jldsave(joinpath(output_dir,replace(filename,".jld2"=>"-instrument.jld2"));instrument_y)
                counts = [last(C.get_dtau_ndoubl(v,precision_cast(frozen[:quad],solveFT);
                    dτ_max_threshold=setup.model.numerics.dτ_max_threshold,
                    dτ_min_floor=setup.model.numerics.dτ_min_floor))
                    for v in precision_cast(frozen[(:optics,0)][1],solveFT)]
                push!(records,Dict("state"=>string(state_mode),"preparation"=>string(prepFT),
                    "rt"=>string(solveFT),"identity_checked"=>identity,
                    "m_used"=>C._LAST_FOURIER_M_USED[],"ndoubl"=>counts,
                    "file"=>filename,"sha256"=>bytes2hex(sha256(read(joinpath(output_dir,filename))))))
                open(io -> TOML.print(io,Dict("records"=>records,
                    "driver_sha256"=>bytes2hex(sha256(read(driver_path))))),
                    joinpath(output_dir,"precision-frozen.toml"),"w")
                GC.gc(); CUDA.reclaim()
            end
            empty!(frozen)
            GC.gc(); CUDA.reclaim()
        end
    end
    println("FROZEN PRECISION COMPLETE")
end
main()
