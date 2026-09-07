# Recompute the actual instrument output from saved high-resolution Stokes.
# No RT rerun or change to the live study. Run from test/ with STUDY_ROOT set.
using JLD2, NCDatasets, LinearAlgebra, Statistics, TOML, SHA
include(joinpath(ENV["STUDY_ROOT"],"inversion/instrument/SyntheticOCO2.jl"))
using .SyntheticOCO2
# Read the actual grid definition without loading RT or constructing models.
common_path = joinpath(ENV["STUDY_ROOT"],"scripts/common.jl")
common_source = read(common_path,String)
start = findfirst("function surface_basis_grids(",common_source).start
stop = findnext("\nend",common_source,start).stop
Base.include_string(Main,common_source[start:stop],"study_grid_definition.jl")
root = ARGS[1]
spec = BAND_SPECS[1]
coefficients = read_representative_coefficients(joinpath(ENV["STUDY_ROOT"],
    "inversion/instrument/representative_stokes_coefficients.nc"))[spec.name]
target = synthetic_grid(spec)
# Include every detector center exactly, plus four intermediate centers.
dense = vcat([target[i] + (target[i+1]-target[i])*j/5
              for i in 1:length(target)-1 for j in 0:4],last(target))
@assert dense[1:5:end] == target
noise = sqrt.(JLD2.load(joinpath(root,
    "state035_corrected_siffalse-reference.jld2"),"variance")[1:length(target)])
grids = Dict(ft => 1e7 ./ Float64.(surface_basis_grids(ft)[1]) for ft in (Float32,Float64))

function output(file,grid)
    saved = JLD2.load(joinpath(root,file*".jld2"))
    stokes = saved["stokes"]
    raw = per_wavenumber_to_per_wavelength(grid,project_oco_analyzer(stokes,coefficients))
    sampled = process_stokes_spectrum(grid,stokes,coefficients,spec)
    @assert sampled == saved["y"] # Exact proof that earlier gaps were post-convolution.
    convolved = gaussian_convolve_resample(grid,raw,dense,spec.fwhm_nm)
    @assert convolved[1:5:end] == sampled
    constant = gaussian_convolve_resample(grid,ones(length(grid)),target,spec.fwhm_nm)
    @assert constant == ones(length(target))
    (;file,grid,raw,sampled,convolved,K=saved["K"],state=saved["state"])
end
function metrics(a,b,coordinates; detector=false)
    delta = b-a
    imax = argmax(abs.(delta))
    result = Dict("n_samples"=>length(a),"max_absolute"=>maximum(abs,delta),
        "rms_absolute"=>sqrt(mean(abs2,delta)),"mean_signed"=>mean(delta),
        "max_relative_to_reference_peak"=>maximum(abs,delta)/maximum(abs,a),
        "relative_l2"=>norm(delta)/norm(a),"worst_wavelength_nm"=>coordinates[imax])
    if detector
        normalized = delta./noise
        result["max_noise_sigma"] = maximum(abs,normalized)
        result["rms_noise_sigma"] = sqrt(mean(abs2,normalized))
        result["samples_above_0p01_noise_sigma"] = count(>(0.01),abs.(normalized))
        result["passes_0p01_noise_sigma"] = maximum(abs,normalized) <= 0.01
    end
    result
end

# Independently accumulate the same sampled trapezoidal Gaussian integral in
# 256-bit arithmetic. This checks arithmetic, not unresolved spectral features.
function big_convolution(grid,raw,center)
    order = sortperm(grid)
    x,y = BigFloat.(grid[order]),BigFloat.(raw[order])
    sigma = BigFloat(spec.fwhm_nm)/(2sqrt(2log(BigFloat(2))))
    # Keep exactly the same finite support membership as the production map.
    radius = 6fwhm_to_sigma(spec.fwhm_nm)
    selected = findall(t -> abs(t-center)<=radius,grid[order])
    numerator,denominator = zero(BigFloat),zero(BigFloat)
    for i in selected
        width = i == 1 ? (x[2]-x[1])/2 : i == length(x) ?
            (x[end]-x[end-1])/2 : (x[i+1]-x[i-1])/2
        weight = width*exp(-((x[i]-BigFloat(center))/sigma)^2/2)
        numerator += weight*y[i]
        denominator += weight
    end
    Float64(numerator/denominator)
end

function main()
    records, arithmetic, remainder = Dict[], Dict[], Dict[]
    outputs = Dict()
    for state in ("reference","optimized")
        println(stderr,"CONVOLUTION CHECK state=$state"); flush(stderr)
        for (label,file,FT) in (
            ("p32r32","frozen-$state-prep32-rt32",Float32),
            ("p32r64","frozen-$state-prep32-rt64",Float32),
            ("p64r32","frozen-$state-prep64-rt32",Float64),
            ("p64r64","frozen-$state-prep64-rt64",Float64),
            ("grid32","precision64-grid32-$state",Float32))
            outputs[(state,label)] = output(file,grids[FT])
        end
        for (label,left,right) in (
            ("RT_only_preparation32","p32r32","p32r64"),
            ("RT_only_preparation64","p64r32","p64r64"),
            ("preparation_matched_grid_RT64","p32r64","grid32"),
            ("full_workflow_matched_grid","p32r32","grid32"),
            ("full_workflow_native_grids","p32r32","p64r64"))
            a,b = outputs[(state,left)],outputs[(state,right)]
            @assert a.state == b.state
            record = Dict("state"=>state,"comparison"=>label,
                "convolved_dense"=>metrics(a.convolved,b.convolved,dense),
                "convolved_detector"=>metrics(a.sampled,b.sampled,target;detector=true))
            if a.grid == b.grid
                inside = findall(x -> first(target)<=x<=last(target),a.grid)
                record["before_convolution"] = metrics(a.raw[inside],b.raw[inside],a.grid[inside])
                convolved_delta = gaussian_convolve_resample(a.grid,b.raw-a.raw,target,spec.fwhm_nm)
                linear_error = maximum(abs,convolved_delta-(b.sampled-a.sampled))
                @assert linear_error <= 1e-12*max(maximum(abs,a.sampled),maximum(abs,b.sampled))
                record["convolution_linearity_max_absolute"] = linear_error
            end
            push!(records,record)
        end
        # Endpoints, representative interior points, and the two worst residuals.
        base = outputs[(state,"p32r32")]
        worst = [argmax(abs.((outputs[(state,k)].sampled-base.sampled)./noise))
                 for k in ("p32r64","grid32")]
        indices = unique(vcat(1,length(target),collect(1:100:length(target)),worst))
        for label in ("p32r32","p32r64","p64r32","p64r64","grid32")
            a = outputs[(state,label)]
            differences = [abs(big_convolution(a.grid,a.raw,target[i])-a.sampled[i])/noise[i]
                           for i in indices]
            @assert maximum(differences) < 1e-9
            push!(arithmetic,Dict("state"=>state,"configuration"=>label,
                "n_checked"=>length(indices),"max_noise_sigma"=>maximum(differences)))
        end
    end
    for label in ("p32r32","p32r64","p64r32","p64r64","grid32")
        a,b = outputs[("reference",label)],outputs[("optimized",label)]
        prediction = a.K*(b.state-a.state)
        residual = b.sampled-a.sampled-prediction
        push!(remainder,Dict("configuration"=>label,
            "max_noise_sigma"=>maximum(abs,residual./noise),
            "rms_noise_sigma"=>sqrt(mean(abs2,residual./noise))))
    end
    a,b,c = [outputs[("reference",label)] for label in ("p32r32","p32r64","grid32")]
    peak = maximum(abs,a.raw)
    JLD2.jldsave(joinpath(root,"convolution-stages.jld2");
        raw_wavelength_nm=a.grid,dense_wavelength_nm=dense,detector_wavelength_nm=target,
        raw_rt_ppm=(b.raw-a.raw).*1e6./peak,raw_matched_ppm=(c.raw-a.raw).*1e6./peak,
        dense_rt_ppm=(b.convolved-a.convolved).*1e6./peak,
        dense_matched_ppm=(c.convolved-a.convolved).*1e6./peak,
        detector_rt_noise=(b.sampled-a.sampled)./noise,
        detector_matched_noise=(c.sampled-a.sampled)./noise)
    input_hashes = Dict(a.file*".jld2"=>bytes2hex(sha256(read(joinpath(root,a.file*".jld2"))))
                        for a in values(outputs))
    result = Dict("band"=>string(spec.name),"fwhm_nm"=>spec.fwhm_nm,
        "support_sigma"=>6,"detector_spacing_nm"=>spec.sampling_interval_nm,
        "detector_samples"=>length(target),"dense_samples"=>length(dense),
        "exact_saved_measurement_checks"=>length(outputs),
        "records"=>records,"bigfloat_arithmetic_checks"=>arithmetic,
        "post_convolution_jacobian_remainders"=>remainder,"input_sha256"=>input_hashes)
    open(io -> TOML.print(io,result),joinpath(root,"convolution-verified.toml"),"w")
    println(stderr,"CONVOLUTION CHECK COMPLETE")
end
setprecision(main,BigFloat,256)
