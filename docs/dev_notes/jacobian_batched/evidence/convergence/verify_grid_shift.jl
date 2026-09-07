# Diagnose whether independently generated grids resemble a spectral shift.
# Fits below are offline residual projections, not added retrieval parameters.
using JLD2, NCDatasets, LinearAlgebra, Statistics, TOML
include(joinpath(ENV["STUDY_ROOT"],"inversion/instrument/SyntheticOCO2.jl"))
using .SyntheticOCO2
source=read(joinpath(ENV["STUDY_ROOT"],"scripts/common.jl"),String)
start=findfirst("function surface_basis_grids(",source).start
stop=findnext("\nend",source,start).stop
Base.include_string(Main,source[start:stop],"study_grid_definition.jl")
root=ARGS[1]
spec=BAND_SPECS[1]
target=synthetic_grid(spec)
coeff=read_representative_coefficients(joinpath(ENV["STUDY_ROOT"],
    "inversion/instrument/representative_stokes_coefficients.nc"))[spec.name]
nu32=Float64.(surface_basis_grids(Float32)[1]);nu64=surface_basis_grids(Float64)[1]
lambda32=1e7./nu32;lambda64=1e7./nu64
variance=JLD2.load(joinpath(root,"state035_corrected_siffalse-reference.jld2"),"variance")[1:length(target)]
noise=sqrt.(variance)

function coordinate_fit(delta,coordinates)
    x=coordinates.-mean(coordinates)
    A=hcat(ones(length(x)),x)
    coefficients=A\delta
    remainder=delta-A*coefficients
    Dict("max_absolute"=>maximum(abs,delta),"mean"=>mean(delta),
        "rms_about_mean"=>sqrt(mean(abs2,delta.-mean(delta))),
        "affine_offset_at_mean"=>coefficients[1],"affine_slope"=>coefficients[2],
        "affine_remainder_max"=>maximum(abs,remainder),
        "affine_remainder_rms"=>sqrt(mean(abs2,remainder)))
end
function sampled(file,grid,centers)
    stokes=JLD2.load(joinpath(root,file*".jld2"),"stokes")
    raw=per_wavenumber_to_per_wavelength(grid,project_oco_analyzer(stokes,coeff))
    gaussian_convolve_resample(grid,raw,centers,spec.fwhm_nm)
end
function describe(d)
    Dict("max_noise_sigma"=>maximum(abs,d./noise),
        "rms_noise_sigma"=>sqrt(mean(abs2,d./noise)))
end
function fit_shift(reference,grid,delta)
    h=1e-4
    derivative=(sampled(reference,grid,target.+h)-sampled(reference,grid,target.-h))/(2h)
    fine=(sampled(reference,grid,target.+h/2)-sampled(reference,grid,target.-h/2))/h
    @assert norm(fine-derivative)/norm(fine)<1e-4
    x=target.-mean(target)
    result=Dict("before"=>describe(delta),
        "derivative_half_step_relative_l2"=>norm(fine-derivative)/norm(fine))
    for affine in (false,true)
        A=affine ? hcat(fine,fine.*x) : reshape(fine,:,1)
        Aw=A./noise
        coefficients=Aw\(delta./noise)
        residual=delta-A*coefficients
        record=describe(residual)
        record["fitted_shift_nm"]=coefficients[1]
        affine && (record["fitted_stretch"]=coefficients[2])
        record["weighted_squared_difference_fraction_removed"]=
            1-sum(abs2,residual./noise)/sum(abs2,delta./noise)
        result[affine ? "shift_and_stretch" : "shift_only"]=record
    end
    result
end
function main()
    records=Dict[]
    for state in ("reference","optimized")
        native="frozen-$state-prep64-rt64"
        matched="precision64-grid32-$state"
        lo="frozen-$state-prep32-rt32"
        y_native=sampled(native,lambda64,target)
        y_matched=sampled(matched,lambda32,target)
        y_lo=sampled(lo,lambda32,target)
        @assert y_native==JLD2.load(joinpath(root,native*".jld2"),"y")
        @assert y_matched==JLD2.load(joinpath(root,matched*".jld2"),"y")
        @assert y_lo==JLD2.load(joinpath(root,lo*".jld2"),"y")
        for (label,ref,grid,delta) in (
            ("grid_only_Float64",native,lambda64,y_matched-y_native),
            ("full_native_grids",lo,lambda32,y_native-y_lo))
            record=fit_shift(ref,grid,delta)
            record["state"]=state;record["comparison"]=label
            push!(records,record)
        end
    end
    open(io->TOML.print(io,Dict("coordinate_difference_Float32_minus_Float64"=>Dict(
        "wavenumber_cm_inverse"=>coordinate_fit(nu32-nu64,nu64),
        "wavelength_nm"=>coordinate_fit(lambda32-lambda64,lambda64)),"records"=>records)),
        joinpath(root,"grid-shift-verified.toml"),"w")
    println("GRID SHIFT CHECK COMPLETE")
end
main()
