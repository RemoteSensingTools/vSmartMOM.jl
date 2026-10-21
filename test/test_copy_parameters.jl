using Test, vSmartMOM, Interpolations, Logging
using vSmartMOM.CoreRT

function _copy_parameter_fixture()
    params = parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
    params.scattering_params.r_max = 3.0
    ν = range(first(params.spec_bands[1])-5,last(params.spec_bands[1])+5;length=8)
    p = range(0.0,1100.0;length=5)
    T = range(180.0,330.0;length=4)
    coefficients = [1e-26 * (1+i/8+j/5+k/4) for i in 1:8,j in 1:5,k in 1:4]
    lut = vSmartMOM.Absorption.InterpolationModel(
        interpolate(coefficients,BSpline(Linear())),7,1,ν,p,T)
    params.absorption_params.luts = [Any[lut]]
    # Cover shared objects across regular/H₂O lists without requiring H₂O
    # spectroscopy in this dry fixture. Sentinels are tested separately below.
    params.absorption_params.h2o_lut = Any[lut]
    params.absorption_params.vmr["O2"] = fill(0.21,length(params.T))
    return params,lut
end

@testset "Trial parameter copies and read-only LUTs" begin
    with_logger(NullLogger()) do
        params,lut = _copy_parameter_fixture()
        params.vaz = params.vza # deepcopy must preserve this internal alias
        independent = copy_parameters(params)
        first_trial = copy_parameters(params;share_luts=true)
        second_trial = copy_parameters(params;share_luts=true)
        @test independent.absorption_params.luts[1][1].itp.coefs !== lut.itp.coefs
        for trial in (first_trial,second_trial)
            @test trial !== params
            @test trial.absorption_params !== params.absorption_params
            @test trial.absorption_params.luts !== params.absorption_params.luts
            @test trial.absorption_params.luts[1] !== params.absorption_params.luts[1]
            @test trial.absorption_params.h2o_lut !== params.absorption_params.h2o_lut
            @test trial.absorption_params.luts[1][1].itp.coefs === lut.itp.coefs
            @test trial.absorption_params.h2o_lut[1].itp.coefs === lut.itp.coefs
            @test trial.vza === trial.vaz
            @test trial.vza !== params.vza
        end
        first_trial.p[end] -= 1
        first_trial.T[1] += 1
        first_trial.absorption_params.vmr["O2"][1] = 0.2
        first_trial.scattering_params.rt_aerosols[1].τ_ref *= 2
        first_trial.brdf[1] = LambertianSurfaceScalar(0.2)
        first_trial.spec_bands[1][1] += 1
        first_trial.vza[1] += 1
        first_trial.absorption_params.luts[1][1] = nothing
        first_trial.absorption_params.h2o_lut[1] = :disabled
        for other in (params,second_trial)
            @test other.p[end] == independent.p[end]
            @test other.T == independent.T
            @test other.absorption_params.vmr == independent.absorption_params.vmr
            @test other.scattering_params.rt_aerosols[1].τ_ref == independent.scattering_params.rt_aerosols[1].τ_ref
            @test other.brdf == independent.brdf
            @test other.spec_bands == independent.spec_bands
            @test other.vza == independent.vza
            @test other.absorption_params.luts[1][1] === lut
            @test other.absorption_params.h2o_lut[1] === lut
        end
        # Ordinary deepcopy remains independent, including after an opt-in copy.
        @test deepcopy(params).absorption_params.luts[1][1].itp.coefs !== lut.itp.coefs
        params.absorption_params.h2o_lut = Any[nothing,:disabled]
        @test copy_parameters(params;share_luts=true).absorption_params.h2o_lut == [nothing,:disabled]
        empty!(params.absorption_params.luts)
        @test isempty(copy_parameters(params;share_luts=true).absorption_params.luts)
        params.absorption_params = nothing
        bare = copy_parameters(params;share_luts=true)
        @test bare.absorption_params === nothing
        @test bare.p == params.p && bare.p !== params.p
    end
end

@testset "Shared LUT model construction preserves tables and Jacobians" begin
    with_logger(NullLogger()) do
        params,lut = _copy_parameter_fixture()
        coefficients = copy(lut.itp.coefs)
        outputs = map((false,true,true)) do share_luts
            trial = copy_parameters(params;share_luts)
            model,lin = model_from_parameters(LinMode(),trial;
                external_solar=true,compute_aerosol_microphysics_jacobians=false)
            @test lut.itp.coefs == coefficients
            @test params.p[end] == 1005.0
            rt_run(model,lin,1,size(lin.τ̇_abs[1],1),1;
                   jacobian_basis=:local,jacobian_adding=:source)
        end
        for result in outputs[2:end]
            @test result.toa == outputs[1].toa
            @test result.toa_jacobian == outputs[1].toa_jacobian
        end
        @test lut.itp.coefs == coefficients
    end
end
