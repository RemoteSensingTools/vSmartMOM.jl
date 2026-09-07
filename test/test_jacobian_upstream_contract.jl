using Test, vSmartMOM
using vSmartMOM.CoreRT

# Exercise extension-defined plans, independent of OCO's parameter names.
struct _ContractFlavor <: AbstractJacobianFlavor end
function _contract_plan(columns)
    keys = [ParameterKey(:custom, Symbol("parameter_$i")) for i in eachindex(columns)]
    names = String[string(k.field) for k in keys]
    layout = ActiveParameterLayout(keys, names, keys, names, columns, length(columns))
    JacobianPlan(_ContractFlavor(), keys, names, [layout])
end
CoreRT.requires_aerosol_microphysics_jacobians(::_ContractFlavor) = false
CoreRT.jacobian_plan(::_ContractFlavor, params, model, lin) = _contract_plan([3])
@testset "Selected derivatives require upstream availability" begin
    validate(columns; micro=false, h2o=false) = CoreRT._validate_plan_upstream(
        _contract_plan(columns), 1, 2, [4];
        compute_aerosol_microphysics_jacobians=micro,
        compute_h2o_jacobians=h2o)
    # Native: pressure 1; aerosol 2:8; q-H2O 9:10; variable gas 11:12.
    @test validate([1,2,7,8,11,12]) === nothing
    @test validate(Int[]) === nothing
    for micro_column in 3:6
        @test_throws ArgumentError validate([micro_column])
        @test validate([micro_column]; micro=true) === nothing
    end
    for water_column in 9:10
        @test_throws ArgumentError validate([water_column])
        @test validate([water_column]; h2o=true) === nothing
    end
    @test validate(collect(1:12); micro=true,h2o=true) === nothing
    @test_throws ArgumentError validate([13]; micro=true,h2o=true)
    @test_throws DimensionMismatch CoreRT._validate_plan_upstream(
        _contract_plan([1]), 1, 2, [4,4];
        compute_aerosol_microphysics_jacobians=true,compute_h2o_jacobians=true)
end

struct _UnimplementedJacobianSource <: AbstractSource end
@testset "Forward source declarations do not imply tangent support" begin
    @test source_ad_mode(_UnimplementedJacobianSource()) isa AnalyticSourceJacobian
    @test !CoreRT._linearized_source_supported(_UnimplementedJacobianSource())
    @test !CoreRT._linearized_source_supported(ThermalEmission())
    @test CoreRT._linearized_source_supported(SolarBeam() + SurfaceSIF())
    params = parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
    params.architecture = vSmartMOM.CPU()
    @test_throws ArgumentError model_from_parameters(_ContractFlavor(), params;
        compute_h2o_jacobians=false)
    model, lin = model_from_parameters(LinMode(), params;
        compute_aerosol_microphysics_jacobians=false,compute_h2o_jacobians=false)
    na = CoreRT.n_aerosols(model)
    ng = size(lin.τ̇_abs[1],1)
    ns = CoreRT.surface_parameter_count(CoreRT.get_surface(model,1))
    for source in (ThermalEmission(),SolarBeam()+ThermalEmission(),
                   _UnimplementedJacobianSource()), adding in (:matrix,:source)
        @test_throws ArgumentError rt_run_lin(model,lin,na,ng,ns;
            sources=source,jacobian_basis=:local,jacobian_adding=adding)
    end
    thermal_model = model_from_parameters(params; sources=ThermalEmission())
    @test_throws ArgumentError rt_run_lin(thermal_model,lin,na,ng,ns)
end
