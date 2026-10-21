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
    validate(columns; micro=false, h2o=false, pressure=true) = CoreRT._validate_plan_upstream(
        _contract_plan(columns), 1, 2, [4];
        compute_aerosol_microphysics_jacobians=micro,
        compute_h2o_jacobians=h2o, compute_pressure_jacobians=pressure)
    # Native: pressure 1; aerosol 2:8; q-H2O 9:10; variable gas 11:12.
    @test validate([1,2,7,8,11,12]) === nothing
    @test validate(Int[]) === nothing
    @test validate([2,7,8,11,12]; pressure=false) === nothing
    @test_throws ArgumentError validate([1]; pressure=false)
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

# A plan with surface pressure and the surface coefficient only.
struct _PressureSurfaceFlavor <: AbstractJacobianFlavor end
CoreRT.requires_aerosol_microphysics_jacobians(::_PressureSurfaceFlavor) = false
function CoreRT.jacobian_plan(::_PressureSurfaceFlavor, params, model, lin)
    keys = [ParameterKey(:atmosphere, :surface_pressure),
            ParameterKey(:surface, :albedo; component=1, band=1)]
    names = ["surface_pressure", "albedo"]
    layout = ActiveParameterLayout(keys, names, keys, names, [1], 1;
                                   surface_columns=2:2)
    JacobianPlan(_PressureSurfaceFlavor(), keys, names, [layout])
end

@testset "Raman, thermal and TMS tangents are rejected at every linearized entry" begin
    params = parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
    params.architecture = vSmartMOM.CPU()
    model, lin = model_from_parameters(LinMode(), params;
        compute_aerosol_microphysics_jacobians=false)
    na, ng = CoreRT.n_aerosols(model), size(lin.τ̇_abs[1], 1)
    ns = CoreRT.surface_parameter_count(CoreRT.get_surface(model, 1))

    # Inelastic (RRS and both VS directions) on the explicit Raman-type entry.
    FT = Float64
    IS = vSmartMOM.InelasticScattering
    ν = model.atmosphere.spec_bands[1]
    n2, o2 = IS.getRamanAtmoConstants(ν[1], FT(250))
    greek = IS.GreekCoefs(ones(FT, 1), ones(FT, 1), ones(FT, 1),
                                           ones(FT, 1), ones(FT, 1), ones(FT, 1))
    shared = (n2=n2, o2=o2, greek_raman=greek, ϖ_λ₁λ₀=zeros(FT, 1),
              i_λ₁λ₀=zeros(Int, 1), Z⁻⁺_λ₁λ₀=zeros(FT, 1, 1),
              Z⁺⁺_λ₁λ₀=zeros(FT, 1, 1), n_Raman=0,
              F₀=zeros(FT, 1, 1), SIF₀=zeros(FT, 1, 1))
    raman_types = (
        IS.RRS(; shared..., fscattRayl=[one(FT)],
            ϖ_Cabannes=[one(FT)], i_ref=1),
        IS.VS_0to1(; shared..., fscattRayl=one(FT),
            ϖ_Cabannes=one(FT), dτ₀=zero(FT), dτ₀_λ=zero(FT), k_Rayl_scatt=one(FT)),
        IS.VS_1to0(; shared..., fscattRayl=one(FT),
            ϖ_Cabannes=one(FT), dτ₀=zero(FT), dτ₀_λ=zero(FT), k_Rayl_scatt=one(FT)))
    for rs in raman_types
        @test_throws r"Linearized Raman-active RT is intentionally unsupported" rt_run(
            rs, model, lin, na, ng, ns, 1)
    end

    # Retrieval-selected solves share the source and TMS guards.
    planned_model, planned = model_from_parameters(_PressureSurfaceFlavor(), params)
    @test all(isfinite, rt_run(planned_model, planned).toa_jacobian)
    for source in (ThermalEmission(), SolarBeam() + ThermalEmission()),
            adding in (:matrix, :source)
        @test_throws r"requires end-to-end source tangent propagation" rt_run(
            planned_model, planned; sources=source, jacobian_basis=:local,
            jacobian_adding=adding)
    end
    thermal_model, thermal_planned = model_from_parameters(_PressureSurfaceFlavor(), params;
        sources=ThermalEmission())
    @test_throws ArgumentError rt_run(thermal_model, thermal_planned)

    tms = deepcopy(params)
    tms.numerics = CoreRT.RTNumericalParameters{FT}(ss_correction=CoreRT.TMSCorrection())
    tms_model, tms_lin = model_from_parameters(LinMode(), tms;
        compute_aerosol_microphysics_jacobians=false)
    tms_planned_model, tms_planned = model_from_parameters(_PressureSurfaceFlavor(), tms)
    for basis in (:physical, :local)
        @test_throws r"TMSCorrection is not implemented in the linearized driver" rt_run(
            tms_model, tms_lin, na, ng, ns; jacobian_basis=basis)
        @test_throws r"TMSCorrection is not implemented in the linearized driver" rt_run(
            tms_planned_model, tms_planned; jacobian_basis=basis)
    end
end
