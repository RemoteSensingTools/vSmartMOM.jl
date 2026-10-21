using Test, YAML, LinearAlgebra, Logging, vSmartMOM
using vSmartMOM.CoreRT
if get(ENV,"VSMARTMOM_SOURCE_GPU_TEST","false") == "true"
    using CUDA
    CUDA.allowscalar(false)
end
isdefined(@__MODULE__, :local_jacobian_fixture) || include("local_jacobian_fixture.jl")

@testset "Equivalent-source adding" begin
    with_logger(NullLogger()) do
        gpu = get(ENV,"VSMARTMOM_SOURCE_GPU_TEST","false") == "true"
        for external in (false,true), (FT,polarized) in ((Float64,true),(Float32,true),(Float64,false))
            model,lin = local_jacobian_fixture(FT,polarized,external;gpu,n_aerosols=2)
            # Embedded μ₀ is already an operator ordinate. Observe it too,
            # exercising the BOA collimated carrier and its attenuation tangent.
            external || (model.obs_geom.vza[2] = model.obs_geom.sza)
            @test all(isfinite,model.τ_aer[1])
            @test minimum(model.τ_aer[1]) >= 0
            ng = size(lin.τ̇_abs[1],1)
            matrix = rt_run(model,lin,2,ng,1;jacobian_basis=:local)
            source = rt_run(model,lin,2,ng,1;jacobian_basis=:local,jacobian_adding=:source)
            tol = FT === Float64 ? 2e-10 : 3e-4
            for (a,b) in zip(matrix,source)
                if a === nothing || b === nothing
                    @test a === b
                else
                    @test a ≈ b rtol=tol atol=10eps(FT)
                end
            end
            # Distinct azimuth tests the U component, not just a principal plane.
            polarized && @test maximum(abs,source.toa[:,3,:]) > 100eps(FT)
            @test_throws ArgumentError rt_run(model,lin,2,ng,1;
                jacobian_basis=:physical,jacobian_adding=:source)
            if FT === Float64
                z = 2; h = 1e-5
                original = copy(model.τ_abs[1][:,z])
                direction = copy(lin.τ̇_abs[1][z,:,z])
                model.τ_abs[1][:,z] .= original .+ h.*direction
                plus = local_jacobian_forward(model)
                model.τ_abs[1][:,z] .= original .- h.*direction
                minus = local_jacobian_forward(model)
                model.τ_abs[1][:,z] .= original
                for (field,jac) in ((:toa,:toa_jacobian),(:boa,:boa_jacobian))
                    getproperty(plus,field) === nothing && continue
                    fd = (getproperty(plus,field)-getproperty(minus,field))/(2h)
                    @test getproperty(source,jac)[:,:,:,15+z] ≈ fd rtol=2e-6 atol=1e-9
                end
                # The lower boundary contributes an albedo derivative even
                # at zero reflectance. Test it with a nonnegative one-sided FD.
                model.surfaces[1] = CoreRT.LambertianSurfaceScalar(0.0)
                zero_surface = rt_run(model,lin,2,ng,1;
                    jacobian_basis=:local,jacobian_adding=:source)
                base = local_jacobian_forward(model)
                model.surfaces[1] = CoreRT.LambertianSurfaceScalar(h)
                one_step = local_jacobian_forward(model)
                model.surfaces[1] = CoreRT.LambertianSurfaceScalar(2h)
                two_steps = local_jacobian_forward(model)
                for (field,jac) in ((:toa,:toa_jacobian),(:boa,:boa_jacobian))
                    getproperty(base,field) === nothing && continue
                    fd = (-3getproperty(base,field)+4getproperty(one_step,field)-getproperty(two_steps,field))/(2h)
                    @test getproperty(zero_surface,jac)[:,:,:,end] ≈ fd rtol=2e-6 atol=1e-9
                end
                for surface in (CoreRT.LambertianSurfaceLegendre([0.1,0.02,-0.01]),)
                    model.surfaces[1] = surface
                    nsurf = CoreRT.surface_parameter_count(surface)
                    dense = rt_run(model,lin,2,ng,nsurf;jacobian_basis=:local)
                    vectors = rt_run(model,lin,2,ng,nsurf;
                        jacobian_basis=:local,jacobian_adding=:source)
                    for (a,b) in zip(dense,vectors)
                        a === nothing ? (@test b === nothing) : (@test a ≈ b rtol=tol atol=1e-12)
                    end
                    forward = local_jacobian_forward(model)
                    model.τ_abs[1][:,z] .= original .+ h.*direction
                    plus = local_jacobian_forward(model)
                    model.τ_abs[1][:,z] .= original .- h.*direction
                    minus = local_jacobian_forward(model)
                    model.τ_abs[1][:,z] .= original
                    for (field,jac) in ((:toa,:toa_jacobian),(:boa,:boa_jacobian))
                        getproperty(forward,field) === nothing && continue
                        @test getproperty(vectors,field) ≈ getproperty(forward,field) rtol=tol atol=1e-12
                        fd = (getproperty(plus,field)-getproperty(minus,field))/(2h)
                        @test getproperty(vectors,jac)[:,:,:,15+z] ≈ fd rtol=2e-6 atol=1e-9
                    end
                end
            end
        end
    end
end
