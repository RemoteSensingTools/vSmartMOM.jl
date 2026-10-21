using Test, YAML, LinearAlgebra, Logging, vSmartMOM
using vSmartMOM.CoreRT
if get(ENV,"VSMARTMOM_SOURCE_GPU_TEST","false") == "true"
    using CUDA
    CUDA.allowscalar(false)
end
isdefined(@__MODULE__, :local_jacobian_fixture) || include("local_jacobian_fixture.jl")

"Forward reference with the same emission and endpoint contract as the tangent solve."
source_sif_forward(model, sources) = model.quad_points.external_solar ?
    (; toa=rt_run_toa(model; sources), boa=nothing) : rt_run(model; sources)

@testset "Equivalent-source adding with SIF" begin
    with_logger(NullLogger()) do
        gpu = get(ENV,"VSMARTMOM_SOURCE_GPU_TEST","false") == "true"
        for external in (false,true), (FT,polarized) in ((Float64,true),(Float32,true),(Float64,false))
            model,lin = local_jacobian_fixture(FT,polarized,external;gpu)
            ν = CoreRT.get_spec_bands(model)[1]
            # An explicit non-unit solar spectrum keeps physical SIF units.
            F₀ = zeros(FT,polarized ? 3 : 1,length(ν))
            F₀[1,:] .= FT.(range(1.5,2.5;length=length(ν)))
            solar = SolarBeam(; F₀)
            # Center this synthetic spectrum in its test band. The reference
            # coordinate is configurable independently of the O2-band default.
            ν_ref = sum(ν)/length(ν)
            emission(a,b) = SurfaceSIF(SIF760=a,mSIF=b,wavenumber_cm1=ν;ν_ref)
            ng = size(lin.τ̇_abs[1],1)
            tol = FT === Float64 ? 2e-10 : 3e-4
            for brdf in (CoreRT.LambertianSurfaceScalar(FT(0.05)),
                         CoreRT.LambertianSurfaceLegendre(FT[0.1,0.02,-0.01]))
                model.surfaces[1] = brdf
                nsurf = CoreRT.surface_parameter_count(brdf)
                layout = CoreRT.ParameterLayout(n_aerosols=1,n_gases=ng,n_surface=nsurf,n_sif=2)
                for amplitude in (0.0,0.2)
                    slope = iszero(amplitude) ? 0.0 : 1e-4
                    sources = solar + emission(amplitude,slope)
                    matrix = rt_run(model,lin,1,ng,nsurf;sources,jacobian_basis=:local)
                    source = rt_run(model,lin,1,ng,nsurf;sources,
                        jacobian_basis=:local,jacobian_adding=:source)
                    forward = source_sif_forward(model,sources)
                    for (a,b) in zip(matrix,source)
                        a === nothing ? (@test b === nothing) :
                            (@test a ≈ b rtol=tol atol=10eps(FT))
                    end
                    for (field,jac) in ((:toa,:toa_jacobian),(:boa,:boa_jacobian))
                        getproperty(forward,field) === nothing && continue
                        @test getproperty(source,field) ≈ getproperty(forward,field) rtol=tol atol=10eps(FT)
                        # The amplitude derivative remains nonzero at zero SIF.
                        @test maximum(abs,getproperty(source,jac)[:,:,:,CoreRT.sif755_index(layout)]) > 0.01
                    end
                    FT === Float64 || continue
                    h = 1e-5
                    for (k,p) in enumerate(CoreRT.sif_range(layout))
                        plus = source_sif_forward(model,solar + emission(amplitude+(k==1 ? h : 0),slope+(k==2 ? h : 0)))
                        minus = source_sif_forward(model,solar + emission(amplitude-(k==1 ? h : 0),slope-(k==2 ? h : 0)))
                        for (field,jac) in ((:toa,:toa_jacobian),(:boa,:boa_jacobian))
                            getproperty(plus,field) === nothing && continue
                            fd = (getproperty(plus,field)-getproperty(minus,field))/(2h)
                            @test getproperty(source,jac)[:,:,:,p] ≈ fd rtol=2e-6 atol=1e-9
                        end
                    end
                    # Absorption perturbs transport of SIF, but never its
                    # emitted boundary amplitude. This catches an erroneous
                    # -j_SIF dτ_column/μ₀ in surface-source expansion.
                    z = 2
                    original = copy(model.τ_abs[1][:,z])
                    direction = copy(lin.τ̇_abs[1][z,:,z])
                    model.τ_abs[1][:,z] .= original .+ h.*direction
                    plus = source_sif_forward(model,sources)
                    model.τ_abs[1][:,z] .= original .- h.*direction
                    minus = source_sif_forward(model,sources)
                    model.τ_abs[1][:,z] .= original
                    p = CoreRT.gas_range(layout)[z]
                    for (field,jac) in ((:toa,:toa_jacobian),(:boa,:boa_jacobian))
                        getproperty(plus,field) === nothing && continue
                        fd = (getproperty(plus,field)-getproperty(minus,field))/(2h)
                        @test getproperty(source,jac)[:,:,:,p] ≈ fd rtol=2e-6 atol=1e-9
                    end
                    # Surface reflectance also changes the return of emitted
                    # photons after they scatter back down from the atmosphere.
                    perturb_surface(δ) = brdf isa CoreRT.LambertianSurfaceScalar ?
                        CoreRT.LambertianSurfaceScalar(FT(0.05)+δ) :
                        CoreRT.LambertianSurfaceLegendre(FT[0.1+δ,0.02,-0.01])
                    model.surfaces[1] = perturb_surface(h)
                    plus = source_sif_forward(model,sources)
                    model.surfaces[1] = perturb_surface(-h)
                    minus = source_sif_forward(model,sources)
                    model.surfaces[1] = brdf
                    for (field,jac) in ((:toa,:toa_jacobian),(:boa,:boa_jacobian))
                        getproperty(plus,field) === nothing && continue
                        fd = (getproperty(plus,field)-getproperty(minus,field))/(2h)
                        @test getproperty(source,jac)[:,:,:,first(CoreRT.surface_range(layout))] ≈ fd rtol=2e-6 atol=1e-9
                    end
                end
                # Prescribed emission has no SIF columns, but still changes
                # the atmospheric and albedo derivatives through illumination.
                prescribed = solar + SurfaceSIF(SIF₀=FT(0.1).*F₀)
                matrix = rt_run(model,lin,1,ng,nsurf;sources=prescribed,jacobian_basis=:local)
                source = rt_run(model,lin,1,ng,nsurf;sources=prescribed,
                    jacobian_basis=:local,jacobian_adding=:source)
                @test source.toa_jacobian ≈ matrix.toa_jacobian rtol=tol atol=10eps(FT)
                if FT === Float64 && polarized && external && nsurf == 1
                    # Multiple source terms retain separate parameter slots;
                    # the legacy wavelength coordinate shares the same dispatch.
                    legacy = SurfaceSIF(SIF755=0.1,slope=1e-4,wavelength_nm=1e7 ./ ν)
                    sources = legacy + solar + emission(0.2,1e-4)
                    matrix = rt_run(model,lin,1,ng,nsurf;sources,jacobian_basis=:local)
                    source = rt_run(model,lin,1,ng,nsurf;sources,
                        jacobian_basis=:local,jacobian_adding=:source)
                    @test size(source.toa_jacobian,4) == CoreRT.n_total(layout)+2
                    @test source.toa_jacobian ≈ matrix.toa_jacobian rtol=tol atol=10eps(FT)
                end
            end
        end
    end
end
