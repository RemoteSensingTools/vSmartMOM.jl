module CanopyIncrementTests

using Test, TOML, LinearAlgebra, vSmartMOM, CanopyOptics
using vSmartMOM.CoreRT, vSmartMOM.InelasticScattering

const SETTINGS = TOML.parsefile(joinpath(@__DIR__,"canopy_fixtures/increment.toml"))
const OPERATORS = (:r⁻⁺,:r⁺⁻,:t⁺⁺,:t⁻⁻,:j₀⁺,:j₀⁻)
const ORIGINAL_BLAS_THREADS = BLAS.get_num_threads()
BLAS.set_num_threads(SETTINGS["numerics"]["blas_threads"])

function context(::Type{FT}, npol, nspec, specular, architecture) where {FT}
    g,o=SETTINGS["geometry"],SETTINGS["optics"]
    pol=npol==1 ? CoreRT.Stokes_I{FT}() : npol==3 ? CoreRT.Stokes_IQU{FT}() : CoreRT.Stokes_IQUV{FT}()
    observation=CoreRT.ObsGeometry(FT(g["solar_zenith_degrees"]),FT[0],FT[0],FT[0])
    arr=vSmartMOM.Architectures.array_type(architecture)
    qp=CoreRT.rt_set_streams(CoreRT.GaussLegQuad(),2*g["streams"]-1,observation,pol,arr)
    dim=length(qp.qp_μ)*npol
    added=CoreRT.make_added_layer(noRS{FT}(),FT,arr,(dim,dim),nspec)
    diffuse=BiLambertianCanopyScattering(R=FT(o["leaf_reflectance"]),T=FT(o["leaf_transmittance"]))
    scattering=specular ? diffuse+SpecularCanopyScattering(
        nᵣ=FT(o["specular_refractive_index"]),κ=FT(o["specular_extinction_index"])) : diffuse
    canopy=CoreRT.CanopySurface(;soil=CoreRT.LambertianSurfaceScalar(FT(o["soil_albedo"])),
        LAI=FT(first(g["lai"])),LAD=beta_leaves(FT(g["leaf_beta_alpha"]),FT(g["leaf_beta_beta"])),
        canopy_scattering=scattering,leaf_reflectance=FT(o["leaf_reflectance"]),
        leaf_transmittance=FT(o["leaf_transmittance"]),
        canopy_quadrature=CanopyQuadrature(n_leaf=g["leaf_nodes"],n_azimuth=g["azimuth_nodes"]))
    CoreRT._init_canopy_cache!(canopy,added,pol,qp,architecture;m_max=maximum(g["fourier_orders"]))
    F0=FT.(o["solar_stokes"][1:npol])*transpose(FT.(o["solar_spectral_scale"][1:nspec]))
    (;pol,qp,dim,nspec,arr,architecture,cache=canopy._cache,F0=arr(F0))
end

function seed(ctx,::Type{FT},lai,m,floor;SFI=true) where {FT}
    core=CoreRT.CoreDirectionalScatteringOpticalProperties(
        ctx.arr(fill(FT(lai),ctx.nspec)),
        ctx.arr(FT.(SETTINGS["optics"]["single_scattering_albedos"][1:ctx.nspec])),
        ctx.arr(ctx.cache.Z_coarse_pp[:,:,1,m+1]),ctx.arr(ctx.cache.Z_coarse_mp[:,:,1,m+1]),ctx.cache.G)
    dt,nd,beam=CoreRT.init_layer(core,ctx.qp,ctx.pol,ctx.architecture;
        dτ_min_floor=FT(floor),dτ_max_threshold=FT(SETTINGS["numerics"]["thin_threshold"]))
    added=CoreRT.make_added_layer(noRS{FT}(),FT,ctx.arr,(ctx.dim,ctx.dim),ctx.nspec)
    CoreRT.elemental!(ctx.pol,SFI,ctx.arr(zeros(FT,ctx.nspec)),dt,ctx.F0,
        core,m,nd,true,ctx.qp,added,ctx.architecture)
    (;added,core,dt,nd,beam)
end

function result(ctx,::Type{FT},lai,m,floor;increment=true,SFI=true) where {FT}
    s=seed(ctx,FT,lai,m,floor;SFI)
    if increment
        CoreRT.canopy_increment_doubling!(ctx.pol,SFI,s.dt,s.core,m,s.nd,
            ctx.qp,s.added,ctx.cache.I_static,ctx.architecture)
    else
        CoreRT.doubling!(ctx.pol,SFI,s.beam,s.nd,s.added,ctx.cache.I_static,ctx.architecture)
    end
    NamedTuple{OPERATORS}(Tuple(Array(getproperty(s.added,key)) for key in OPERATORS))
end

function compare(actual,reference,::Type{FT};sources=true) where {FT}
    v=SETTINGS["validation"]
    prefix=FT===Float32 ? "float32" : "float64"
    for key in (sources ? OPERATORS : OPERATORS[1:4])
        x,y=getproperty(actual,key),getproperty(reference,key)
        @test all(isfinite,x)
        @test isapprox(x,y;rtol=v[prefix*"_rtol"],atol=v[prefix*"_atol"])
    end
end

@testset "Canopy transmission increments: CPU" begin
    g,n=SETTINGS["geometry"],SETTINGS["numerics"]
    cpu=vSmartMOM.Architectures.CPU()
    for npol in g["stokes_components"], nspec in g["spectral_counts"], specular in (false,true)
        c32=context(Float32,npol,nspec,specular,cpu)
        c64=context(Float64,npol,nspec,specular,cpu)
        for lai in g["lai"], m in g["fourier_orders"]
            fine32=result(c32,Float32,lai,m,n["fine_floor"])
            fine64=result(c64,Float64,lai,m,n["fine_floor"])
            compare(fine32,fine64,Float32)
            @test eltype(fine32.t⁺⁺)===Float32
            original=result(c64,Float64,lai,m,n["algebra_floor"];increment=false)
            rewritten=result(c64,Float64,lai,m,n["algebra_floor"])
            compare(rewritten,original,Float64)
        end
        # Diffuse operators must not depend on whether the solar source is used.
        dark=result(c64,Float64,last(g["lai"]),first(g["fourier_orders"]),n["fine_floor"];SFI=false)
        lit=result(c64,Float64,last(g["lai"]),first(g["fourier_orders"]),n["fine_floor"])
        compare(dark,lit,Float64;sources=false)
    end
end

if get(ENV,"VSMART_CANOPY_INCREMENT_GPU","false")=="true"
    using CUDA
    CUDA.functional() || error("CUDA requested for canopy increment validation but is unavailable")
    CUDA.allowscalar(false)
    @testset "Canopy transmission increments: CUDA" begin
        g,n=SETTINGS["geometry"],SETTINGS["numerics"]
        for FT in (Float32,Float64), npol in g["stokes_components"],
                nspec in g["spectral_counts"], specular in (false,true)
            cpu=context(FT,npol,nspec,specular,vSmartMOM.Architectures.CPU())
            gpu=context(FT,npol,nspec,specular,vSmartMOM.Architectures.GPU())
            for lai in g["lai"],m in g["fourier_orders"]
                reference=result(cpu,FT,lai,m,n["fine_floor"])
                device=result(gpu,FT,lai,m,n["fine_floor"])
                compare(device,reference,FT)
            end
        end
        println("Canopy increment kernels verified on ",CUDA.name(CUDA.device()))
    end
end

BLAS.set_num_threads(ORIGINAL_BLAS_THREADS)
end
