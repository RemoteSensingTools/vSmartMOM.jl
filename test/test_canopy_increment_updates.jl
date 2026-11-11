module CanopyIncrementUpdateTests

using Test, TOML, LinearAlgebra, vSmartMOM
using vSmartMOM.CoreRT

const SETTINGS = TOML.parsefile(joinpath(@__DIR__,"canopy_fixtures/increment_updates.toml"))

function check_updates(::Type{FT},N,nspec,architecture) where {FT}
    c,v = SETTINGS["cases"],SETTINGS["validation"]
    scale = FT(c["matrix_scale"])
    indices = reshape(FT.(1:N*N*nspec),N,N,nspec)
    R = scale.*sin.(indices)
    delta = scale.*sin.(indices.^2)
    C = scale.*cos.(permutedims(indices,(2,1,3)))
    source_indices = reshape(FT.(1:N*nspec),N,1,nspec)
    jd,ju = sin.(source_indices),cos.(source_indices)
    beta = FT.(c["beam_increments"][1:nspec])
    expected_R,expected_delta = similar(R),similar(delta)
    expected_down,expected_up = similar(jd),similar(ju)
    # Independent dense adding equations, deliberately forming I+C/I+Delta
    # for these moderate synthetic increments. Inputs need not commute or be
    # physically positive; this also exercises signed Stokes-like components.
    for k in 1:nspec
        r,d,a = R[:,:,k],delta[:,:,k],I+C[:,:,k]
        down,up = jd[:,1,k],ju[:,1,k]
        beam = one(FT)+beta[k]
        expected_up[:,1,k] = up+a*(r*down+beam*up)
        expected_down[:,1,k] = beam*down+a*(r*(beam*up)+down)
        expected_R[:,:,k] = r+a*r*(I+d)
        expected_delta[:,:,k] = (I+C[:,:,k])*(I+d)-I
    end
    expected_beta = (one(FT).+beta).^2 .- one(FT)
    arr = vSmartMOM.Architectures.array_type(architecture)
    r,d,a,b,down,up = arr.((R,delta,C,beta,jd,ju))
    backend = vSmartMOM.Architectures.devi(architecture)
    CoreRT.ka_canopy_increment_sources!(down,up,r,a,b,backend)
    CoreRT.ka_canopy_increment_operators!(r,d,a,b,backend)
    CoreRT.synchronize_if_gpu()
    prefix = FT===Float32 ? "float32" : "float64"
    rtol,atol = FT(v[prefix*"_rtol"]),FT(v[prefix*"_atol"])
    for (actual,expected) in zip((r,d,down,up,b),
            (expected_R,expected_delta,expected_down,expected_up,expected_beta))
        @test eltype(actual)===FT
        @test all(isfinite,Array(actual))
        @test isapprox(Array(actual),expected;rtol,atol)
    end
end

@testset "Fused canopy updates: independent dense equations" begin
    c = SETTINGS["cases"]
    for FT in (Float32,Float64),N in c["dimensions"],nspec in c["spectral_counts"]
        check_updates(FT,N,nspec,vSmartMOM.Architectures.CPU())
    end
end

if get(ENV,"VSMART_CANOPY_INCREMENT_GPU","false")=="true"
    using CUDA
    CUDA.functional() || error("CUDA required for requested fused-update test")
    CUDA.allowscalar(false)
    @testset "Fused canopy updates: CUDA dense-equation comparison" begin
        c = SETTINGS["cases"]
        for FT in (Float32,Float64),N in c["dimensions"],nspec in c["spectral_counts"]
            check_updates(FT,N,nspec,vSmartMOM.Architectures.GPU())
        end
    end
end

end
