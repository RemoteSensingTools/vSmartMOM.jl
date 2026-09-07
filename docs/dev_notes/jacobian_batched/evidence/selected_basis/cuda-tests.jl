# Run from test/ with CUDA_VISIBLE_DEVICES=0 and VSMARTMOM_SOURCE_GPU_TEST=true.
using Test, CUDA, Logging
CUDA.allowscalar(false)
source = read("test_local_jacobian.jl",String)
source = replace(source,"include(\"local_jacobian_fixture.jl\")" =>
    "include(joinpath(pwd(),\"local_jacobian_fixture.jl\"))")
Base.include_string(Main,first(split(source,"@testset \"Local optical basis CPU\"")),
                    "local_jacobian_definitions.jl")
@testset "Local optical basis CUDA" begin
    with_logger(NullLogger()) do
        check_local_jacobian(Float64,true,true;gpu=true)
        check_local_jacobian(Float32,true,false;gpu=true)
        check_local_jacobian(Float32,true,true;gpu=true,n_aerosols=3)
        check_local_jacobian(Float32,true,true;gpu=true,n_aerosols=3,
            selected_columns=[1,2,7,9,14,16,21,24],expected_basis=5)
        check_local_jacobian(Float64,true,true;gpu=true,n_aerosols=3,
            selected_columns=[1,2,7,9,10,13,14,16,19,21,24],expected_basis=8)
        # No pressure or aerosol parameter is selected: mixture directions
        # remain structurally present, but their coefficients are exact zero.
        check_local_jacobian(Float64,true,true;gpu=true,n_aerosols=3,
            selected_columns=[24],expected_basis=5)
    end
end
include(joinpath(pwd(),"test_source_adding_sif.jl"))
include(joinpath(pwd(),"test_source_adding.jl"))
