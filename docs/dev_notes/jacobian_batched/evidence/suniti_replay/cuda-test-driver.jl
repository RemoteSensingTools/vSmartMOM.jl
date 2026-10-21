using Test, CUDA, Logging
CUDA.allowscalar(false)
source=read("test_local_jacobian.jl",String)
source=replace(source, "include(\"local_jacobian_fixture.jl\")"=>"include(joinpath(pwd(),\"local_jacobian_fixture.jl\"))")
Base.include_string(Main,first(split(source,"@testset \"Local optical basis CPU\"")),"local_jacobian_definitions.jl")
@testset "Local optical basis CUDA" begin
 with_logger(NullLogger()) do
  check_local_jacobian(Float64,true,true;gpu=true)
  check_local_jacobian(Float32,true,false;gpu=true)
  check_local_jacobian(Float32,true,true;gpu=true,n_aerosols=3)
 end
end
include(joinpath(pwd(),"test_source_adding_sif.jl"))
include(joinpath(pwd(),"test_source_adding.jl"))
