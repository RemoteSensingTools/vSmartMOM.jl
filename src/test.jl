using KernelAbstractions
using CUDA
using Plots

@kernel function line_shape!(A, @Const(grid),γ_d,ν  )
    cSqrtLn2divSqrtPi  = Float32(0.469718639319144059835)
    cLn2               = Float32(0.6931471805599)
    #γ_d = Float32(0.08);
    #ν = Float32(6050);
    I = @index(Global)
    @inbounds A[I] = cSqrtLn2divSqrtPi*exp(-cLn2*((grid[I] - ν) /γ_d) ^2) /γ_d
end

@kernel function line_shape!(A, @Const(grid))
    I = @index(Global)
    @inbounds A[I] = exp(-Float32(0.6931471805599)*((grid[I]-Float32(6050))/Float32(0.1))^2)/Float32(0.1)
end

kernelCPU = line_shape!(CPU(), 1)
@time kernelCPU(out, grid, ndrange=length(grid))


grid = collect(Float32,6000:0.001:6300);
out = similar(grid)
const γ_d = Float32(0.08);
const ν = Float32(6050);
const S = Float32(1e-20)
gridC = CuArray(grid);
outC = CuArray(out);
const y = Float32(1.0)
kernelGPU = line_shape!(CUDADevice(), 256)

@time event = tester!(outC, gridC); wait(event)
plot(grid, out)

function testerVoigt!(A::CuArray, B::CuArray)
    @assert size(A) == size(B)
    voigt_shape32!(CUDADevice(), 256)(A, B,γ_d,ν,y,S, ndrange=length(A))
end

function testerVoigt!(A::Array, B::Array)
    @assert size(A) == size(B)
    voigt_shape32!(CPU(), 4)(A, B,γ_d,ν,y,S, ndrange=length(A))
end

function tester!(A::CuArray, B::CuArray)
    @assert size(A) == size(B)
    line_shape!(CUDADevice(), 256)(A, B,γ_d,ν, ndrange=length(A))
end

function tester!(A::Array, B::Array)
    @assert size(A) == size(B)
    kernel = line_shape!(CPU(), 8)
    kernel(A, B, ndrange=length(A))
end