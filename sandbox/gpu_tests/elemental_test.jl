using KernelAbstractions
using CUDA
# One can write the weights as a matrix (μⱼ/(μⱼ + μⱼ)); needs to be checked!
# qp_matrix = qp_μ ./ (qp_μ .+ qp_μ');

# qp_matrix2 = (1 ./ qp_μ .+ 1 ./ qp_μ')
# Random tau:
# τ_nSpec = 4 .+ randn(1000, )

# and exp(1/μᵢ + 1/μⱼ)
# mat2   = repeat(qp_matrix2, outer=[1, 1, size(τ_nSpec, 1)])
# tau_3d = repeat(τ_nSpec, inner=size(qp_matrix, 1), outer=size(qp_matrix, 1))

@kernel function get_r!(r⁻⁺, t⁺⁺, ϖ, dτ, Z⁻⁺, Z⁺⁺, μ, w, ndoubl)
    i, j, n = @index(Global, NTuple)
    # @show i, j, n
    r⁻⁺[i,j,n] = ϖ[n] * Z⁻⁺[i,j] * (μ[j] / (μ[i] + μ[j])) * (1 - exp(-dτ[n] * ((1 / μ[i]) + (1 / μ[j])))) * w[j] 
    if μ[i] == μ[j]
        # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
        if i == j
            t⁺⁺[i,j,n] = exp(-dτ[n] / μ[i]) + ϖ[n] * Z⁺⁺[i,i] * (dτ[n] / μ[i]) * exp(-dτ[n] / μ[i]) * w[i]
        else
            t⁺⁺[i,j,n] = ϖ[n] * Z⁺⁺[i,i] * (dτ[n] / μ[i]) * exp(-dτ[n] / μ[i]) .* w[i]
        end
    else
        # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
        # (𝑖 ≠ 𝑗)
        t⁺⁺[i,j,n] = ϖ[n] * Z⁺⁺[i,j] .* (μ[j] / (μ[i] - μ[j])) * (exp(-dτ[n] / μ[i]) - exp.(-dτ[n] / μ[j])) * w[j]
    end
end


n = 40
nSpec = 10000
r⁻⁺_CPU, t⁺⁺_CPU = rand(n, n, nSpec), rand(n, n, nSpec);
ϖ, dτ    = rand(nSpec), rand(nSpec);
μ, w     = rand(n, n), rand(n, n);
Z⁻⁺, Z⁺⁺ = rand(n, n), rand(n, n);

device = KernelAbstractions.CPU()
kernel! = get_r!(device)
kernel!(r⁻⁺_CPU, t⁺⁺_CPU, ϖ, dτ, Z⁻⁺, Z⁺⁺, μ, w,0, ndrange=size(r⁻⁺_CPU));

function test2(r⁻⁺, t⁺⁺, ϖ, dτ, Z⁻⁺, Z⁺⁺, μ, w)
    event = kernel!(r⁻⁺, t⁺⁺, ϖ, dτ, Z⁻⁺, Z⁺⁺, μ, w, 0, ndrange=size(r⁻⁺));
    #wait(device, event)
    synchronize();
end
@time test2(r⁻⁺_CPU, t⁺⁺_CPU, ϖ, dτ, Z⁻⁺, Z⁺⁺, μ, w)
@time test2(r⁻⁺_CPU, t⁺⁺_CPU, ϖ, dτ, Z⁻⁺, Z⁺⁺, μ, w)

r⁻⁺, t⁺⁺ = CuArray(rand(n, n, nSpec)), CuArray(rand(n, n, nSpec));
ϖ, dτ    = CuArray(rand(nSpec)), CuArray(rand(nSpec));
μ, w     = CuArray(rand(n, n)), CuArray(rand(n, n));
Z⁻⁺, Z⁺⁺ = CuArray(rand(n, n)), CuArray(rand(n, n));

device = KernelAbstractions.CUDADevice()
kernel! = get_r!(device)
kernel!(r⁻⁺, t⁺⁺, ϖ, dτ, Z⁻⁺, Z⁺⁺, μ, w,0, ndrange=size(r⁻⁺));
function test2(r⁻⁺, t⁺⁺, ϖ, dτ, Z⁻⁺, Z⁺⁺, μ, w)
    event = kernel!(r⁻⁺, t⁺⁺, ϖ, dτ, Z⁻⁺, Z⁺⁺, μ, w, 0, ndrange=size(r⁻⁺));
    #wait(device, event)
    synchronize();
end
@time test2(r⁻⁺, t⁺⁺, ϖ, dτ, Z⁻⁺, Z⁺⁺, μ, w)
@time test2(r⁻⁺, t⁺⁺, ϖ, dτ, Z⁻⁺, Z⁺⁺, μ, w)
