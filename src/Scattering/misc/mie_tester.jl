# Quick testing function:
using Revise
using RadiativeTransfer.Scattering
using ForwardDiff, DiffResults
using Distributions

x = [0.3, 6.82, 1.33,0.00001]
x2 = copy(x);
δ = 1e-8
x2[3] = x[3]+δ

result = DiffResults.JacobianResult(zeros(4568),x);

result = ForwardDiff.jacobian!(result, Scattering.f_test, x);

a1= Scattering.f_test(x)


@time a2= Scattering.f_test(x2)

plot(DiffResults.jacobian(result)[:,3], label="ForwardDiff")
plot!((a2 - a1)/δ, label="Finite Diff")


# Test CUDA Matrix stuff:
elty = Float32
m = 32
A = [rand(elty,m,m) for i in 1:10000]
# move to device
d_A = CuArray{elty, 2}[]
for i in 1:length(A)
    push!(d_A,CuArray(A[i]))
end

function invCUDA(d_A)
    pivot, info = CUBLAS.getrf_batched!(d_A, true)
    pivot, info, d_C = CUBLAS.getri_batched(d_A, pivot)
    return d_C
end






