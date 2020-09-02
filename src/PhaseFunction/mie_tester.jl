# Quick testing function:
using Revise
using RadiativeTransfer.PhaseFunction
using ForwardDiff, DiffResults

x = [0.3, 6.82, 1.33,0.00001]
x2 = copy(x);
δ = 1e-8
x2[3] = x[3]+δ

result = DiffResults.JacobianResult(zeros(4568),x);

result = ForwardDiff.jacobian!(result, PhaseFunction.f_test, x);

a1= PhaseFunction.f_test(x)


@time a2= PhaseFunction.f_test(x2)

plot(DiffResults.jacobian(result)[:,3], label="ForwardDiff")
plot!((a2 - a1)/δ, label="Finite Diff")










