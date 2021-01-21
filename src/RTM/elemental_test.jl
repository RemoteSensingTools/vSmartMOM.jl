
# One can write the weights as a matrix (μⱼ/(μⱼ + μⱼ)); needs to be checked!
qp_matrix = qp_μ ./ (qp_μ .+ qp_μ');

qp_matrix2 = (1 ./ qp_μ .+ 1 ./ qp_μ')
# Random tau:
τ_nSpec = 4 .+ randn(1000, )

# and exp(1/μᵢ + 1/μⱼ)
mat2   = repeat(qp_matrix2, outer=[1, 1, size(τ_nSpec, 1)])
# tau_3d = repeat(τ_nSpec, inner=size(qp_matrix, 1), outer=size(qp_matrix, 1))

