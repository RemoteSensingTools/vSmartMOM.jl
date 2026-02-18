using Revise
using vSmartMOM
using CUDA
using Test
using BenchmarkTools
using LinearAlgebra
using NNlib


# Needs some warning if memory is getting too large !
FT = Float32

# Sweep settings (full grid)
n_values = [16, 32, 64]
nSpec_values = [5000, 10000, 30000, 60000]
num_iters = 5

function run_inv_case(n, nSpec)
	# Create CPU Matrices:
	A_ = (randn(FT, n, n, nSpec));

	# And move to GPU as CuArray
	A = CuArray(A_);

	# Batched matrix inversion tests
	A_inv = similar(A);
	A_for_inv = copy(A);
	A_inv_ = similar(A_);
	A_for_inv_ = copy(A_);

	# Pointer-based batch inversion (GPU)
	A_inv_ptrs = CUDA.CUBLAS.unsafe_strided_batch(A_inv);
	A_for_inv_ptrs = CUDA.CUBLAS.unsafe_strided_batch(A_for_inv);

	# Warm-up
	CUDA.@sync vSmartMOM.CoreRT.batch_inv!(A_inv, A_for_inv, A_inv_ptrs, A_for_inv_ptrs)

	# Timed runs
	t_start = time()
	for _ in 1:num_iters
		CUDA.@sync vSmartMOM.CoreRT.batch_inv!(A_inv, A_for_inv, A_inv_ptrs, A_for_inv_ptrs)
	end
	t_elapsed = (time() - t_start) / num_iters

	println("n=$(n), nSpec=$(nSpec) | Julia batched inv: $(round(t_elapsed, digits=6)) s/iter")
	return t_elapsed

end

times = zeros(length(n_values), length(nSpec_values))

for (i, n) in pairs(n_values)
	for (j, nSpec) in pairs(nSpec_values)
		times[i, j] = run_inv_case(n, nSpec)
	end
end

# Heatmap plot
heatmap = UnicodePlots.heatmap(
	nSpec_values,
	n_values,
	times;
	xlabel="nSpec",
	ylabel="n",
	title="Julia batched inv (s/iter)",

	colormap=:viridis,
)

println(heatmap)