import time

import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt

# Match Julia defaults
FT = jnp.float32

# Test ranges (full grid)
n_values = [16, 32, 64]
nSpec_values = [5000, 10000, 30000, 60000]

@jax.jit
def batched_mul(a, b):
    return jnp.matmul(a, b)

@jax.jit
def batched_inv(a):
    return jnp.linalg.inv(a)

num_iters = 10
key = jax.random.PRNGKey(0)

def run_case(n, nSpec):
    # Create CPU arrays, then move to device
    A = jax.random.normal(key, (n, n, nSpec), dtype=FT)
    B = jax.random.normal(key, (n, n, nSpec), dtype=FT)
    A = jax.device_put(A)
    B = jax.device_put(B)

    # Batched matmul for last two dims; transpose to (batch, n, n)
    A_b = jnp.transpose(A, (2, 0, 1))
    B_b = jnp.transpose(B, (2, 0, 1))

    # Warm-up
    batched_mul(A_b, B_b).block_until_ready()

    # Timed runs
    start = time.time()
    for _ in range(num_iters):
        batched_mul(A_b, B_b).block_until_ready()
    end = time.time()

    matmul_time = (end - start) / num_iters

    # Batched inversion (stabilize with small diagonal)
    eps = FT(1e-3)
    I = jnp.eye(n, dtype=FT)[None, :, :]
    A_b_stable = A_b + eps * I

    # Warm-up
    batched_inv(A_b_stable).block_until_ready()

    # Timed runs
    start = time.time()
    for _ in range(num_iters):
        batched_inv(A_b_stable).block_until_ready()
    end = time.time()

    inv_time = (end - start) / num_iters

    print(f"n={n}, nSpec={nSpec} | JAX batched matmul: {matmul_time:.6f} s/iter | JAX batched inv: {inv_time:.6f} s/iter")
    return matmul_time, inv_time


matmul_grid = jnp.zeros((len(n_values), len(nSpec_values)))
inv_grid = jnp.zeros((len(n_values), len(nSpec_values)))

for i, n in enumerate(n_values):
    for j, nSpec in enumerate(nSpec_values):
        matmul_time, inv_time = run_case(n, nSpec)
        matmul_grid = matmul_grid.at[i, j].set(matmul_time)
        inv_grid = inv_grid.at[i, j].set(inv_time)

# Plot heatmaps
fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)

im0 = axes[0].imshow(matmul_grid, origin="lower", aspect="auto")
axes[0].set_title("JAX batched matmul (s/iter)")
axes[0].set_xlabel("nSpec")
axes[0].set_ylabel("n")
axes[0].set_xticks(range(len(nSpec_values)), labels=[str(v) for v in nSpec_values])
axes[0].set_yticks(range(len(n_values)), labels=[str(v) for v in n_values])
fig.colorbar(im0, ax=axes[0])

im1 = axes[1].imshow(inv_grid, origin="lower", aspect="auto")
axes[1].set_title("JAX batched inv (s/iter)")
axes[1].set_xlabel("nSpec")
axes[1].set_ylabel("n")
axes[1].set_xticks(range(len(nSpec_values)), labels=[str(v) for v in nSpec_values])
axes[1].set_yticks(range(len(n_values)), labels=[str(v) for v in n_values])
fig.colorbar(im1, ax=axes[1])

plot_path = "/home/cfranken/code/gitHub/vSmartMOM.jl/test/gpu_tests/jax_batched_grid.png"
fig.savefig(plot_path, dpi=150)
print(f"Saved plot to {plot_path}")
