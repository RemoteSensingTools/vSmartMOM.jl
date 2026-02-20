import time
import os
os.environ["JAX_PLATFORMS"] = "cuda"

import jax
import jax.numpy as jnp

# Enable float64 support
jax.config.update("jax_enable_x64", True)

@jax.jit
def batched_mul(a, b):
    return jnp.matmul(a, b)

@jax.jit
def batched_inv(a):
    return jnp.linalg.inv(a)

matrix_sizes = [12, 24, 36]
batch_sizes  = [6837, 10000, 50000]
num_iters    = 20
key = jax.random.PRNGKey(0)

print(f"JAX batched operations benchmark")
print(f"JAX version: {jax.__version__}")
print(f"JAX backend: {jax.default_backend()}")
devs = jax.devices()
for d in devs:
    print(f"Device: {d}")
print("=" * 75)

for FT, label in [(jnp.float64, "Float64"), (jnp.float32, "Float32")]:
    print(f"\n── {label} ──")
    print()
    print(f"  {'size':<8s}  {'batch':>8s}  │ {'matmul (ms)':>12s}  {'μs/mat':>10s}  │ {'inv (ms)':>12s}  {'μs/mat':>10s}")
    print("  " + "─" * 71)

    for n in matrix_sizes:
        for nBatch in batch_sizes:
            A = jax.random.normal(key, (nBatch, n, n), dtype=FT)
            B = jax.random.normal(key, (nBatch, n, n), dtype=FT)

            # Warmup matmul
            C = batched_mul(A, B)
            C.block_until_ready()

            mul_times = []
            for _ in range(num_iters):
                start = time.perf_counter()
                C = batched_mul(A, B)
                C.block_until_ready()
                mul_times.append(time.perf_counter() - start)

            mul_med = sorted(mul_times)[num_iters // 2] * 1000

            # Stabilize for inversion
            eye = jnp.eye(n, dtype=FT)[None, :, :]
            A_stable = A + float(n) * eye
            A_stable.block_until_ready()

            # Warmup inv
            X = batched_inv(A_stable)
            X.block_until_ready()

            inv_times = []
            for _ in range(num_iters):
                start = time.perf_counter()
                X = batched_inv(A_stable)
                X.block_until_ready()
                inv_times.append(time.perf_counter() - start)

            inv_med = sorted(inv_times)[num_iters // 2] * 1000

            print(f"  {n:2d}×{n:<4d}  {nBatch:8d}  │ {mul_med:10.3f} ms  {mul_med*1000/nBatch:8.3f} μs  │ {inv_med:10.3f} ms  {inv_med*1000/nBatch:8.3f} μs")

            del A, B, A_stable
