#=
KernelAbstractions-based Mie scattering kernels.

These kernels are backend-agnostic: they run on CPU() or CUDABackend() via
KernelAbstractions.jl.  The precision tier for the Dn recursion is selected
by the MiePrecisionPolicy parameter.

Kernel 1 -- Mie coefficients (an, bn) per radius point
Kernel 2 -- Scattering amplitudes S1, S2 per (angle, radius) point
Kernel 3 -- Phase matrix elements f11, f33, f12, f34  (fused with Kernel 2)
Kernel 4 -- Size-distribution reduction (bulk cross-sections + bulk phase matrix)
Kernel 5 -- Greek coefficients via Legendre integration
=#

using KernelAbstractions

# ============================================================================
# Kernel 1: Mie coefficients -- one thread per radius point
# ============================================================================

"""
    mie_coefficients_kernel!(an, bn, x_params, m_re, m_im, nmax_per_r, nmx_max)

Compute Mie coefficients an[i, 1:nmax_i] and bn[i, 1:nmax_i] for each radius
index i.  The Dn downward recursion uses DoubleSingle emulated FP64 (Tier 1).
Riccati-Bessel forward recursion and an/bn computation use native FT (Tier 3).

Layout: an, bn are (nquad_radius, nmax_global) -- transposed for coalescing.
"""
@kernel function mie_coefficients_kernel_ds!(
    an, bn,                     # output: (nquad_radius, nmax_global) Complex{FT}
    @Const(x_params),           # input:  (nquad_radius,) FT -- size parameters
    m_re, m_im,                 # scalars: real/imag refractive index parts
    @Const(nmax_per_r),         # input:  (nquad_radius,) Int -- per-radius nmax
    # FOLLOW-UP (dead parameter, pre-existing): despite the name, the body
    # below never reads `nmx_max` -- it recomputes its own per-thread
    # `nmx_i` from `n_max_i`/`y_re`/`y_im` instead. This predates the batched
    # seam (compute_NAI2_nodes_batched.jl) and is now propagated, unchanged,
    # into its `_batched` sibling kernel (mie_coefficients_kernel_ds_batched!)
    # for signature parity. Candidate cleanup: drop this parameter from all
    # four kernel variants (this one, `_f64!`, and both `_batched!` siblings)
    # in one pass, verifying bit-compat is preserved.
    nmx_max::Int                # scalar: max Dn recursion depth
)
    i = @index(Global)
    FT = eltype(x_params)

    x = x_params[i]
    n_max_i = nmax_per_r[i]

    # Complex refractive index and y = x * m
    # Using DS for the recursion
    y_re = x * m_re
    y_im = x * m_im
    nmx_i = max(n_max_i, round(Int, sqrt(y_re^2 + y_im^2)) + 1) + 50

    # --- Dn downward recursion in DoubleSingle precision ---
    # Dn[nmx] = 0, then Dn[n] = (n+1)/y - 1/(Dn[n+1] + (n+1)/y)
    Dn_re = DoubleSingle{FT}(zero(FT))
    Dn_im = DoubleSingle{FT}(zero(FT))
    Dn_prev = ComplexDS(Dn_re, Dn_im)

    # We store Dn values that we need (1:n_max_i) in registers/local memory
    # Since we're going downward, we accumulate and only store when n <= n_max_i
    # We'll use the output arrays to store intermediate Dn values
    # First pass: compute all Dn[1:n_max_i] via downward recursion

    # Convert y to ComplexDS
    y_ds = cds_complex(FT(y_re), FT(y_im))

    @inbounds for n = (nmx_i - 1):-1:1
        # ratio = (n+1) / y
        n_plus_1 = DoubleSingle{FT}(FT(n + 1))
        ratio = cds_div(cds_from_real(n_plus_1), y_ds)

        # Dn_prev = ratio - 1 / (Dn_prev + ratio)
        sum_val = cds_add(Dn_prev, ratio)
        inv_sum = cds_inv(sum_val)
        Dn_prev = cds_sub(ratio, inv_sum)

        # Store if within needed range
        if n <= n_max_i
            # Convert DS to native FT for storage
            dn_native = to_complex(Dn_prev)
            # Store in a temporary location -- we'll use an array for Dn
            # Since we need Dn[n] for the forward pass, store in bn temporarily
            # Actually, store directly -- we need an intermediate buffer
            # Use the bn array as temporary Dn storage (will be overwritten)
            bn[i, n] = dn_native
        end
    end

    # --- Riccati-Bessel forward recursion (native FT, Tier 3) ---
    ψ₀ = cos(x)
    ψ₁ = sin(x)
    χ₀ = -sin(x)
    χ₁ = cos(x)
    ξ₁_re = ψ₁
    ξ₁_im = χ₁

    @inbounds for n = 1:n_max_i
        # Forward recursion for ψ and χ
        coef = FT(2n - 1) / x
        ψ = coef * ψ₁ - ψ₀
        χ = coef * χ₁ - χ₀
        ξ_re = ψ
        ξ_im = χ

        # Retrieve Dn[n] from temporary storage
        Dn_n = bn[i, n]  # stored during downward recursion
        Dn_re_n = real(Dn_n)
        Dn_im_n = imag(Dn_n)

        # t_a = Dn[n] / m + n/x  (complex / complex + real)
        # Dn/m: (Dn_re + i*Dn_im) / (m_re + i*m_im)
        inv_m_denom = one(FT) / (m_re^2 + m_im^2)
        ta_re = (Dn_re_n * m_re + Dn_im_n * m_im) * inv_m_denom + FT(n) / x
        ta_im = (Dn_im_n * m_re - Dn_re_n * m_im) * inv_m_denom

        # t_b = Dn[n] * m + n/x
        tb_re = Dn_re_n * m_re - Dn_im_n * m_im + FT(n) / x
        tb_im = Dn_re_n * m_im + Dn_im_n * m_re

        # an[n] = (t_a * ψ - ψ₁) / (t_a * ξ - ξ₁)
        # Numerator: (ta_re + i*ta_im) * ψ - ψ₁
        num_a_re = ta_re * ψ - ψ₁
        num_a_im = ta_im * ψ

        # Denominator: (ta_re + i*ta_im) * (ξ_re + i*ξ_im) - (ξ₁_re + i*ξ₁_im)
        den_a_re = ta_re * ξ_re - ta_im * ξ_im - ξ₁_re
        den_a_im = ta_re * ξ_im + ta_im * ξ_re - ξ₁_im

        # Complex division
        inv_den_a = one(FT) / (den_a_re^2 + den_a_im^2)
        an[i, n] = Complex{FT}(
            (num_a_re * den_a_re + num_a_im * den_a_im) * inv_den_a,
            (num_a_im * den_a_re - num_a_re * den_a_im) * inv_den_a
        )

        # bn[n] = (t_b * ψ - ψ₁) / (t_b * ξ - ξ₁)
        num_b_re = tb_re * ψ - ψ₁
        num_b_im = tb_im * ψ

        den_b_re = tb_re * ξ_re - tb_im * ξ_im - ξ₁_re
        den_b_im = tb_re * ξ_im + tb_im * ξ_re - ξ₁_im

        inv_den_b = one(FT) / (den_b_re^2 + den_b_im^2)
        bn[i, n] = Complex{FT}(
            (num_b_re * den_b_re + num_b_im * den_b_im) * inv_den_b,
            (num_b_im * den_b_re - num_b_re * den_b_im) * inv_den_b
        )

        # Shift for next iteration
        ψ₀ = ψ₁; ψ₁ = ψ
        χ₀ = χ₁; χ₁ = χ
        ξ₁_re = ψ₁; ξ₁_im = χ₁
    end
end

"""
    mie_coefficients_kernel_f64!(an, bn, x_params, m_re, m_im, nmax_per_r, nmx_max)

Native floating-point Mie coefficient kernel used by BOTH the `NativeFloat64`
policy (`Float64` arrays) and the `NativeFloat32` policy (`Float32` arrays):
despite the `_f64!` name (historical -- it predates `NativeFloat32`), every
scalar operation inside uses `eltype(x_params)`, so the kernel body never
hardcodes a floating-point type and is reused VERBATIM for plain native-Float32
arithmetic. The host launcher (`_mie_kernel1` in `gpu_precision.jl`) restricts
`NativeFloat64` to `Float64` arrays and `NativeFloat32` to `Float32` arrays.
"""
@kernel function mie_coefficients_kernel_f64!(
    an, bn,                     # output: (nquad_radius, nmax_global) Complex{FT}
    @Const(x_params),           # input:  (nquad_radius,) FT -- size parameters
    m_re, m_im,                 # scalars: real/imag refractive index parts
    @Const(nmax_per_r),         # input:  (nquad_radius,) Int -- per-radius nmax
    # FOLLOW-UP (dead parameter, pre-existing): see mie_coefficients_kernel_ds!'s
    # identical note -- `nmx_max` is never read in the body below either, and
    # is now propagated, unchanged, into mie_coefficients_kernel_f64_batched!.
    nmx_max::Int                # scalar: max Dn recursion depth
)
    i = @index(Global)
    FT = eltype(x_params)

    x = x_params[i]
    n_max_i = nmax_per_r[i]

    # Complex y in native FT. The host launcher only selects this kernel for
    # NativeFloat64 policy with Float64 arrays.
    y_re = FT(x) * FT(m_re)
    y_im = FT(x) * FT(m_im)
    nmx_i = max(n_max_i, round(Int, sqrt(y_re^2 + y_im^2)) + 1) + 50

    # --- Dn downward recursion in native FT ---
    Dn_prev_re = zero(FT)
    Dn_prev_im = zero(FT)

    @inbounds for n = (nmx_i - 1):-1:1
        # ratio = (n+1) / y  [complex division]
        n1 = FT(n + 1)
        denom = y_re^2 + y_im^2
        ratio_re = n1 * y_re / denom
        ratio_im = -n1 * y_im / denom

        # sum = Dn_prev + ratio
        sum_re = Dn_prev_re + ratio_re
        sum_im = Dn_prev_im + ratio_im

        # inv_sum = 1 / sum
        inv_denom = sum_re^2 + sum_im^2
        inv_re = sum_re / inv_denom
        inv_im = -sum_im / inv_denom

        # Dn_prev = ratio - inv_sum
        Dn_prev_re = ratio_re - inv_re
        Dn_prev_im = ratio_im - inv_im

        if n <= n_max_i
            # Store as native FT complex in bn (temporary)
            bn[i, n] = Complex{FT}(FT(Dn_prev_re), FT(Dn_prev_im))
        end
    end

    # --- Riccati-Bessel forward recursion (native FT) ---
    ψ₀ = cos(x)
    ψ₁ = sin(x)
    χ₀ = -sin(x)
    χ₁ = cos(x)
    ξ₁_re = ψ₁
    ξ₁_im = χ₁

    @inbounds for n = 1:n_max_i
        coef = FT(2n - 1) / x
        ψ = coef * ψ₁ - ψ₀
        χ = coef * χ₁ - χ₀
        ξ_re = ψ
        ξ_im = χ

        Dn_n = bn[i, n]
        Dn_re_n = real(Dn_n)
        Dn_im_n = imag(Dn_n)

        inv_m_denom = one(FT) / (m_re^2 + m_im^2)
        ta_re = (Dn_re_n * m_re + Dn_im_n * m_im) * inv_m_denom + FT(n) / x
        ta_im = (Dn_im_n * m_re - Dn_re_n * m_im) * inv_m_denom

        tb_re = Dn_re_n * m_re - Dn_im_n * m_im + FT(n) / x
        tb_im = Dn_re_n * m_im + Dn_im_n * m_re

        num_a_re = ta_re * ψ - ψ₁
        num_a_im = ta_im * ψ
        den_a_re = ta_re * ξ_re - ta_im * ξ_im - ξ₁_re
        den_a_im = ta_re * ξ_im + ta_im * ξ_re - ξ₁_im
        inv_den_a = one(FT) / (den_a_re^2 + den_a_im^2)
        an[i, n] = Complex{FT}(
            (num_a_re * den_a_re + num_a_im * den_a_im) * inv_den_a,
            (num_a_im * den_a_re - num_a_re * den_a_im) * inv_den_a
        )

        num_b_re = tb_re * ψ - ψ₁
        num_b_im = tb_im * ψ
        den_b_re = tb_re * ξ_re - tb_im * ξ_im - ξ₁_re
        den_b_im = tb_re * ξ_im + tb_im * ξ_re - ξ₁_im
        inv_den_b = one(FT) / (den_b_re^2 + den_b_im^2)
        bn[i, n] = Complex{FT}(
            (num_b_re * den_b_re + num_b_im * den_b_im) * inv_den_b,
            (num_b_im * den_b_re - num_b_re * den_b_im) * inv_den_b
        )

        ψ₀ = ψ₁; ψ₁ = ψ
        χ₀ = χ₁; χ₁ = χ
        ξ₁_re = ψ₁; ξ₁_im = χ₁
    end
end

"""
    _mie_kernel1(policy::MiePrecisionPolicy, backend)

Select the Kernel-1 (Mie coefficients) KA kernel for `policy`, used by the
caller-node GPU preamble (`_prepare_node_mie_gpu`). `NativeFloat64` and
`NativeFloat32` both launch the SAME generic native kernel
(`mie_coefficients_kernel_f64!` -- see its docstring for why the name doesn't
restrict it to Float64); `DSEmulated` launches the double-single kernel.
"""
_mie_kernel1(::NativeFloat64, backend) = mie_coefficients_kernel_f64!(backend)
_mie_kernel1(::NativeFloat32, backend) = mie_coefficients_kernel_f64!(backend)
_mie_kernel1(::DSEmulated,    backend) = mie_coefficients_kernel_ds!(backend)

# ============================================================================
# Kernel 1 (BATCHED sibling): per-ensemble refractive index via node_to_ensemble
# ============================================================================
# These are SIBLINGS of mie_coefficients_kernel_f64!/_ds! (added for the
# grouped caller-node batched API, compute_NAI2_nodes_batched.jl), not
# modifications of them: the single-ensemble kernels are untouched so their
# existing bit-compat guarantees carry zero risk from this addition. Each
# batched kernel's body is IDENTICAL to its single-ensemble sibling except
# that the scalar `m_re, m_im` parameters become `m_re_ens[e], m_im_ens[e]`
# where `e = node_to_ensemble[i]` -- since every other quantity (`x`, `n_max_i`,
# and the per-thread-local `nmx_i` recursion depth) is already indexed per
# node/thread `i`, this substitution does not change a single FLOP for any
# node: a node's aₙ/bₙ depend only on ITS OWN `x`/`m_re`/`m_im`/`n_max_i`,
# never on which OTHER nodes/ensembles share the launch. This is what makes
# exact-`nmax` GROUPED batching bitwise-identical to the single-ensemble
# kernel launch, by construction (see the batched module's docstring).

"""
    mie_coefficients_kernel_ds_batched!(an, bn, x_params, node_to_ensemble,
        m_re_ens, m_im_ens, nmax_per_r, nmx_max)

Batched sibling of [`mie_coefficients_kernel_ds!`](@ref) (identical body,
per-ensemble refractive index via `node_to_ensemble`). See the "Kernel 1
(BATCHED sibling)" section docstring for the bit-compat argument.
"""
@kernel function mie_coefficients_kernel_ds_batched!(
    an, bn,                     # output: (n_group_nodes, nmax_global) Complex{FT}
    @Const(x_params),           # input:  (n_group_nodes,) FT -- size parameters
    @Const(node_to_ensemble),   # input:  (n_group_nodes,) Int -- LOCAL (1:E_group) ensemble id per node
    @Const(m_re_ens), @Const(m_im_ens),  # input: (E_group,) FT -- per-ensemble refractive index
    @Const(nmax_per_r),         # input:  (n_group_nodes,) Int -- per-node nmax
    nmx_max::Int                # scalar: unused inside (kept for signature parity, see single-ensemble sibling)
)
    i = @index(Global)
    FT = eltype(x_params)

    x = x_params[i]
    n_max_i = nmax_per_r[i]
    e = node_to_ensemble[i]
    m_re = m_re_ens[e]
    m_im = m_im_ens[e]

    y_re = x * m_re
    y_im = x * m_im
    nmx_i = max(n_max_i, round(Int, sqrt(y_re^2 + y_im^2)) + 1) + 50

    Dn_re = DoubleSingle{FT}(zero(FT))
    Dn_im = DoubleSingle{FT}(zero(FT))
    Dn_prev = ComplexDS(Dn_re, Dn_im)

    y_ds = cds_complex(FT(y_re), FT(y_im))

    @inbounds for n = (nmx_i - 1):-1:1
        n_plus_1 = DoubleSingle{FT}(FT(n + 1))
        ratio = cds_div(cds_from_real(n_plus_1), y_ds)

        sum_val = cds_add(Dn_prev, ratio)
        inv_sum = cds_inv(sum_val)
        Dn_prev = cds_sub(ratio, inv_sum)

        if n <= n_max_i
            dn_native = to_complex(Dn_prev)
            bn[i, n] = dn_native
        end
    end

    ψ₀ = cos(x)
    ψ₁ = sin(x)
    χ₀ = -sin(x)
    χ₁ = cos(x)
    ξ₁_re = ψ₁
    ξ₁_im = χ₁

    @inbounds for n = 1:n_max_i
        coef = FT(2n - 1) / x
        ψ = coef * ψ₁ - ψ₀
        χ = coef * χ₁ - χ₀
        ξ_re = ψ
        ξ_im = χ

        Dn_n = bn[i, n]
        Dn_re_n = real(Dn_n)
        Dn_im_n = imag(Dn_n)

        inv_m_denom = one(FT) / (m_re^2 + m_im^2)
        ta_re = (Dn_re_n * m_re + Dn_im_n * m_im) * inv_m_denom + FT(n) / x
        ta_im = (Dn_im_n * m_re - Dn_re_n * m_im) * inv_m_denom

        tb_re = Dn_re_n * m_re - Dn_im_n * m_im + FT(n) / x
        tb_im = Dn_re_n * m_im + Dn_im_n * m_re

        num_a_re = ta_re * ψ - ψ₁
        num_a_im = ta_im * ψ
        den_a_re = ta_re * ξ_re - ta_im * ξ_im - ξ₁_re
        den_a_im = ta_re * ξ_im + ta_im * ξ_re - ξ₁_im
        inv_den_a = one(FT) / (den_a_re^2 + den_a_im^2)
        an[i, n] = Complex{FT}(
            (num_a_re * den_a_re + num_a_im * den_a_im) * inv_den_a,
            (num_a_im * den_a_re - num_a_re * den_a_im) * inv_den_a
        )

        num_b_re = tb_re * ψ - ψ₁
        num_b_im = tb_im * ψ
        den_b_re = tb_re * ξ_re - tb_im * ξ_im - ξ₁_re
        den_b_im = tb_re * ξ_im + tb_im * ξ_re - ξ₁_im
        inv_den_b = one(FT) / (den_b_re^2 + den_b_im^2)
        bn[i, n] = Complex{FT}(
            (num_b_re * den_b_re + num_b_im * den_b_im) * inv_den_b,
            (num_b_im * den_b_re - num_b_re * den_b_im) * inv_den_b
        )

        ψ₀ = ψ₁; ψ₁ = ψ
        χ₀ = χ₁; χ₁ = χ
        ξ₁_re = ψ₁; ξ₁_im = χ₁
    end
end

"""
    mie_coefficients_kernel_f64_batched!(an, bn, x_params, node_to_ensemble,
        m_re_ens, m_im_ens, nmax_per_r, nmx_max)

Batched sibling of [`mie_coefficients_kernel_f64!`](@ref) (identical body,
per-ensemble refractive index via `node_to_ensemble`) -- backs BOTH
`NativeFloat64` and `NativeFloat32` in the batched path, exactly as its
single-ensemble sibling does for the non-batched path (see
[`_mie_kernel1`](@ref) / [`_mie_kernel1_batched`](@ref)).
"""
@kernel function mie_coefficients_kernel_f64_batched!(
    an, bn,                     # output: (n_group_nodes, nmax_global) Complex{FT}
    @Const(x_params),           # input:  (n_group_nodes,) FT
    @Const(node_to_ensemble),   # input:  (n_group_nodes,) Int -- LOCAL (1:E_group) ensemble id per node
    @Const(m_re_ens), @Const(m_im_ens),  # input: (E_group,) FT
    @Const(nmax_per_r),         # input:  (n_group_nodes,) Int
    nmx_max::Int                # scalar: unused inside (kept for signature parity)
)
    i = @index(Global)
    FT = eltype(x_params)

    x = x_params[i]
    n_max_i = nmax_per_r[i]
    e = node_to_ensemble[i]
    m_re = m_re_ens[e]
    m_im = m_im_ens[e]

    y_re = FT(x) * FT(m_re)
    y_im = FT(x) * FT(m_im)
    nmx_i = max(n_max_i, round(Int, sqrt(y_re^2 + y_im^2)) + 1) + 50

    Dn_prev_re = zero(FT)
    Dn_prev_im = zero(FT)

    @inbounds for n = (nmx_i - 1):-1:1
        n1 = FT(n + 1)
        denom = y_re^2 + y_im^2
        ratio_re = n1 * y_re / denom
        ratio_im = -n1 * y_im / denom

        sum_re = Dn_prev_re + ratio_re
        sum_im = Dn_prev_im + ratio_im

        inv_denom = sum_re^2 + sum_im^2
        inv_re = sum_re / inv_denom
        inv_im = -sum_im / inv_denom

        Dn_prev_re = ratio_re - inv_re
        Dn_prev_im = ratio_im - inv_im

        if n <= n_max_i
            bn[i, n] = Complex{FT}(FT(Dn_prev_re), FT(Dn_prev_im))
        end
    end

    ψ₀ = cos(x)
    ψ₁ = sin(x)
    χ₀ = -sin(x)
    χ₁ = cos(x)
    ξ₁_re = ψ₁
    ξ₁_im = χ₁

    @inbounds for n = 1:n_max_i
        coef = FT(2n - 1) / x
        ψ = coef * ψ₁ - ψ₀
        χ = coef * χ₁ - χ₀
        ξ_re = ψ
        ξ_im = χ

        Dn_n = bn[i, n]
        Dn_re_n = real(Dn_n)
        Dn_im_n = imag(Dn_n)

        inv_m_denom = one(FT) / (m_re^2 + m_im^2)
        ta_re = (Dn_re_n * m_re + Dn_im_n * m_im) * inv_m_denom + FT(n) / x
        ta_im = (Dn_im_n * m_re - Dn_re_n * m_im) * inv_m_denom

        tb_re = Dn_re_n * m_re - Dn_im_n * m_im + FT(n) / x
        tb_im = Dn_re_n * m_im + Dn_im_n * m_re

        num_a_re = ta_re * ψ - ψ₁
        num_a_im = ta_im * ψ
        den_a_re = ta_re * ξ_re - ta_im * ξ_im - ξ₁_re
        den_a_im = ta_re * ξ_im + ta_im * ξ_re - ξ₁_im
        inv_den_a = one(FT) / (den_a_re^2 + den_a_im^2)
        an[i, n] = Complex{FT}(
            (num_a_re * den_a_re + num_a_im * den_a_im) * inv_den_a,
            (num_a_im * den_a_re - num_a_re * den_a_im) * inv_den_a
        )

        num_b_re = tb_re * ψ - ψ₁
        num_b_im = tb_im * ψ
        den_b_re = tb_re * ξ_re - tb_im * ξ_im - ξ₁_re
        den_b_im = tb_re * ξ_im + tb_im * ξ_re - ξ₁_im
        inv_den_b = one(FT) / (den_b_re^2 + den_b_im^2)
        bn[i, n] = Complex{FT}(
            (num_b_re * den_b_re + num_b_im * den_b_im) * inv_den_b,
            (num_b_im * den_b_re - num_b_re * den_b_im) * inv_den_b
        )

        ψ₀ = ψ₁; ψ₁ = ψ
        χ₀ = χ₁; χ₁ = χ
        ξ₁_re = ψ₁; ξ₁_im = χ₁
    end
end

"""
    _mie_kernel1_batched(policy::MiePrecisionPolicy, backend)

Batched-path counterpart of [`_mie_kernel1`](@ref): same policy → kernel
mapping, but selecting the `_batched!` sibling kernels (per-ensemble
refractive index via `node_to_ensemble`, see the "Kernel 1 (BATCHED sibling)"
section).
"""
_mie_kernel1_batched(::NativeFloat64, backend) = mie_coefficients_kernel_f64_batched!(backend)
_mie_kernel1_batched(::NativeFloat32, backend) = mie_coefficients_kernel_f64_batched!(backend)
_mie_kernel1_batched(::DSEmulated,    backend) = mie_coefficients_kernel_ds_batched!(backend)

# ============================================================================
# Kernel 2+3 (fused): Amplitude functions S1,S2 + phase matrix elements
# ============================================================================

"""
    amplitude_phase_kernel!(S1, S2, f11, f33, f12, f34, an, bn, leg_pi, leg_tau,
                            x_params, nmax_per_r)

Fused Kernel 2+3: Compute S1, S2 amplitude functions with Neumaier-compensated
summation over l, then immediately compute phase matrix elements.

Grid: (n_mu, nquad_radius) -- one thread per (angle, radius) pair.
"""
@kernel function amplitude_phase_kernel!(
    f11, f33, f12, f34,         # output: (n_mu, nquad_radius) FT
    @Const(an), @Const(bn),     # input:  (nquad_radius, nmax_global) Complex{FT}
    @Const(leg_pi), @Const(leg_tau),  # input: (n_mu, nmax_global) FT
    @Const(x_params),           # input:  (nquad_radius,) FT
    @Const(nmax_per_r)          # input:  (nquad_radius,) Int
)
    iμ, ir = @index(Global, NTuple)
    FT = eltype(f11)

    n_max_i = nmax_per_r[ir]
    x = x_params[ir]

    # Neumaier-compensated accumulation of S1, S2
    S1_acc = ComplexNeumaier{FT}()
    S2_acc = ComplexNeumaier{FT}()

    @inbounds for l = 1:n_max_i
        prefac = FT(2l + 1) / FT(l * (l + 1))

        an_l = an[ir, l]
        bn_l = bn[ir, l]
        pi_l = leg_pi[iμ, l]
        tau_l = leg_tau[iμ, l]

        # S1 += prefac * (an * tau + bn * pi)
        term1 = prefac * (an_l * tau_l + bn_l * pi_l)
        S1_acc = cneumaier_add(S1_acc, real(term1), imag(term1))

        # S2 += prefac * (an * pi + bn * tau)
        term2 = prefac * (an_l * pi_l + bn_l * tau_l)
        S2_acc = cneumaier_add(S2_acc, real(term2), imag(term2))
    end

    s1 = cneumaier_sum(S1_acc)
    s2 = cneumaier_sum(S2_acc)

    # Phase matrix elements (Tier 3 -- native FT)
    inv_x2 = convert(FT, 1//2) / (x * x)
    abs2_s1 = real(s1)^2 + imag(s1)^2
    abs2_s2 = real(s2)^2 + imag(s2)^2
    cross_re = real(s1) * real(s2) + imag(s1) * imag(s2)  # Re(s1 * conj(s2))
    cross_im = imag(s1) * real(s2) - real(s1) * imag(s2)  # Im(s1 * conj(s2))

    @inbounds f11[iμ, ir] =  inv_x2 * (abs2_s1 + abs2_s2)
    @inbounds f33[iμ, ir] =  inv_x2 * FT(2) * cross_re
    @inbounds f12[iμ, ir] = -inv_x2 * (abs2_s1 - abs2_s2)
    @inbounds f34[iμ, ir] = -inv_x2 * FT(2) * cross_im
end

# ============================================================================
# Kernel 4: Cross-sections + size-distribution reduction
# ============================================================================

"""
    cross_sections_kernel!(C_sca, C_ext, an, bn, x_params, k_wavenum, nmax_per_r)

Compute per-radius extinction and scattering cross-sections with Neumaier
compensation. One thread per radius point.
"""
@kernel function cross_sections_kernel!(
    C_sca, C_ext,               # output: (nquad_radius,) FT
    @Const(an), @Const(bn),     # input:  (nquad_radius, nmax_global) Complex{FT}
    k_wavenum,                  # scalar: 2pi/lambda
    @Const(nmax_per_r)          # input:  (nquad_radius,) Int
)
    i = @index(Global)
    FT = eltype(C_sca)

    n_max_i = nmax_per_r[i]

    sca_acc = NeumaierAccum{FT}()
    ext_acc = NeumaierAccum{FT}()

    @inbounds for n = 1:n_max_i
        w = FT(2n + 1)
        an_n = an[i, n]
        bn_n = bn[i, n]

        sca_acc = neumaier_add(sca_acc, w * (real(an_n)^2 + imag(an_n)^2 + real(bn_n)^2 + imag(bn_n)^2))
        ext_acc = neumaier_add(ext_acc, w * (real(an_n) + real(bn_n)))
    end

    prefac = FT(2) * FT(pi) / (k_wavenum^2)
    @inbounds C_sca[i] = prefac * neumaier_sum(sca_acc)
    @inbounds C_ext[i] = prefac * neumaier_sum(ext_acc)
end

"""
    size_reduction_kernel!(bulk_f, f_matrix, wr, n_mu)

Reduce phase matrix over size distribution: bulk_f[imu] = sum_i(f[imu, i] * wr[i]).
Uses Neumaier compensation for the 10-order-of-magnitude weight range.
One thread per angle point.
"""
@kernel function size_reduction_kernel!(
    bulk_f,                     # output: (n_mu,) FT
    @Const(f_matrix),           # input:  (n_mu, nquad_radius) FT
    @Const(wr)                  # input:  (nquad_radius,) FT -- 4pi*r^2*wx weights
)
    iμ = @index(Global)
    FT = eltype(bulk_f)
    nquad = size(f_matrix, 2)

    acc = NeumaierAccum{FT}()
    @inbounds for i = 1:nquad
        acc = neumaier_add(acc, f_matrix[iμ, i] * wr[i])
    end
    @inbounds bulk_f[iμ] = neumaier_sum(acc)
end

"""
    weighted_sum_kernel!(result, values, weights, n)

Compute result[1] = sum(values .* weights) with Neumaier compensation.
Single-thread kernel (or could be a parallel reduction).
"""
@kernel function weighted_sum_kernel!(
    result,                     # output: (1,) FT
    @Const(values),             # input:  (n,) FT
    @Const(weights)             # input:  (n,) FT
)
    FT = eltype(result)
    n = length(values)
    acc = NeumaierAccum{FT}()
    @inbounds for i = 1:n
        acc = neumaier_add(acc, values[i] * weights[i])
    end
    @inbounds result[1] = neumaier_sum(acc)
end

# ============================================================================
# Kernel 5: Greek coefficients
# ============================================================================

"""
    greek_coefficients_kernel!(greek_out, bulk_f11, bulk_f33, bulk_f12, bulk_f34,
                                P, P2, R2, T2, w_mu, l_max)

Compute Greek coefficients (alpha, beta, gamma, delta, epsilon, zeta) for each
Legendre order l. Neumaier-compensated dot products.

Grid: (l_max,) -- one thread per l value.
"""
@kernel function greek_coefficients_kernel!(
    alpha, beta, gamma, delta, epsilon, zeta,  # output: (l_max,) FT
    @Const(bulk_f11), @Const(bulk_f33),        # input:  (n_mu,) FT
    @Const(bulk_f12), @Const(bulk_f34),        # input:  (n_mu,) FT
    @Const(P), @Const(P2),                     # input:  (n_mu, l_max) FT
    @Const(R2), @Const(T2),                    # input:  (n_mu, l_max) FT
    @Const(w_mu),                              # input:  (n_mu,) FT
    n_mu::Int
)
    l_idx = @index(Global)  # 1-based, maps to l = l_idx - 1
    FT = eltype(alpha)
    l = l_idx - 1  # actual l value (0-based)

    # Pre-factor for l >= 2
    fac = l >= 2 ? FT(2l + 1) / FT(2) * sqrt(one(FT) / FT((l - 1) * l * (l + 1) * (l + 2))) : zero(FT)
    fac_beta = FT(2l + 1) / FT(2)

    # Neumaier-compensated dot products
    acc_beta  = NeumaierAccum{FT}()
    acc_delta = NeumaierAccum{FT}()
    acc_alpha = NeumaierAccum{FT}()
    acc_gamma = NeumaierAccum{FT}()
    acc_eps   = NeumaierAccum{FT}()
    acc_zeta  = NeumaierAccum{FT}()

    @inbounds for iμ = 1:n_mu
        w = w_mu[iμ]
        f11 = bulk_f11[iμ]
        f33 = bulk_f33[iμ]
        f12 = bulk_f12[iμ]
        f34 = bulk_f34[iμ]
        p   = P[iμ, l_idx]
        p2  = P2[iμ, l_idx]
        r2  = R2[iμ, l_idx]
        t2  = T2[iμ, l_idx]

        acc_beta  = neumaier_add(acc_beta,  w * f11 * p)
        acc_delta = neumaier_add(acc_delta, w * f33 * p)
        acc_gamma = neumaier_add(acc_gamma, w * f12 * p2)
        acc_eps   = neumaier_add(acc_eps,   w * f34 * p2)
        acc_alpha = neumaier_add(acc_alpha, w * (f11 * r2 + f33 * t2))
        acc_zeta  = neumaier_add(acc_zeta,  w * (f33 * r2 + f11 * t2))
    end

    @inbounds beta[l_idx]    = fac_beta * neumaier_sum(acc_beta)
    @inbounds delta[l_idx]   = fac_beta * neumaier_sum(acc_delta)
    @inbounds gamma[l_idx]   = fac * neumaier_sum(acc_gamma)
    @inbounds epsilon[l_idx] = fac * neumaier_sum(acc_eps)
    @inbounds alpha[l_idx]   = fac * neumaier_sum(acc_alpha)
    @inbounds zeta[l_idx]    = fac * neumaier_sum(acc_zeta)
end

# ============================================================================
# Kernel 4/5 (SEGMENTED siblings): per-ensemble reduction within a group
# ============================================================================
# Used by the grouped caller-node batched API (compute_NAI2_nodes_batched.jl).
# Each kernel below reduces one ensemble's node RANGE within a group's
# concatenated arrays with a thread PER (ensemble[, angle/l]) -- the exact same
# serial ascending-index Neumaier loop as the single-ensemble kernel it mirrors
# (weighted_sum_kernel!/size_reduction_kernel!/greek_coefficients_kernel!),
# just over `offsets[e]+1:offsets[e+1]` (LOCAL, group-relative offsets)
# instead of `1:n`. Segmentation changes ONLY the loop bounds -- the visited
# (value, weight) PAIRS, in order, are identical to running the single-
# ensemble kernel on that ensemble's own node sub-array, so bitwise equality
# to the single-ensemble path follows by the same argument as Kernel 1's
# batched siblings above.

"""
    segmented_weighted_sum_kernel!(result, values, weights, offsets)

Batched/segmented sibling of [`weighted_sum_kernel!`](@ref): one thread per
ensemble (LOCAL index `e` within the group, `1:length(offsets)-1`), each
computing `result[e] = sum(values[lo:hi] .* weights[lo:hi])`
(`lo, hi = offsets[e]+1, offsets[e+1]`) with the SAME Neumaier compensation
and ascending-index loop order as the single-ensemble kernel.
"""
@kernel function segmented_weighted_sum_kernel!(
    result,                     # output: (n_ensembles_in_group,) FT
    @Const(values),             # input:  (n_group_nodes,) FT
    @Const(weights),            # input:  (n_group_nodes,) FT
    @Const(offsets)             # input:  (n_ensembles_in_group+1,) Int -- LOCAL, offsets[1]=0
)
    e = @index(Global)
    FT = eltype(result)

    @inbounds lo = offsets[e] + 1
    @inbounds hi = offsets[e + 1]

    acc = NeumaierAccum{FT}()
    @inbounds for i = lo:hi
        acc = neumaier_add(acc, values[i] * weights[i])
    end
    @inbounds result[e] = neumaier_sum(acc)
end

"""
    segmented_size_reduction_kernel!(bulk_f, f_matrix, wr, offsets)

Batched/segmented sibling of [`size_reduction_kernel!`](@ref): one thread per
`(ensemble, angle)` pair, `bulk_f[iμ, e] = sum(f_matrix[iμ, lo:hi] .* wr[lo:hi])`
(`lo, hi = offsets[e]+1, offsets[e+1]`), same Neumaier compensation and
ascending-index loop order as the single-ensemble kernel.
"""
@kernel function segmented_size_reduction_kernel!(
    bulk_f,                     # output: (n_mu, n_ensembles_in_group) FT
    @Const(f_matrix),           # input:  (n_mu, n_group_nodes) FT
    @Const(wr),                 # input:  (n_group_nodes,) FT -- 4pi*r^2*w weights
    @Const(offsets)             # input:  (n_ensembles_in_group+1,) Int -- LOCAL, offsets[1]=0
)
    iμ, e = @index(Global, NTuple)
    FT = eltype(bulk_f)

    @inbounds lo = offsets[e] + 1
    @inbounds hi = offsets[e + 1]

    acc = NeumaierAccum{FT}()
    @inbounds for i = lo:hi
        acc = neumaier_add(acc, f_matrix[iμ, i] * wr[i])
    end
    @inbounds bulk_f[iμ, e] = neumaier_sum(acc)
end

"""
    greek_coefficients_kernel_batched!(alpha, beta, gamma, delta, epsilon, zeta,
        bulk_f11, bulk_f33, bulk_f12, bulk_f34, P, P2, R2, T2, w_mu, n_mu)

Batched sibling of [`greek_coefficients_kernel!`](@ref): one thread per
`(l, ensemble)` pair. `P`/`P2`/`R2`/`T2`/`w_mu` (the Legendre/quadrature
tables) are shared across the WHOLE group -- valid because grouping is by
exact `n_max_global` (hence exact shared `n_mu`), so every ensemble in the
group has the identical angular grid the single-ensemble path would compute
for it. `bulk_f11`/`bulk_f33`/`bulk_f12`/`bulk_f34` are per-ensemble bulk
phase-matrix arrays (output of [`segmented_size_reduction_kernel!`](@ref));
outputs are `(l_max, n_ensembles_in_group)`. Per-`(l, e)` arithmetic is
IDENTICAL to the single-ensemble kernel's per-`l` arithmetic.
"""
@kernel function greek_coefficients_kernel_batched!(
    alpha, beta, gamma, delta, epsilon, zeta,  # output: (l_max, n_ensembles_in_group) FT
    @Const(bulk_f11), @Const(bulk_f33),        # input:  (n_mu, n_ensembles_in_group) FT
    @Const(bulk_f12), @Const(bulk_f34),        # input:  (n_mu, n_ensembles_in_group) FT
    @Const(P), @Const(P2),                     # input:  (n_mu, l_max) FT -- shared across the group
    @Const(R2), @Const(T2),                    # input:  (n_mu, l_max) FT -- shared across the group
    @Const(w_mu),                              # input:  (n_mu,) FT -- shared across the group
    n_mu::Int
)
    l_idx, e = @index(Global, NTuple)  # l_idx 1-based -> l = l_idx - 1
    FT = eltype(alpha)
    l = l_idx - 1

    fac = l >= 2 ? FT(2l + 1) / FT(2) * sqrt(one(FT) / FT((l - 1) * l * (l + 1) * (l + 2))) : zero(FT)
    fac_beta = FT(2l + 1) / FT(2)

    acc_beta  = NeumaierAccum{FT}()
    acc_delta = NeumaierAccum{FT}()
    acc_alpha = NeumaierAccum{FT}()
    acc_gamma = NeumaierAccum{FT}()
    acc_eps   = NeumaierAccum{FT}()
    acc_zeta  = NeumaierAccum{FT}()

    @inbounds for iμ = 1:n_mu
        w = w_mu[iμ]
        f11 = bulk_f11[iμ, e]
        f33 = bulk_f33[iμ, e]
        f12 = bulk_f12[iμ, e]
        f34 = bulk_f34[iμ, e]
        p   = P[iμ, l_idx]
        p2  = P2[iμ, l_idx]
        r2  = R2[iμ, l_idx]
        t2  = T2[iμ, l_idx]

        acc_beta  = neumaier_add(acc_beta,  w * f11 * p)
        acc_delta = neumaier_add(acc_delta, w * f33 * p)
        acc_gamma = neumaier_add(acc_gamma, w * f12 * p2)
        acc_eps   = neumaier_add(acc_eps,   w * f34 * p2)
        acc_alpha = neumaier_add(acc_alpha, w * (f11 * r2 + f33 * t2))
        acc_zeta  = neumaier_add(acc_zeta,  w * (f33 * r2 + f11 * t2))
    end

    @inbounds beta[l_idx, e]    = fac_beta * neumaier_sum(acc_beta)
    @inbounds delta[l_idx, e]   = fac_beta * neumaier_sum(acc_delta)
    @inbounds gamma[l_idx, e]   = fac * neumaier_sum(acc_gamma)
    @inbounds epsilon[l_idx, e] = fac * neumaier_sum(acc_eps)
    @inbounds alpha[l_idx, e]   = fac * neumaier_sum(acc_alpha)
    @inbounds zeta[l_idx, e]    = fac * neumaier_sum(acc_zeta)
end

# ============================================================================
# Kernel 0 (batched-only): per-node Mie-series truncation order, device-side
# ============================================================================
# Used by the grouped batched API (compute_NAI2_nodes_batched.jl) so that,
# after a device-side permutation gather assembles a group's `x_size_param`
# array directly on-device (see that module's docstring for why -- avoiding a
# host round-trip for the float node payload), the per-node `n_max_i` each
# node needs (for Kernel 1's per-node loop bound and for sizing later
# kernels) can ALSO be computed on-device, without copying `x_size_param`
# back to the host at all. Calls the existing scalar `get_n_max` (defined in
# mie_helper_functions.jl) elementwise -- same one function, same values,
# whether invoked from host code or from this kernel.

"""
    nmax_per_node_kernel!(nmax_out, x_params)

`nmax_out[i] = get_n_max(x_params[i])`, one thread per node. Device-side
counterpart of the host-side broadcast `get_n_max.(x_size_param)` used
throughout the rest of this module.
"""
@kernel function nmax_per_node_kernel!(
    nmax_out,           # output: (n,) Int
    @Const(x_params)    # input:  (n,) FT
)
    i = @index(Global)
    @inbounds nmax_out[i] = get_n_max(x_params[i])
end
