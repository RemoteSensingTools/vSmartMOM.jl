#=
GPU-accelerated NAI2 computation of aerosol optical properties.

Entry point: compute_aerosol_optical_properties_gpu(model, backend; precision_policy, FT2)

This does NOT modify existing CPU code paths.  It provides a parallel GPU
implementation that produces the same AerosolOptics output as the CPU NAI2
path, using the kernel decomposition from gpu_mie_kernels.jl and precision
primitives from gpu_precision.jl.
=#

using KernelAbstractions
using FastGaussQuadrature

"""
    compute_aerosol_optical_properties_gpu(model::MieModel{<:NAI2}, backend;
                                           precision_policy=NativeFloat64(),
                                           FT2=Float64)

GPU-accelerated version of `compute_aerosol_optical_properties` for the NAI2 method.

# Arguments
- `model`: MieModel configured for NAI2 computation
- `backend`: KernelAbstractions backend (e.g., `CPU()`, `CUDABackend()`)
- `precision_policy`: `NativeFloat64()` for A100/V100, `DSEmulated()` for L40S
- `FT2`: output float type (default Float64)

# Returns
- `AerosolOptics` with greek_coefs, SSA, extinction, truncation factor
"""
function compute_aerosol_optical_properties_gpu(
    model::MieModel{FDT}, backend;
    precision_policy::MiePrecisionPolicy = NativeFloat64(),
    FT2::Type = Float64
) where FDT <: NAI2

    # Unpack the model
    (; aerosol, λ, r_max, nquad_radius) = model
    (; size_distribution, nᵣ, nᵢ) = aerosol

    @assert nᵢ >= 0

    FT = eltype(nᵣ)

    # --- Setup (same as CPU path) ---
    r_min = max(quantile(size_distribution, 1e-8), 1e-6 * r_max)
    r, wᵣ = gauleg_log(nquad_radius, r_min, r_max; norm=false)

    k_wavenum = FT(2π / λ)
    x_size_param = k_wavenum .* r

    n_max_global = get_n_max(maximum(x_size_param))
    n_mu = 2n_max_global - 1

    μ, w_μ = gausslegendre(n_mu)
    leg_π, leg_τ = compute_mie_π_τ(FT.(μ), n_max_global)

    wₓ = compute_wₓ(size_distribution, wᵣ, r, r_max)

    m_ref_re = FT(nᵣ)
    m_ref_im = FT(-nᵢ)  # Convention: n = nr - i*ni, so m_im is negative

    # Per-radius nmax values
    nmax_per_r = [get_n_max(x) for x in x_size_param]

    # nmx_max for Dn recursion depth
    y_max = maximum(x_size_param) * sqrt(m_ref_re^2 + m_ref_im^2)
    nmx_max = round(Int, max(n_max_global, y_max) + 51)

    # --- Allocate device arrays ---
    AT = KernelAbstractions.allocate  # shorthand

    x_dev       = AT(backend, FT, nquad_radius)
    nmax_dev    = AT(backend, Int, nquad_radius)
    an_dev      = AT(backend, Complex{FT}, (nquad_radius, n_max_global))
    bn_dev      = AT(backend, Complex{FT}, (nquad_radius, n_max_global))
    f11_dev     = AT(backend, FT, (n_mu, nquad_radius))
    f33_dev     = AT(backend, FT, (n_mu, nquad_radius))
    f12_dev     = AT(backend, FT, (n_mu, nquad_radius))
    f34_dev     = AT(backend, FT, (n_mu, nquad_radius))
    C_sca_dev   = AT(backend, FT, nquad_radius)
    C_ext_dev   = AT(backend, FT, nquad_radius)
    leg_pi_dev  = AT(backend, FT, (n_mu, n_max_global))
    leg_tau_dev = AT(backend, FT, (n_mu, n_max_global))

    # Copy data to device
    KernelAbstractions.copyto!(backend, x_dev, FT.(x_size_param))
    KernelAbstractions.copyto!(backend, nmax_dev, nmax_per_r)
    KernelAbstractions.copyto!(backend, leg_pi_dev, FT.(leg_π))
    KernelAbstractions.copyto!(backend, leg_tau_dev, FT.(leg_τ))

    # Zero output arrays
    fill!(an_dev, zero(Complex{FT}))
    fill!(bn_dev, zero(Complex{FT}))

    # --- Kernel 1: Mie coefficients ---
    if precision_policy isa NativeFloat64
        @assert FT === Float64 "NativeFloat64 Mie precision policy requires Float64 model inputs; use DSEmulated() for Float32."
        kernel1 = mie_coefficients_kernel_f64!(backend)
    else
        kernel1 = mie_coefficients_kernel_ds!(backend)
    end
    kernel1(an_dev, bn_dev, x_dev, m_ref_re, m_ref_im, nmax_dev, nmx_max;
            ndrange=nquad_radius)
    KernelAbstractions.synchronize(backend)

    # --- Kernel 2+3: Amplitude functions + phase matrix ---
    kernel23 = amplitude_phase_kernel!(backend)
    kernel23(f11_dev, f33_dev, f12_dev, f34_dev,
             an_dev, bn_dev, leg_pi_dev, leg_tau_dev,
             x_dev, nmax_dev;
             ndrange=(n_mu, nquad_radius))
    KernelAbstractions.synchronize(backend)

    # --- Kernel 4a: Cross-sections ---
    kernel4a = cross_sections_kernel!(backend)
    kernel4a(C_sca_dev, C_ext_dev, an_dev, bn_dev, k_wavenum, nmax_dev;
             ndrange=nquad_radius)
    KernelAbstractions.synchronize(backend)

    # --- Copy back to CPU for reduction + Greek coefficients ---
    # (These are relatively small arrays, so the transfer is cheap)
    f11_host = Array(f11_dev)
    f33_host = Array(f33_dev)
    f12_host = Array(f12_dev)
    f34_host = Array(f34_dev)
    C_sca_host = Array(C_sca_dev)
    C_ext_host = Array(C_ext_dev)

    # --- Size distribution reduction (CPU, Neumaier-compensated) ---
    bulk_C_sca = neumaier_dot(C_sca_host, FT.(wₓ))
    bulk_C_ext = neumaier_dot(C_ext_host, FT.(wₓ))

    wr = FT.(4π .* r.^2 .* wₓ)

    # Compute bulk phase matrix via matrix-vector multiply with Neumaier
    bulk_f11 = neumaier_matvec(f11_host, wr)
    bulk_f33 = neumaier_matvec(f33_host, wr)
    bulk_f12 = neumaier_matvec(f12_host, wr)
    bulk_f34 = neumaier_matvec(f34_host, wr)

    # Normalize
    inv_bulk_C_sca = one(FT) / bulk_C_sca
    bulk_f11 .*= inv_bulk_C_sca
    bulk_f33 .*= inv_bulk_C_sca
    bulk_f12 .*= inv_bulk_C_sca
    bulk_f34 .*= inv_bulk_C_sca

    # --- Greek coefficients (CPU, Neumaier-compensated) ---
    l_max = n_mu
    P, P², R², T² = compute_legendre_poly(FT.(μ), l_max)

    α = zeros(FT, l_max)
    β = zeros(FT, l_max)
    γ = zeros(FT, l_max)
    δ = zeros(FT, l_max)
    ϵ = zeros(FT, l_max)
    ζ = zeros(FT, l_max)

    w_μ_ft = FT.(w_μ)
    for l = 0:l_max-1
        fac = l >= 2 ? FT(2l + 1) / FT(2) * sqrt(one(FT) / FT((l - 1) * l * (l + 1) * (l + 2))) : zero(FT)
        fac_beta = FT(2l + 1) / FT(2)

        β[l+1] = fac_beta * neumaier_dot_3(w_μ_ft, bulk_f11, view(P, :, l+1))
        δ[l+1] = fac_beta * neumaier_dot_3(w_μ_ft, bulk_f33, view(P, :, l+1))
        if l >= 2
            γ[l+1] = fac * neumaier_dot_3(w_μ_ft, bulk_f12, view(P², :, l+1))
            ϵ[l+1] = fac * neumaier_dot_3(w_μ_ft, bulk_f34, view(P², :, l+1))
            α[l+1] = fac * (neumaier_dot_3(w_μ_ft, bulk_f11, view(R², :, l+1)) +
                            neumaier_dot_3(w_μ_ft, bulk_f33, view(T², :, l+1)))
            ζ[l+1] = fac * (neumaier_dot_3(w_μ_ft, bulk_f33, view(R², :, l+1)) +
                            neumaier_dot_3(w_μ_ft, bulk_f11, view(T², :, l+1)))
        end
    end

    # --- Package output ---
    if FT <: AbstractFloat && FT2 <: AbstractFloat
        greek_coefs = GreekCoefs(convert.(FT2, α), convert.(FT2, β),
                                 convert.(FT2, γ), convert.(FT2, δ),
                                 convert.(FT2, ϵ), convert.(FT2, ζ))
        return AerosolOptics(greek_coefs=greek_coefs,
                             ω̃=FT2(bulk_C_sca / bulk_C_ext),
                             k=FT2(bulk_C_ext), fᵗ=FT2(1))
    else
        greek_coefs = GreekCoefs(α, β, γ, δ, ϵ, ζ)
        return AerosolOptics(greek_coefs=greek_coefs,
                             ω̃=(bulk_C_sca / bulk_C_ext),
                             k=(bulk_C_ext), fᵗ=one(eltype(α)))
    end
end

# ============================================================================
# CPU helper functions for Neumaier-compensated reductions
# ============================================================================

"""Neumaier-compensated weighted sum: sum(a .* b)."""
function neumaier_dot(a::AbstractVector{T}, b::AbstractVector{T}) where {T}
    acc = NeumaierAccum{T}()
    @inbounds for i in eachindex(a, b)
        acc = neumaier_add(acc, a[i] * b[i])
    end
    neumaier_sum(acc)
end

"""Neumaier-compensated triple product sum: sum(a .* b .* c)."""
function neumaier_dot_3(a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector) where {T}
    acc = NeumaierAccum{T}()
    @inbounds for i in eachindex(a, b, c)
        acc = neumaier_add(acc, a[i] * b[i] * T(c[i]))
    end
    neumaier_sum(acc)
end

"""Neumaier-compensated matrix-vector product: result[i] = sum_j(M[i,j] * v[j])."""
function neumaier_matvec(M::AbstractMatrix{T}, v::AbstractVector{T}) where {T}
    n_rows = size(M, 1)
    result = zeros(T, n_rows)
    @inbounds for i = 1:n_rows
        acc = NeumaierAccum{T}()
        for j in eachindex(v)
            acc = neumaier_add(acc, M[i, j] * v[j])
        end
        result[i] = neumaier_sum(acc)
    end
    result
end
