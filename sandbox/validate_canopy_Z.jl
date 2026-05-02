# =================================================================
# Validate the CanopyOptics v0.2 Z-matrix convention used by vSmartMOM.
#
# This checks the public analytic Fourier stack for
# BiLambertianCanopyScattering:
#   1. m=0 non-negativity.
#   2. single-moment API slices agree with the stacked API.
#   3. spherical-LAD m=0 agrees with the direct isotropic-kernel reference.
#   4. vSmartMOM normalization: sum_i w_i (Z++ + Z-+)[:, j] ~= 2.
#   5. Stokes expansion populates only I->I for diffuse bi-Lambertian leaves.
# =================================================================

using CanopyOptics
using LinearAlgebra
using Printf

const CO = CanopyOptics
const FT = Float64

R_leaf = FT(0.5)
T_leaf = FT(0.5)
ω = R_leaf + T_leaf
m_max = 6

scatter = CO.BiLambertianCanopyScattering(R = R_leaf, T = T_leaf)
LAD = CO.spherical_leaves(FT)
quadrature = CO.CanopyQuadrature(n_leaf = 64, n_azimuth = 48)

μ, w = CO.gauleg(20, zero(FT), one(FT))
μ = collect(μ)
w = collect(w)
nμ = length(μ)

println("BiLambertianCanopyScattering: R=$(R_leaf), T=$(T_leaf), ω=$(ω)")
println("LAD: spherical leaves (G(μ) = 0.5)")
println("μ grid: $(nμ) Gauss-Legendre points on (0, 1]")
println("Fourier moments: 0:$(m_max)")

Zpp, Zmp = CO.compute_Z_matrices_aniso_analytic(
    scatter, μ, LAD, m_max; quadrature)

println("\n--- m=0 non-negativity ---")
@printf("Z++ min = %.6e\n", minimum(Zpp[:, :, 1]))
@printf("Z-+ min = %.6e\n", minimum(Zmp[:, :, 1]))
@assert minimum(Zpp[:, :, 1]) >= -1e-12
@assert minimum(Zmp[:, :, 1]) >= -1e-12

println("\n--- Stacked API vs single-moment API ---")
for m in 0:m_max
    Zpp_m, Zmp_m = CO.compute_Z_matrices(scatter, μ, LAD, m; quadrature)
    err_pp = norm(Zpp[:, :, m + 1] - Zpp_m) / max(norm(Zpp_m), eps(FT))
    err_mp = norm(Zmp[:, :, m + 1] - Zmp_m) / max(norm(Zmp_m), eps(FT))
    @printf("m=%d  relerr(Z++) = %.3e  relerr(Z-+) = %.3e\n", m, err_pp, err_mp)
    @assert err_pp < 1e-12
    @assert err_mp < 1e-12
end

println("\n--- m=0 direct isotropic-kernel reference ---")
nφ = 256
φ = collect(range(zero(FT), FT(π); length = nφ + 1))[1:end-1]
dφ = FT(π) / nφ
wφ = (FT(2) / FT(π)) * dφ

Zpp_ref = zeros(FT, nμ, nμ)
Zmp_ref = zeros(FT, nμ, nμ)
scale = FT(2) / (ω * FT(0.5)) # spherical LAD has G(μ_in) = 0.5
for j in 1:nμ, i in 1:nμ
    Ωin = CO.dirVector_μ(μ[j], zero(FT))
    Γdown = zero(FT)
    Γup = zero(FT)
    for ϕ in φ
        Γdown += CO.compute_Γ_isotropic(scatter, Ωin, CO.dirVector_μ( μ[i], ϕ)) * wφ
        Γup   += CO.compute_Γ_isotropic(scatter, Ωin, CO.dirVector_μ(-μ[i], ϕ)) * wφ
    end
    Zpp_ref[i, j] = scale * Γdown
    Zmp_ref[i, j] = scale * Γup
end

err_ref_pp = norm(Zpp[:, :, 1] - Zpp_ref) / norm(Zpp_ref)
err_ref_mp = norm(Zmp[:, :, 1] - Zmp_ref) / norm(Zmp_ref)
@printf("relerr(Z++) = %.3e  relerr(Z-+) = %.3e\n", err_ref_pp, err_ref_mp)
@assert err_ref_pp < 5e-3
@assert err_ref_mp < 5e-3

println("\n--- vSmartMOM m=0 flux normalization ---")
fluxes = [sum(w .* (Zpp[:, j, 1] .+ Zmp[:, j, 1])) for j in eachindex(μ)]
@printf("flux range = [%.6f, %.6f], mean = %.6f, target = 2\n",
        minimum(fluxes), maximum(fluxes), sum(fluxes) / length(fluxes))
@assert all(abs.(fluxes .- FT(2)) .< FT(0.03))

println("\n--- Stokes expansion smoke ---")
npol = 4
Zpp4, Zmp4 = CO.compute_Z_matrices_aniso_analytic(
    scatter, μ, LAD, 2; quadrature, npol)
@assert size(Zpp4) == (npol * nμ, npol * nμ, 3)
@assert size(Zmp4) == size(Zpp4)
for si in 1:npol, sj in 1:npol
    (si, sj) == (1, 1) && continue
    @assert all(iszero, Zpp4[si:npol:end, sj:npol:end, :])
    @assert all(iszero, Zmp4[si:npol:end, sj:npol:end, :])
end

println("\nDone.")
