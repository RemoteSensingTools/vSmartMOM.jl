"""
exact_ss_reference.jl — standalone reference implementation of exact single
scattering through a plane-parallel atmosphere with a Lambertian surface.

Companion to derivations.md.

Purpose:
- Validate the architectural choices proposed for vSmartMOM v0.6
  (analytic-SS replacement, AbstractPostprocessCorrection unification).
- Reference for testing future ExactAnalyticSS implementations against an
  independent computation of the same physics.
- Demonstrates the contributor pattern's `exact_phase_function` trait
  without the full vSmartMOM type hierarchy.

Scope:
- Scalar Stokes I (no polarization).
- Lambertian surface (no BRDF).
- Single solar beam, azimuth 0.
- All four SS paths:
    1. sun → atmosphere → sensor
    2. sun → surface → sensor
    3. sun → atmosphere → surface → sensor
    4. sun → surface → atmosphere → sensor

Limitations (deliberate):
- No δ-M truncation: this *is* the reference, computed without truncation.
- Plane-parallel; no sphericity.
- Lambertian only; BRDF generalization left for future work.

Run from the REPL:
    julia> include("exact_ss_reference.jl")
    julia> ExactSSReference.run_all_tests()
"""
module ExactSSReference

using Test
using LinearAlgebra

# ============================================================================
# Contributor abstraction
# ============================================================================

abstract type AbstractContributor end

"""
    exact_phase_function(c::AbstractContributor, cosΘ) -> Float64

Phase function P(cosΘ) for scattering angle Θ.
Normalization: ∫ P(cosΘ) dΩ / 4π = 1, i.e.
∫_{-1}^{1} P(μ) dμ / 2 = 1.
"""
function exact_phase_function end

"""
    optical_depth(c::AbstractContributor, layer_idx) -> Float64

Optical depth contribution to layer `layer_idx`.
"""
function optical_depth end

"""
    single_scattering_albedo(c::AbstractContributor) -> Float64

ϖ. Pure scattering: 1. Pure absorption: 0.
"""
function single_scattering_albedo end

# ----------------------------------------------------------------------------
# Concrete contributors
# ----------------------------------------------------------------------------

"""Rayleigh: P = (3/4)(1 + cos²Θ), ϖ = 1."""
struct RayleighContributor <: AbstractContributor
    τ::Vector{Float64}
end

exact_phase_function(::RayleighContributor, cosΘ) = 0.75 * (1.0 + cosΘ^2)
optical_depth(c::RayleighContributor, iz) = c.τ[iz]
single_scattering_albedo(::RayleighContributor) = 1.0

"""Henyey-Greenstein aerosol: P = (1-g²)/(1+g²-2g·cosΘ)^(3/2)."""
struct HGAerosolContributor <: AbstractContributor
    g::Float64
    ϖ::Float64
    τ::Vector{Float64}
end

function exact_phase_function(c::HGAerosolContributor, cosΘ)
    g = c.g
    return (1.0 - g^2) / (1.0 + g^2 - 2g * cosΘ)^1.5
end
optical_depth(c::HGAerosolContributor, iz) = c.τ[iz]
single_scattering_albedo(c::HGAerosolContributor) = c.ϖ

"""Pure absorption."""
struct AbsorptionContributor <: AbstractContributor
    τ::Vector{Float64}
end

exact_phase_function(::AbsorptionContributor, cosΘ) = 0.0
optical_depth(c::AbsorptionContributor, iz) = c.τ[iz]
single_scattering_albedo(::AbsorptionContributor) = 0.0

# ============================================================================
# Layer-effective optical properties
# ============================================================================

"""Total layer optical depth (sum over contributors)."""
function layer_total_τ(contributors, iz)
    sum(optical_depth(c, iz) for c in contributors)
end

"""
Layer-effective single-scattering albedo:
    ϖ_eff = (∑ τ_i ϖ_i) / (∑ τ_i)
"""
function layer_effective_ϖ(contributors, iz)
    τ_total = layer_total_τ(contributors, iz)
    τ_total == 0 && return 0.0
    τ_scat = sum(optical_depth(c, iz) * single_scattering_albedo(c)
                 for c in contributors)
    return τ_scat / τ_total
end

"""
Layer-effective phase function at scattering angle Θ:
    P_eff(cosΘ) = (∑ τ_i ϖ_i P_i(cosΘ)) / (∑ τ_i ϖ_i)

Same physics as the existing CoreScatteringOpticalProperties `+` overload
in vSmartMOM types.jl, expressed in closed form.
"""
function layer_effective_phase(contributors, iz, cosΘ)
    weight_sum = 0.0
    weighted_phase = 0.0
    for c in contributors
        ϖ = single_scattering_albedo(c)
        ϖ == 0 && continue
        w = optical_depth(c, iz) * ϖ
        weight_sum += w
        weighted_phase += w * exact_phase_function(c, cosΘ)
    end
    weight_sum == 0 && return 0.0
    return weighted_phase / weight_sum
end

# ============================================================================
# Geometry
# ============================================================================

"""
Scattering angle cosine for sun → view geometry.

Sun direction (downward at zenith θ₀, azimuth ϕ₀=0): (sin θ₀, 0, -μ₀).
View direction (upward to sensor at θ_v, ϕ_v): (sin θ_v cos ϕ_v, sin θ_v sin ϕ_v, +μ_v).

cos Θ = -μ₀ μ_v + √(1-μ₀²) √(1-μ_v²) cos(Δϕ)
"""
function scattering_angle_cosine(μ₀, μ_v, Δϕ)
    return -μ₀ * μ_v + sqrt(1 - μ₀^2) * sqrt(1 - μ_v^2) * cos(Δϕ)
end

"""
Cumulative optical depth from TOA at top of each layer.
τ_cum[1] = 0 (TOA). τ_cum[N+1] = τ_total (surface).
"""
function cumulative_τ(contributors)
    n_layers = length(contributors[1].τ)
    τ_cum = zeros(n_layers + 1)
    for iz in 1:n_layers
        τ_cum[iz + 1] = τ_cum[iz] + layer_total_τ(contributors, iz)
    end
    return τ_cum
end

# ============================================================================
# PATH 1: sun → atmosphere → sensor
# ============================================================================

"""
Atmospheric SS radiance at TOA, summed over layers.

For each homogeneous layer between cumulative optical depths τ_top and τ_bot:

    ΔI_layer = (ϖ_eff I₀ / 4π) P_eff(cosΘ) ·
               [exp(-τ_top·a) - exp(-τ_bot·a)] / a

with a = 1/μ₀ + 1/μ_v.

Closed form. No quadratures.
"""
function path1_atmospheric_ss(contributors, μ₀, μ_v, Δϕ, I₀)
    @assert μ₀ > 0 && μ_v > 0 "Plane-parallel upward-viewing only"
    cosΘ = scattering_angle_cosine(μ₀, μ_v, Δϕ)
    τ_cum = cumulative_τ(contributors)
    n_layers = length(τ_cum) - 1
    a = 1.0/μ₀ + 1.0/μ_v
    # Factor 1/μ_v from the RTE solution: dτ_slant = dτ/μ_v along upward path.
    prefactor = I₀ / (4π * μ_v) / a

    L = 0.0
    for iz in 1:n_layers
        ϖ_eff = layer_effective_ϖ(contributors, iz)
        ϖ_eff == 0 && continue
        P_eff = layer_effective_phase(contributors, iz, cosΘ)
        layer_factor = exp(-τ_cum[iz] * a) - exp(-τ_cum[iz+1] * a)
        L += prefactor * ϖ_eff * P_eff * layer_factor
    end
    return L
end

# ============================================================================
# PATH 2: sun → surface → sensor
# ============================================================================

"""
Lambertian surface direct-beam contribution at TOA:

    L = (μ₀ I₀ albedo / π) · exp(-τ_total/μ₀) · exp(-τ_total/μ_v)
"""
function path2_surface_direct(τ_total, μ₀, μ_v, albedo, I₀)
    return (μ₀ * I₀ * albedo / π) * exp(-τ_total/μ₀) * exp(-τ_total/μ_v)
end

# ============================================================================
# PATHS 3 and 4: atm-surface coupled
# ============================================================================

"""
Gauss-Legendre nodes and weights on (0, 1].
Returns (nodes, weights) for ∫₀¹ f(x) dx ≈ ∑ w_i f(x_i).
"""
function gauss_legendre_01(n)
    # Standard Gauss-Legendre on [-1, 1], then map to [0, 1].
    # We use the eigenvalue-of-Jacobi-matrix approach for arbitrary n.
    β = [k / sqrt(4*k^2 - 1) for k in 1:(n-1)]
    T = SymTridiagonal(zeros(n), β)
    eig = eigen(T)
    x_std = eig.values
    w_std = 2 .* eig.vectors[1, :].^2
    # Map [-1, 1] -> [0, 1]:
    x = 0.5 .* (x_std .+ 1.0)
    w = 0.5 .* w_std
    return x, w
end

"""
Azimuthally-averaged phase function:
    P̄(μ', μ₀) = (1/2π) ∫₀^{2π} P(cos Θ(μ', μ₀, ϕ')) dϕ'

For Rayleigh, closed form:
    P̄_Ray(μ', μ₀) = 0.75(1 + (μ₀ μ')² + 0.5(1-μ₀²)(1-μ'²))

For HG and other phase functions, computed by 1D ϕ' quadrature.

Uses N_phi-point trapezoidal rule on [0, 2π] which converges exponentially
fast for periodic smooth integrands.
"""
function azimuthal_average_phase(c::AbstractContributor, μ_d, μ₀; N_phi=64)
    # μ_d > 0 is downward cosine (μ' = -μ_d in our convention).
    # cos Θ = -μ₀·μ' + √(1-μ₀²)√(1-μ'²) cos ϕ'
    #        = μ₀·μ_d + √(1-μ₀²)√(1-μ_d²) cos ϕ'
    a = μ₀ * μ_d
    b = sqrt(1.0 - μ₀^2) * sqrt(1.0 - μ_d^2)
    s = 0.0
    for k in 0:(N_phi-1)
        ϕ = 2π * k / N_phi
        cosΘ = a + b * cos(ϕ)
        # Clamp to avoid numerical |cosΘ| > 1 from roundoff
        cosΘ = clamp(cosΘ, -1.0, 1.0)
        s += exact_phase_function(c, cosΘ)
    end
    return s / N_phi
end

"""Closed-form Rayleigh azimuthal average (sanity-check shortcut)."""
function azimuthal_average_phase_rayleigh_analytic(μ_d, μ₀)
    return 0.75 * (1.0 + (μ₀ * μ_d)^2 + 0.5 * (1.0 - μ₀^2) * (1.0 - μ_d^2))
end

"""
Layer-effective azimuthally-averaged phase, scattering-weighted across
contributors, for downward direction μ_d at sun cosine μ₀.
"""
function layer_effective_azimuthal_phase(contributors, iz, μ_d, μ₀; N_phi=64)
    weight_sum = 0.0
    weighted_avg = 0.0
    for c in contributors
        ϖ = single_scattering_albedo(c)
        ϖ == 0 && continue
        w = optical_depth(c, iz) * ϖ
        weight_sum += w
        P_avg = azimuthal_average_phase(c, μ_d, μ₀; N_phi=N_phi)
        weighted_avg += w * P_avg
    end
    weight_sum == 0 && return 0.0
    return weighted_avg / weight_sum
end

"""
Per-layer τ-integral for path 3:
    ∫_{τ_top}^{τ_bot} exp(-τ/μ₀) exp(-(τ_total - τ)/μ_d) dτ

Numerically stable form: combine exponentials at each endpoint before
subtracting, to avoid overflow when μ_d is small.
"""
function path3_layer_τ_integral(τ_top, τ_bot, τ_total, μ₀, μ_d)
    b = 1.0/μ₀ - 1.0/μ_d
    # Endpoint values of integrand (without 1/b factor):
    f_top = exp(-τ_top/μ₀ - (τ_total - τ_top)/μ_d)
    f_bot = exp(-τ_bot/μ₀ - (τ_total - τ_bot)/μ_d)
    if abs(b) < 1e-10
        return 0.5 * (f_top + f_bot) * (τ_bot - τ_top)
    else
        return (f_top - f_bot) / b
    end
end

"""
Path 3: sun → atm → surface → sensor.

The downward irradiance on the surface from path-3 photons is

    F_path3 = ϖ I₀ ∫ exp(-τ_s/μ₀) [∫ P̄(μ_d, μ₀) exp(-(τ_total - τ_s)/μ_d) μ_d (1/2) dμ_d] dτ_s

After Lambertian reflection and upward attenuation:

    L_path3 = (albedo/π) · F_path3 · exp(-τ_total/μ_v)

Inner integration: Gauss-Legendre on (0, 1] for μ_d, with N_quad points.
Outer integration: per-layer closed form (path3_layer_τ_integral).
"""
function path3_atm_to_surface(contributors, μ₀, μ_v, albedo, I₀; N_quad=16, N_phi=64)
    τ_cum = cumulative_τ(contributors)
    τ_total = τ_cum[end]
    n_layers = length(τ_cum) - 1

    μ_nodes, μ_weights = gauss_legendre_01(N_quad)

    F_surface = 0.0
    for iz in 1:n_layers
        ϖ_eff_layer = layer_effective_ϖ(contributors, iz)
        ϖ_eff_layer == 0 && continue
        τ_top = τ_cum[iz]
        τ_bot = τ_cum[iz+1]

        layer_sum = 0.0
        for k in 1:N_quad
            μ_d = μ_nodes[k]
            w_d = μ_weights[k]
            P̄ = layer_effective_azimuthal_phase(contributors, iz, μ_d, μ₀; N_phi=N_phi)
            τ_int = path3_layer_τ_integral(τ_top, τ_bot, τ_total, μ₀, μ_d)
            # Factor 0.5 = 2π/4π from azimuthal normalization.
            # No μ_d: the 1/μ_d (RTE attenuation along slant) and μ_d (radiance →
            # irradiance conversion) cancel.
            layer_sum += w_d * P̄ * τ_int * 0.5
        end
        F_surface += ϖ_eff_layer * I₀ * layer_sum
    end

    L = (albedo / π) * F_surface * exp(-τ_total/μ_v)
    return L
end

"""
Path 4: sun → surface → atm → sensor.

Direct beam reaches surface: F_direct_to_surface = μ₀ I₀ exp(-τ_total/μ₀).
Lambertian reflection produces upward isotropic radiance:
    L_upward(μ_up, ϕ_up) = (albedo/π) F_direct_to_surface

This isotropic radiance experiences one atmospheric scatter. For each upward
direction (μ_u, ϕ_u), photons travel upward through (τ_total - τ_s) of optical
depth from surface to scattering depth τ_s, then scatter into the view direction
(μ_v, ϕ_v), then exit to TOA through τ_s of optical depth.

By azimuthal symmetry, the integral reduces to a 1D quadrature over μ_u,
weighted by the azimuthally-averaged phase function P̄(μ_u, μ_v).

For Lambertian + plane-parallel, this is equivalent to path 3 with μ₀ ↔ μ_v.
We compute it directly to verify reciprocity numerically.

L_path4 = (albedo/π) · μ₀ I₀ exp(-τ_total/μ₀) · ∫ ∫ ... [analogous structure]

The geometry: an upward-traveling photon scatters at depth τ_s within layer iz
and emerges along view direction. Source per unit dτ_s and per unit incoming
solid angle, integrated over the up-hemisphere of incoming directions:

    G_4(τ_s; μ_v) = ∫_{up hemi} (P̄(μ_u, μ_v)/2) · exp(-(τ_total-τ_s)/μ_u) μ_u dμ_u

(same structure as G in path 3, with μ_v replacing μ₀).

Outer τ_s integration: photon then exits up at angle μ_v through depth τ_s,
attenuated by exp(-τ_s/μ_v).

So:

    L_path4 = (albedo/π) · μ₀ I₀ exp(-τ_total/μ₀) · (ϖ/4π) ·
              ∫ exp(-τ_s/μ_v) · G_4(τ_s; μ_v) · dτ_s · 4π · ...

Wait, I need to redo this cleanly.

Photons leave the surface upward in direction (μ_u, ϕ_u) with radiance
(albedo/π)·μ₀ I₀ exp(-τ_total/μ₀). The irradiance at depth (τ_total - h) above
the surface in this direction is L · exp(-(τ_total - h)/μ_u). The scattering
source per unit volume from this beam into the view direction (μ_v, ϕ_v) is
(ϖ/4π) P(cos Θ_4) times the beam's radiance, where Θ_4 is the angle between
(μ_u, ϕ_u) and (μ_v, ϕ_v). The contribution to TOA radiance from scattering at
depth (τ_total - h) within layer iz, integrated over the up-hemisphere of
(μ_u, ϕ_u):

L_path4_per_layer(iz) = (albedo/π) · μ₀ I₀ exp(-τ_total/μ₀) · ϖ_eff(iz) ·
                        ∫_{τ_s ∈ layer} exp(-τ_s/μ_v) ·
                          [∫∫_up P̄(μ_u, μ_v) exp(-(τ_total - τ_s)/μ_u) μ_u dμ_u dϕ_u/(2π)] · (1/2) ·
                          dτ_s

By μ₀ ↔ μ_v symmetry with path 3, the inner integral is the same shape.
"""
function path4_surface_to_atm(contributors, μ₀, μ_v, albedo, I₀; N_quad=16, N_phi=64)
    τ_cum = cumulative_τ(contributors)
    τ_total = τ_cum[end]
    n_layers = length(τ_cum) - 1

    μ_nodes, μ_weights = gauss_legendre_01(N_quad)

    L = 0.0
    L_surface = (albedo / π) * μ₀ * I₀ * exp(-τ_total/μ₀)

    for iz in 1:n_layers
        ϖ_eff_layer = layer_effective_ϖ(contributors, iz)
        ϖ_eff_layer == 0 && continue
        τ_top = τ_cum[iz]
        τ_bot = τ_cum[iz+1]

        layer_sum = 0.0
        for k in 1:N_quad
            μ_u = μ_nodes[k]
            w_u = μ_weights[k]
            # Azimuthally averaged phase. By symmetry of the scattering-angle
            # geometry under (μ_u, μ_v) ↔ (μ_d, μ₀), we can re-use the same
            # function with arguments swapped.
            P̄ = layer_effective_azimuthal_phase(contributors, iz, μ_u, μ_v; N_phi=N_phi)

            b = 1.0/μ_v - 1.0/μ_u
            f_top = exp(-τ_top/μ_v - (τ_total - τ_top)/μ_u)
            f_bot = exp(-τ_bot/μ_v - (τ_total - τ_bot)/μ_u)
            if abs(b) < 1e-10
                τ_int = 0.5 * (f_top + f_bot) * (τ_bot - τ_top)
            else
                τ_int = (f_top - f_bot) / b
            end

            # No μ_u: source-into-view from radiance input has no μ_u weight.
            layer_sum += w_u * P̄ * τ_int * 0.5
        end
        # 1/μ_v from RTE solution along upward view path.
        L += L_surface * ϖ_eff_layer * layer_sum / μ_v
    end

    return L
end

# ============================================================================
# Top-level: total exact-SS radiance at TOA
# ============================================================================

"""
    exact_ss_toa(contributors, μ₀, μ_v, Δϕ, albedo, I₀; N_quad=16, N_phi=64)

Total TOA radiance from exact single scattering through a plane-parallel
atmosphere with Lambertian surface.

Returns NamedTuple with fields:
- `path1`: sun → atm → sensor
- `path2`: sun → surface → sensor
- `path3`: sun → atm → surface → sensor
- `path4`: sun → surface → atm → sensor
- `total`: sum of all four paths
"""
function exact_ss_toa(contributors, μ₀, μ_v, Δϕ, albedo, I₀;
                       N_quad=16, N_phi=64)
    τ_cum = cumulative_τ(contributors)
    τ_total = τ_cum[end]

    p1 = path1_atmospheric_ss(contributors, μ₀, μ_v, Δϕ, I₀)
    p2 = path2_surface_direct(τ_total, μ₀, μ_v, albedo, I₀)
    p3 = path3_atm_to_surface(contributors, μ₀, μ_v, albedo, I₀;
                              N_quad=N_quad, N_phi=N_phi)
    p4 = path4_surface_to_atm(contributors, μ₀, μ_v, albedo, I₀;
                              N_quad=N_quad, N_phi=N_phi)

    return (path1=p1, path2=p2, path3=p3, path4=p4,
            total = p1 + p2 + p3 + p4)
end

# ============================================================================
# Tests
# ============================================================================

function run_all_tests()
    @testset "ExactSSReference" begin

        @testset "scattering geometry" begin
            # μ₀=μ_v=0.5, Δϕ=0: cos Θ = -0.25 + 0.75 = 0.5
            @test scattering_angle_cosine(0.5, 0.5, 0.0) ≈ 0.5
            # μ₀=μ_v=0.5, Δϕ=π: cos Θ = -0.25 - 0.75 = -1
            @test scattering_angle_cosine(0.5, 0.5, π) ≈ -1.0  rtol=1e-12
            # μ₀=μ_v=1, any Δϕ: cos Θ = -1
            @test scattering_angle_cosine(1.0, 1.0, 0.0) ≈ -1.0  rtol=1e-12
            @test scattering_angle_cosine(1.0, 1.0, π) ≈ -1.0  rtol=1e-12
        end

        @testset "phase function normalization" begin
            ray = RayleighContributor([0.5])
            # ∫_{-1}^{1} P(μ) dμ / 2 = 1
            n_q = 10000
            s = 0.0
            for i in 1:n_q
                μ = -1.0 + 2.0 * (i - 0.5) / n_q
                s += exact_phase_function(ray, μ)
            end
            @test (2.0 / n_q) * s / 2.0 ≈ 1.0  rtol=1e-3

            for g in [0.0, 0.3, 0.6, -0.4]
                hg = HGAerosolContributor(g, 1.0, [0.1])
                s = 0.0
                for i in 1:n_q
                    μ = -1.0 + 2.0 * (i - 0.5) / n_q
                    s += exact_phase_function(hg, μ)
                end
                @test (2.0 / n_q) * s / 2.0 ≈ 1.0  rtol=1e-3
            end

            # HG with g=0 is isotropic
            hg_iso = HGAerosolContributor(0.0, 1.0, [0.1])
            for cosΘ in [-1.0, -0.3, 0.0, 0.5, 1.0]
                @test exact_phase_function(hg_iso, cosΘ) ≈ 1.0
            end
        end

        @testset "Rayleigh azimuthal average analytic vs numerical" begin
            ray = RayleighContributor([0.5])
            for μ_d in [0.1, 0.4, 0.8], μ₀ in [0.2, 0.5, 0.9]
                p_num = azimuthal_average_phase(ray, μ_d, μ₀; N_phi=128)
                p_ana = azimuthal_average_phase_rayleigh_analytic(μ_d, μ₀)
                @test p_num ≈ p_ana  rtol=1e-10
            end
        end

        @testset "Gauss-Legendre on (0,1] integrates polynomials exactly" begin
            for n in [4, 8, 16]
                x, w = gauss_legendre_01(n)
                # ∫₀¹ 1 dx = 1
                @test sum(w) ≈ 1.0  rtol=1e-12
                # ∫₀¹ x dx = 1/2
                @test sum(w .* x) ≈ 0.5  rtol=1e-12
                # ∫₀¹ x^k dx = 1/(k+1) for k up to 2n-1
                for k in 1:(2*n - 1)
                    @test sum(w .* x.^k) ≈ 1.0/(k+1)  rtol=1e-10
                end
            end
        end

        @testset "path1: single Rayleigh layer at zenith" begin
            # μ₀=μ_v=1, Δϕ=0 → cos Θ = -1, P_Ray = 1.5
            # Single layer τ=0.1, ϖ=1
            # Expected: (1/(4π·μ_v)) · ϖ · P · [1 - exp(-τ·a)]/a
            # NB: μ_v=1 here so this passes either way; not sufficient to test 1/μ_v.
            ray = RayleighContributor([0.1])
            r = path1_atmospheric_ss([ray], 1.0, 1.0, 0.0, 1.0)
            expected = (1.5 / (4π)) * (1.0 - exp(-0.2)) / 2.0
            @test r ≈ expected  rtol=1e-12
        end

        @testset "path1: off-zenith isotropic (catches the 1/μ_v factor)" begin
            # Isotropic phase (HG g=0), ϖ=1, P=1.
            # μ₀=0.5, μ_v=0.25 → 1/μ_v=4, so a missing 1/μ_v would underestimate
            # by a factor of 4. Strong test.
            iso = HGAerosolContributor(0.0, 1.0, [0.1])
            μ₀, μ_v, τ_layer = 0.5, 0.25, 0.1
            a = 1.0/μ₀ + 1.0/μ_v
            expected = (1.0 / (4π * μ_v)) * (1.0 - exp(-τ_layer * a)) / a
            r = path1_atmospheric_ss([iso], μ₀, μ_v, 0.0, 1.0)
            @test r ≈ expected  rtol=1e-12
        end

        @testset "path1: layer subdivision invariance" begin
            # Two-layer τ=(0.04, 0.06) ≡ one-layer τ=0.10 for path 1
            r1 = RayleighContributor([0.10])
            r2 = RayleighContributor([0.04, 0.06])
            for μ₀ in [0.3, 0.7], μ_v in [0.4, 0.9], Δϕ in [0.0, π/2, π]
                v1 = path1_atmospheric_ss([r1], μ₀, μ_v, Δϕ, 1.0)
                v2 = path1_atmospheric_ss([r2], μ₀, μ_v, Δϕ, 1.0)
                @test v1 ≈ v2  rtol=1e-12
            end
        end

        @testset "path1: Rayleigh + HG mixed phase consistency" begin
            ray = RayleighContributor([0.05])
            hg  = HGAerosolContributor(0.6, 0.9, [0.1])
            cosΘ = -0.7
            ϖ_eff = layer_effective_ϖ([ray, hg], 1)
            P_eff = layer_effective_phase([ray, hg], 1, cosΘ)
            ϖ_exp = (0.05*1.0 + 0.1*0.9) / 0.15
            P_ray_at = 0.75*(1 + cosΘ^2)
            g = 0.6
            P_hg_at = (1 - g^2) / (1 + g^2 - 2g*cosΘ)^1.5
            P_exp = (0.05*1.0*P_ray_at + 0.1*0.9*P_hg_at) / (0.05*1.0 + 0.1*0.9)
            @test ϖ_eff ≈ ϖ_exp  rtol=1e-12
            @test P_eff ≈ P_exp  rtol=1e-12
        end

        @testset "path2: closed form sanity" begin
            # No atmosphere: τ_total=0, just Lambertian.
            # L = μ₀ I₀ albedo / π
            empty_ray = RayleighContributor([0.0])
            r = exact_ss_toa([empty_ray], 0.5, 0.7, 0.3, 0.4, 1.0)
            @test r.path2 ≈ 0.5 * 1.0 * 0.4 / π  rtol=1e-12
            # Path 1 should be zero (no scattering).
            @test r.path1 ≈ 0.0
            # Paths 3 and 4 should be zero (no scattering).
            @test r.path3 ≈ 0.0  atol=1e-15
            @test r.path4 ≈ 0.0  atol=1e-15
        end

        @testset "black surface: paths 2,3,4 all zero" begin
            ray = RayleighContributor([0.2])
            r = exact_ss_toa([ray], 0.6, 0.4, π/3, 0.0, 1.0)
            @test r.path2 == 0.0
            @test r.path3 == 0.0
            @test r.path4 == 0.0
            @test r.total ≈ r.path1
        end

        @testset "absorption-only atmosphere: only path 2 nonzero" begin
            # No scattering anywhere; surface direct beam still reflects.
            abs = AbsorptionContributor([0.3])
            r = exact_ss_toa([abs], 0.5, 0.7, 0.0, 0.5, 1.0)
            @test r.path1 ≈ 0.0
            @test r.path3 ≈ 0.0  atol=1e-15
            @test r.path4 ≈ 0.0  atol=1e-15
            # Path 2: μ₀ I₀ albedo/π · exp(-0.3/μ₀) · exp(-0.3/μ_v)
            expected = 0.5 * 1.0 * 0.5 / π * exp(-0.3/0.5) * exp(-0.3/0.7)
            @test r.path2 ≈ expected  rtol=1e-12
        end

        @testset "path 1 BRDF reciprocity" begin
            # With the corrected 1/μ_v prefactor, path 1 satisfies:
            #     path1(μ₀, μ_v) / μ₀ = path1(μ_v, μ₀) / μ_v
            # which is BRDF reciprocity (the natural symmetry that emerges from
            # the RTE Green's function).
            ray = RayleighContributor([0.15])
            hg  = HGAerosolContributor(0.4, 0.95, [0.05])
            for (μ₀, μ_v) in [(0.6, 0.4), (0.3, 0.8), (0.1, 0.9)],
                Δϕ in [0.0, π/3, π]
                p_orig = path1_atmospheric_ss([ray, hg], μ₀, μ_v, Δϕ, 1.0)
                p_swap = path1_atmospheric_ss([ray, hg], μ_v, μ₀, Δϕ, 1.0)
                @test p_orig/μ₀ ≈ p_swap/μ_v  rtol=1e-12
            end
        end

        @testset "path 2 BRDF reciprocity" begin
            τ_total, albedo = 0.2, 0.3
            for (μ₀, μ_v) in [(0.6, 0.4), (0.3, 0.8), (0.5, 0.5)]
                p_orig = path2_surface_direct(τ_total, μ₀, μ_v, albedo, 1.0)
                p_swap = path2_surface_direct(τ_total, μ_v, μ₀, albedo, 1.0)
                @test p_orig/μ₀ ≈ p_swap/μ_v  rtol=1e-12
            end
        end

        @testset "paths 3+4 BRDF reciprocity (under μ₀ ↔ μ_v swap)" begin
            # Path 3 (sun → atm → surface → sensor) maps onto path 4
            # (sun → surface → atm → sensor) when (μ₀ ↔ μ_v) is swapped.
            # Specifically:
            #     path3(μ₀, μ_v) / μ₀ = path4(μ_v, μ₀) / μ_v
            #     path4(μ₀, μ_v) / μ₀ = path3(μ_v, μ₀) / μ_v
            ray = RayleighContributor([0.15])
            hg  = HGAerosolContributor(0.4, 0.95, [0.05])
            μ₀, μ_v = 0.6, 0.4
            albedo = 0.3
            r_orig = exact_ss_toa([ray, hg], μ₀, μ_v, 0.0, albedo, 1.0;
                                  N_quad=64, N_phi=128)
            r_swap = exact_ss_toa([ray, hg], μ_v, μ₀, 0.0, albedo, 1.0;
                                  N_quad=64, N_phi=128)
            @test r_orig.path3 / μ₀ ≈ r_swap.path4 / μ_v  rtol=1e-4
            @test r_orig.path4 / μ₀ ≈ r_swap.path3 / μ_v  rtol=1e-4
            # Full reciprocity: total reflectance is symmetric under (μ₀, μ_v) swap
            @test r_orig.total / μ₀ ≈ r_swap.total / μ_v  rtol=1e-4
        end

        @testset "Quadrature convergence" begin
            # Increasing N_quad and N_phi should converge to a fixed answer.
            ray = RayleighContributor([0.2])
            hg  = HGAerosolContributor(0.5, 0.9, [0.1])
            results = Float64[]
            for N_q in [8, 16, 32, 64]
                r = exact_ss_toa([ray, hg], 0.5, 0.6, π/4, 0.3, 1.0;
                                 N_quad=N_q, N_phi=128)
                push!(results, r.total)
            end
            d1 = abs(results[2] - results[1])
            d2 = abs(results[3] - results[2])
            d3 = abs(results[4] - results[3])
            # Differences should generally decrease (with small tolerance for
            # plateauing once quadrature error is below other noise).
            @test d2 ≤ d1 * 1.5
            @test d3 ≤ d2 * 1.5
            r_ref = exact_ss_toa([ray, hg], 0.5, 0.6, π/4, 0.3, 1.0;
                                 N_quad=128, N_phi=256)
            @test abs(results[end] - r_ref.total) / abs(r_ref.total) < 1e-4
        end

        @testset "Optically thin limit: linear in τ" begin
            # Use truly small τ (0.001) so exponentials are near 1 and we are
            # in the genuinely linear regime. At τ=0.1 the second-order term
            # is ~5% which spoils a strict-linearity test.
            base_τ = 0.001
            ray_unit = RayleighContributor([base_τ])
            r_unit = exact_ss_toa([ray_unit], 0.6, 0.5, 0.0, 0.0, 1.0)  # black surface
            base = r_unit.total
            for scale in [0.5, 0.25, 0.1]
                ray_scaled = RayleighContributor([base_τ * scale])
                r = exact_ss_toa([ray_scaled], 0.6, 0.5, 0.0, 0.0, 1.0)
                @test r.total / base ≈ scale  rtol=0.01
            end
        end

        @testset "Energy conservation: Rayleigh black surface, hemispheric R" begin
            # For thin Rayleigh + black surface, hemispheric reflectance from
            # SS only is well-defined and proportional to τ (in the linear regime).
            μ₀ = 0.7
            I₀ = 1.0
            n_q = 16
            μ_nodes, μ_weights = gauss_legendre_01(n_q)
            n_p = 16

            function hemispheric_R(τ)
                ray = RayleighContributor([τ])
                R = 0.0
                for k in 1:n_q
                    μ_v = μ_nodes[k]
                    I_avg = 0.0
                    for j in 0:(n_p - 1)
                        Δϕ = 2π * j / n_p
                        r = exact_ss_toa([ray], μ₀, μ_v, Δϕ, 0.0, I₀;
                                         N_quad=8, N_phi=32)
                        I_avg += r.total
                    end
                    I_avg /= n_p
                    R += μ_weights[k] * I_avg * μ_v * 2π
                end
                return R / (μ₀ * I₀)
            end

            R_full = hemispheric_R(0.01)
            R_half = hemispheric_R(0.005)
            @test 0.0 < R_full < 1.0
            @test R_half / R_full ≈ 0.5  rtol=0.02   # SS linear in τ at small τ
        end

        @testset "Reflectance limit: only surface contributes when τ=0" begin
            # No atmosphere: total should equal path 2 = μ₀ I₀ albedo/π.
            empty = RayleighContributor([0.0])
            μ₀, μ_v, albedo = 0.5, 0.3, 0.6
            r = exact_ss_toa([empty], μ₀, μ_v, 0.4, albedo, 1.0)
            @test r.total ≈ μ₀ * 1.0 * albedo / π  rtol=1e-12
        end

    end
end

end # module ExactSSReference

# Auto-run if invoked as script
if abspath(PROGRAM_FILE) == @__FILE__
    ExactSSReference.run_all_tests()
end
