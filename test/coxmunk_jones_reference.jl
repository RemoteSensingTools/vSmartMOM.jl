#=
Independent electric-field (Jones) references for polarized reflection,
used by test_coxmunk_polarization.jl. Nothing here calls the Cox–Munk
implementation: Stokes matrices are formed from 3-D field maps and explicit
basis vectors, M_ij = ½ tr(σᵢ A σⱼ A†).

Geometry (same as the atmospheric phase matrix): incident propagation
kᵢ = (sinθᵢ, 0, −μᵢ), reflected/scattered kᵣ = (sinθᵣ cosΔϕ, sinθᵣ sinΔϕ, μᵣ).
Stokes basis (e₁, e₂) = (φ̂, θ̂) of each beam with
Q = |E₁|² − |E₂|², U = 2Re(E₁E₂*), V = 2Im(E₁E₂*); the first test in
test_coxmunk_polarization.jl establishes that this is vSmartMOM's convention
by reproducing the code's own Rayleigh phase matrix.
=#
module JonesReference

using LinearAlgebra

kin(μ)      = [sqrt(1 - μ^2), 0.0, -μ]
kout(μ, Δϕ) = [sqrt(1 - μ^2) * cos(Δϕ), sqrt(1 - μ^2) * sin(Δϕ), μ]

"Meridian basis (φ̂, θ̂) of a beam with azimuth ϕ and propagation k."
function basis(k, ϕ)
    s = hypot(k[1], k[2])
    θ̂ = [k[3] * cos(ϕ), k[3] * sin(ϕ), -s]
    φ̂ = [-sin(ϕ), cos(ϕ), 0.0]
    return φ̂, θ̂
end

const σ = (ComplexF64[1 0; 0 1], ComplexF64[1 0; 0 -1],
           ComplexF64[0 1; 1 0], ComplexF64[0 im; -im 0])

"4×4 Stokes matrix of a 3-D linear field map J for kᵢ(μᵢ) → kᵣ(μᵣ, Δϕ)."
function stokes_matrix(J, μᵢ, μᵣ, Δϕ)
    bᵢ = basis(kin(μᵢ), 0.0)
    bᵣ = basis(kout(μᵣ, Δϕ), Δϕ)
    A = ComplexF64[dot(bᵣ[j], J(bᵢ[k])) for j in 1:2, k in 1:2]
    return [real(tr(σ[i] * A * σ[j] * A')) / 2 for i in 1:4, j in 1:4]
end

"Electric-dipole (Rayleigh, ρ = 0) scattering, phase function normalized to 4π."
function rayleigh(μᵢ, μᵣ, Δϕ)
    kᵣ = kout(μᵣ, Δϕ)
    return 1.5 .* stokes_matrix(E -> E - kᵣ * dot(kᵣ, E), μᵢ, μᵣ, Δϕ)
end

"""
Cox–Munk glint BRDF (sr⁻¹) for per-axis slope variance σ², complex index n:
P(zx, zy)/(4 μᵢ μᵣ cos⁴β) times the facet's Fresnel Stokes matrix. The facet
map uses r_s, r_p in the basis (ŝ, p̂ = k × ŝ), for which a perfect conductor
(r_s = −1, r_p = +1) gives E ↦ −E + 2(E·n̂)n̂.
"""
function coxmunk(μᵢ, μᵣ, Δϕ, n, σ²)
    kᵢ, kᵣ = kin(μᵢ), kout(μᵣ, Δϕ)
    n̂ = normalize(kᵣ - kᵢ)
    cosθ = dot(kᵣ, n̂)
    cost = sqrt(complex(1 - (1 - cosθ^2) / n^2))
    rs = (cosθ - n * cost) / (cosθ + n * cost)
    rp = (n * cosθ - cost) / (n * cosθ + cost)
    s = cross(kᵢ, kᵣ)
    ŝ = norm(s) > 1e-14 ? normalize(s) : [0.0, 1.0, 0.0]
    p̂ᵢ, p̂ᵣ = cross(kᵢ, ŝ), cross(kᵣ, ŝ)
    J(E) = rs * dot(ŝ, E) * ŝ + rp * dot(p̂ᵢ, E) * p̂ᵣ
    cosβ = n̂[3]
    P = exp(-(1 - cosβ^2) / cosβ^2 / (2σ²)) / (2π * σ²)
    return P / (4 * μᵢ * μᵣ * cosβ^4) .* stokes_matrix(J, μᵢ, μᵣ, Δϕ)
end

end # module
