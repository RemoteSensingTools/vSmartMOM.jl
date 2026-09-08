# Shared assembly for analytic BRDFs

"""
    create_surface_layer!(brdf::AbstractSurfaceType, added_layer, SFI, m,
                          pol_type, quad_points, τ_sum, architecture)

Build an [`AddedLayer`](@ref) from any analytic BRDF that implements
`reflectance(brdf, stokes_index, μᵢ, μᵣ, Δϕ)`. RPV and Ross-Li specialize
that scalar interface; this shared method performs the Fourier integration and
operator assembly.
"""
function create_surface_layer!(brdf::AbstractSurfaceType,
                               added_layer::Union{AddedLayer,AddedLayerRS},
                               SFI,
                               m::Int,
                               pol_type,
                               quad_points,
                               τ_sum,
                               architecture;
                               F₀=nothing)
    (; qp_μ, qp_μN, wt_μN) = quad_points
    FT = eltype(qp_μN)
    Nquad = size(added_layer.r⁻⁺, 1) ÷ pol_type.n
    arr_type = array_type(architecture)
    T_surf = arr_type(Diagonal(ones(FT, pol_type.n * Nquad)))

    # The surface-layer convention carries an additional factor two at m=0.
    ρ = (m == 0 ? 2 : 1) * reflectance(brdf, pol_type, collect(qp_μ), m)
    R_surf = arr_type(ρ)
    if SFI
        _surface_source!(added_layer, R_surf, τ_sum, quad_points, pol_type,
                         FT, architecture; F₀)
    end
    _fill_surface_layer!(added_layer,
                         R_surf * Diagonal(qp_μN .* wt_μN), T_surf)
    return nothing
end

"""
    reflectance(brdf::AbstractSurfaceType, pol_type, μ, m)

Fourier moment `m` of an analytic BRDF on quadrature directions `μ`:

```math
R_{ij}^{(m)} = \\frac{f_m}{\\pi}\\int_0^\\pi
    \\rho(i, \\mu_i, \\mu_j, \\phi)\\cos(m\\phi)\\,d\\phi,
\\qquad f_0=1,\\;f_{m>0}=2.
```

Concrete BRDF types provide the five-argument scalar `reflectance` methods;
Julia dispatch selects them inside this common quadrature.
"""
function reflectance(brdf::AbstractSurfaceType, pol_type,
                     μ::AbstractArray{FT}, m::Int) where {FT}
    n_azimuth = 100
    array_constructor = array_type(architecture(μ))
    n = length(μ) * pol_type.n
    R_surf = array_constructor(zeros(FT, n, n))
    fourier_weight = m == 0 ? one(FT) : FT(2)

    for stokes_index in 1:pol_type.n
        integrand(ϕ) = reflectance.((brdf,), stokes_index, μ, μ', ϕ) * cos(m * ϕ)
        ϕ, weights = CanopyOptics.gauleg(n_azimuth, zero(FT), FT(π))
        samples = integrand.(ϕ)
        block = fill!(similar(R_surf[stokes_index:pol_type.n:end,
                                    stokes_index:pol_type.n:end]), zero(FT))
        for i in eachindex(samples)
            block += weights[i] * samples[i]
        end
        R_surf[stokes_index:pol_type.n:end,
               stokes_index:pol_type.n:end] .= block / π
    end
    return fourier_weight * R_surf
end
