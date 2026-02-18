"""
    RawAerosolJacobian{T}

Holds the **un-truncated** derivatives of single-aerosol optical properties with
respect to the aerosol sub-parameters `[τ_ref, nᵣ, nᵢ, rₘ, σ_g, p₀, σ_p]`.

This struct is the **AD boundary for aerosols**: everything upstream (Mie code,
aerosol profile code, or ForwardDiff) produces these raw derivatives. The
[`delta_m_truncation_lin`](@ref) function then maps them through the δ-M
truncation chain rule to produce the `CoreScatteringOpticalPropertiesLin` that
the RT kernel expects.

# Fields
- `τ̇_aer`: `∂τ_aer/∂x  [Nparams_aer × nSpec]` — optical depth derivatives
- `ω̃̇`:     `∂ω̃/∂x     [Nparams_mie × nSpec]` — SSA derivatives (Mie params only)
- `ḟᵗ`:     `∂fᵗ/∂x     [Nparams_mie × nSpec]` — truncation factor derivatives (Mie params only)
- `Ż⁺⁺`:    `∂Z⁺⁺/∂x   [Nparams_mie × nμ × nμ]` — phase matrix derivatives (Mie params only)
- `Ż⁻⁺`:    `∂Z⁻⁺/∂x   [Nparams_mie × nμ × nμ]` — phase matrix derivatives (Mie params only)
"""
struct RawAerosolJacobian{T}
    τ̇_aer::T
    ω̃̇::T
    ḟᵗ::T
    Ż⁺⁺::T
    Ż⁻⁺::T
end

"""
    delta_m_forward(τ_aer, ω̃, fᵗ, Z⁺⁺, Z⁻⁺)

Apply δ-M truncation (Nakajima & Tanaka 1988) to raw aerosol optical properties.

```math
\\tau_\\text{mod} = (1 - f^t \\tilde{\\omega}) \\, \\tau_\\text{aer}, \\quad
\\varpi_\\text{mod} = \\frac{(1 - f^t) \\, \\tilde{\\omega}}{1 - f^t \\tilde{\\omega}}
```

Phase matrices `Z⁺⁺, Z⁻⁺` pass through unchanged (truncation is already applied
to the Greek coefficients during Mie computation).

# Returns
`CoreScatteringOpticalProperties` with δ-M scaled `(τ_mod, ϖ_mod, Z⁺⁺, Z⁻⁺)`.
"""
function delta_m_forward(τ_aer, ω̃, fᵗ, Z⁺⁺, Z⁻⁺)
    τ_mod = (1 .- fᵗ * ω̃) .* τ_aer
    ϖ_mod = (1 .- fᵗ) .* ω̃ ./ (1 .- fᵗ * ω̃)
    return CoreScatteringOpticalProperties(τ_mod, ϖ_mod, Z⁺⁺, Z⁻⁺)
end

"""
    delta_m_truncation_lin(τ_aer, fᵗ, ω̃, τ̇_aer, ω̃̇, ḟᵗ, Ż⁺⁺, Ż⁻⁺, Z⁺⁺, Z⁻⁺,
                           n_aer_params, mie_range, profile_range, arr_type)

Apply the δ-M truncation chain rule analytically to map raw (un-truncated) aerosol
derivatives to δ-M modified derivatives.

The 4 raw variables `(τ_aer, ω̃, fᵗ, Z)` are mapped to 3 modified variables
`(τ_mod, ϖ_mod, Z)` via:

```math
\\dot{\\tau}_\\text{mod}^{(j)} = (1 - f^t \\tilde{\\omega}) \\, \\dot{\\tau}_\\text{aer}^{(j)}
  - \\tau_\\text{aer} \\left( f^t \\dot{\\tilde{\\omega}}^{(j)} 
  + \\tilde{\\omega} \\dot{f}^{t(j)} \\right)
```

```math
\\dot{\\varpi}_\\text{mod}^{(j)} = \\frac{
  \\dot{\\tilde{\\omega}}^{(j)} (1 - f^t)
  - \\dot{f}^{t(j)} \\tilde{\\omega} (1 - \\tilde{\\omega})
}{(1 - f^t \\tilde{\\omega})^2}
```

For profile-only parameters (`p₀, σ_p`), `ω̃̇ = ḟᵗ = Ż = 0`, so only the τ chain
contributes.

# Arguments
- `n_aer_params::Int`: Total aerosol sub-parameters (currently 7).
- `mie_range`: Indices within the 7-param block that have Mie sensitivity (currently `2:5`).
- `profile_range`: Indices that only affect τ via the vertical profile (currently `[1, 6, 7]`).
- `arr_type`: Array constructor (`Array` or `CuArray`).

# Returns
- `CoreScatteringOpticalProperties`: Forward δ-M scaled properties.
- `CoreScatteringOpticalPropertiesLin`: Linearized properties `[n_aer_params × nSpec]`.
"""
function delta_m_truncation_lin(τ_aer, fᵗ, ω̃, 
                                 τ̇_aer, ω̃̇, ḟᵗ, 
                                 Ż⁺⁺, Ż⁻⁺, Z⁺⁺, Z⁻⁺,
                                 n_aer_params::Int,
                                 mie_range, profile_range,
                                 arr_type)
    n = size(τ_aer, 1)

    fwd = delta_m_forward(τ_aer, ω̃, fᵗ, Z⁺⁺, Z⁻⁺)

    τ̇_mod = arr_type(zeros(eltype(τ_aer), n_aer_params, n))
    ϖ̇_mod = arr_type(zeros(eltype(τ_aer), n_aer_params, n))

    nμ = size(Z⁺⁺, 1)
    nμ2 = size(Z⁺⁺, 2)
    full_Ż⁺⁺ = arr_type(zeros(eltype(τ_aer), n_aer_params, nμ, nμ2))
    full_Ż⁻⁺ = arr_type(zeros(eltype(τ_aer), n_aer_params, nμ, nμ2))

    denom  = (1 .- fᵗ * ω̃)
    denom² = denom .^ 2

    # Profile-only params: ω̃̇=0, ḟᵗ=0, Ż=0
    for ip in profile_range
        τ̇_mod[ip, :] .= denom' .* τ̇_aer[ip, :]
    end

    # Mie-sensitive params: full δ-M chain rule
    nm = length(mie_range)
    τ̇_aer_mie = τ̇_aer[mie_range, :]

    tmp = fᵗ * ω̃̇ .+ ḟᵗ .* ω̃'
    τ̇_mod[mie_range, :] .= denom' .* τ̇_aer_mie .- tmp .* τ_aer'
    ϖ̇_mod[mie_range, :] .= (ω̃̇ .* (1 .- fᵗ) .- ḟᵗ .* (ω̃ .* (1 .- ω̃))') ./ denom²'

    full_Ż⁺⁺[mie_range, :, :] .= Ż⁺⁺
    full_Ż⁻⁺[mie_range, :, :] .= Ż⁻⁺

    lin = CoreScatteringOpticalPropertiesLin(τ̇_mod, ϖ̇_mod, full_Ż⁺⁺, full_Ż⁻⁺)
    return fwd, lin
end
