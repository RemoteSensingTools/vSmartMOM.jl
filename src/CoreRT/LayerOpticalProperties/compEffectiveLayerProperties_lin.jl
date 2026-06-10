#=
AD Architecture Boundary
========================
This file is the interface between the "physical parameter → optical property" chain rule
and the "optical property → TOA radiance" adding-doubling propagation.

Physical parameters (VMR, aerosol τ_ref/n_r/n_i/size/profile, surface albedo) produce
derivatives ∂(τ,ϖ,Z)/∂x_phys via Mie theory (ForwardDiff), Beer-Lambert law, and
mixing rules (quotient rule). These derivatives are computed HERE and stored in
`CoreScatteringOpticalPropertiesLin`.

The RT kernel (elemental!, doubling!, interaction!) then only needs to propagate
w.r.t. the 3 core variables (τ, ϖ, Z) using the chain rule:
    ∂R/∂x = ∂R/∂τ · ∂τ/∂x + ∂R/∂ϖ · ∂ϖ/∂x + ∂R/∂Z · ∂Z/∂x
This clean separation makes it straightforward to add new physical parameters
without touching the RT propagation code.
=#

"""
    constructCoreOpticalProperties(RS_type, iBand, m, model, lin_model)

Construct the combined effective layer optical properties and their linearized derivatives
for all atmospheric layers, merging Rayleigh scattering, aerosol scattering (with δ-M
truncation), and trace gas absorption.

For each atmospheric layer, this function:
1. Creates Rayleigh scattering properties (τ_rayl, ϖ=1, Z_rayl).
2. For each aerosol type, computes δ-M scaled optical properties via [`createAero`](@ref),
   combining Mie scattering with the aerosol vertical profile.
3. Adds Rayleigh and aerosol properties using the `+` operator on
   `UmbrellaCoreScatteringOpticalProperties`, which correctly propagates derivatives through
   the mixing formulas for combined ϖ and Z.
4. Adds gas absorption via `UmbrellaCoreAbsorptionOpticalProperties`, updating ϖ and τ.

The output derivative dimension has `Nparams = 7×NAer + NGas` entries per layer
(surface parameters are handled separately in the interaction step).

# Returns
- `layer_opt`: Array of `CoreScatteringOpticalProperties` (one per layer).
- `layer_opt_lin`: Array of `CoreScatteringOpticalPropertiesLin` (derivatives, one per layer).
- `fscat_opt`: Rayleigh scattering fraction per layer (for inelastic scattering weight).
"""
function constructCoreOpticalProperties(RS_type, iBand, m, model, lin_model) #where {FT<:Real}
    (; τ_rayl, τ_aer, τ_abs, aerosol_optics, 
            greek_rayleigh, greek_cabannes, ϖ_Cabannes) = model
    (; τ̇_aer, τ̇_abs, lin_aerosol_optics) = lin_model
    @assert all(iBand .≤ length(τ_rayl)) "iBand exceeded number of bands"
    FT = eltype(τ_rayl[1])
    
    arr_type = CoreRT.array_type(model)

    pol_type = CoreRT.polarization_type(model)
    # Do this in CPU space only first:
    
    # Quadrature points:
    μ = collect(model.quad_points.qp_μ )
    # Number of Aerosols:
    nAero = size(τ_aer[iBand[1]],1)
    nZ    = size(τ_rayl[1],2)
    # Rayleigh Z matrix:
    
    band_layer_props     = [];
    band_layer_props_lin = [];
    band_fScattRayleigh  = [];
    for iB in iBand
        if (typeof(RS_type)<:noRS) #!(typeof(RS_type)<:RRS)
            Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, μ, 
                                                            greek_rayleigh[iB], m, 
                                                            arr_type = arr_type);
        else
            Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, μ, 
                                                            greek_cabannes[iB], m, 
                                                            arr_type = arr_type);
        end

        rayl = [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]), 1.0,
                (Rayl𝐙⁺⁺), (Rayl𝐙⁻⁺)) for i=1:nZ]
        # Initiate combined properties with rayleigh
        combrella = [UmbrellaCoreScatteringOpticalProperties(rayl[i],nothing) for i=1:nZ]
        # Loop over all aerosol types:
        for iaer=1:nAero
            # Precomute Z matrices per type (constant per layer)
            AerZ⁺⁺, AerZ⁻⁺, AerŻ⁺⁺, AerŻ⁻⁺ = Scattering.compute_Z_moments(
                                pol_type, μ, 
                                aerosol_optics[iB][iaer].greek_coefs, 
                                lin_aerosol_optics[iB][iaer].lin_greek_coefs, 
                                m, arr_type=arr_type)
            # Generate Core optical properties for Aerosols iaer
            aer = []
            lin_aer = []
            for iz=1:nZ
                t_aer, t_lin_aer =  createAero(
                            arr_type(τ_aer[iB][iaer,:,iz]), 
                            aerosol_optics[iB][iaer], 
                            arr_type(AerZ⁺⁺), arr_type(AerZ⁻⁺), 
                            arr_type(τ̇_aer[iB][iaer,:,:,iz]), 
                            lin_aerosol_optics[iB][iaer], 
                            arr_type(AerŻ⁺⁺), arr_type(AerŻ⁻⁺), 
                            arr_type) 
                push!(aer, t_aer)
                push!(lin_aer, t_lin_aer)
            end 
            aer_combrella = [UmbrellaCoreScatteringOpticalProperties(aer[i],lin_aer[i]) for i=1:nZ]
            # Mix with previous Core Optical Properties
            combrella = [combrella[i]+aer_combrella[i] for i=1:nZ]
        end


        # fScattRayleigh:
        # Assume ϖ of 1 for Rayleight here:
        # Create Core Optical Properties merged with trace gas absorptions:
        gas = [CoreAbsorptionOpticalProperties(arr_type(τ_abs[iB][:,iz])) for iz=1:nZ]
        lin_gas = [CoreAbsorptionOpticalPropertiesLin(arr_type(collect(τ̇_abs[iB][:,:,iz]'))) for iz=1:nZ] 
        gas_combrella = [UmbrellaCoreAbsorptionOpticalProperties(gas[iz],lin_gas[iz]) for iz=1:nZ]   
        tmp = combrella .+ gas_combrella            
        combrella = tmp
        fScattRayleigh = [collect(rayl[iz].τ  ./ combrella[iz].fwd.τ) for iz=1:nZ] 
        combo =     [combrella[iz].fwd for iz=1:nZ]
        combo_lin = [combrella[iz].lin for iz=1:nZ]

        push!(band_layer_props, combo)
        push!(band_layer_props_lin, combo_lin)
        push!(band_fScattRayleigh,fScattRayleigh)
    end
    layer_opt     = [prod([band_layer_props[i][iz] for i=1:length(iBand)]) for iz=1:nZ]
    layer_opt_lin = [prod([band_layer_props_lin[i][iz] for i=1:length(iBand)]) for iz=1:nZ]
    fscat_opt     = [[band_fScattRayleigh[i][iz] for i=1:length(iBand)] for iz=1:nZ]
    return layer_opt, layer_opt_lin, fscat_opt
end
 

"""
    createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺, τ̇Aer, lin_aerosol_optics, AerŻ⁺⁺, AerŻ⁻⁺, arr_type)

!!! note "Legacy function"
    Consider using [`delta_m_truncation_lin`](@ref) directly for new code.
    It provides the same δ-M chain rule in a cleaner, composable form
    with explicit `RawAerosolJacobian` as the AD boundary.

Compute **δ-M scaled** aerosol optical properties and their derivatives for one aerosol
type in one atmospheric layer.

Applies the δ-M truncation correction (Nakajima & Tanaka, 1988; Sanghavi & Stephens, 2013)
to the aerosol optical depth and single-scattering albedo, and computes derivatives with
respect to 7 aerosol sub-parameters: `[τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]`.

# δ-M scaling formulas

```math
\\tau_\\text{mod} = (1 - f^t \\tilde{\\omega}) \\cdot \\tau_\\text{aer}
```
```math
\\varpi_\\text{mod} = \\frac{(1 - f^t) \\tilde{\\omega}}{1 - f^t \\tilde{\\omega}}
```

# Derivative chain rule

For each physical parameter ``p_j``:
```math
\\frac{\\partial \\tau_\\text{mod}}{\\partial p_j} = 
  (1 - f^t \\tilde{\\omega}) \\frac{\\partial \\tau_\\text{aer}}{\\partial p_j}
  - \\tau_\\text{aer} \\left(f^t \\frac{\\partial \\tilde{\\omega}}{\\partial p_j} + 
    \\tilde{\\omega} \\frac{\\partial f^t}{\\partial p_j}\\right)
```

For `τ_ref, p₀, σp`: only the `τ_aer` chain contributes (Mie properties are independent).

!!! note "Bug 18 fix"
    Corrected index mapping: `τ̇Aer[k,:]` now correctly maps to parameter `k`
    (was previously off-by-one, mixing τ_ref with Mie parameter derivatives).

# Returns
- `CoreScatteringOpticalProperties`: Forward δ-M scaled properties.
- `CoreScatteringOpticalPropertiesLin`: Linearized properties (7 sub-params).
"""
# ----------------------------------------------------------------------
# Helper: normalise the Mie-parameter linearisation arrays into a
# uniform `(n_spec, 4)` block for the v0.6+ Jacobian assembly path.
#
# `lin_aerosol_optics.ω̃̇` and `.ḟᵗ` carry derivatives w.r.t. the four Mie
# parameters (n_r, n_i, r_m, sigma_r). Their concrete shape depends on
# the Fourier truncation chosen at parameters_from_yaml time:
#   ()                 -- scalar (NoTruncation: f^t = 0 everywhere)
#   (nparams,)         -- spectrally flat (legacy delta-BGE pre-v0.7)
#   (nparams, 1)       -- degenerate matrix form
#   (nparams, n_spec)  -- fully resolved (v0.7 auto + Mie Dual chain)
# `createAero` then writes columns 2:5 of the 7-parameter layout
# [tau_ref, n_r, n_i, r_m, sigma_r, p_0, sigma_p] uniformly.
function _lift_mie_param_to_n_x_4(x, n, arr_type)
    eltyp = eltype(x)
    if x isa Number
        return arr_type(fill(eltyp(x), n, 4))
    elseif ndims(x) == 1
        # (nparams,) -> broadcast to (n, nparams), pad/truncate to 4 cols
        x4 = length(x) == 4 ? collect(x) : _pad4(collect(x), eltyp)
        return arr_type(repeat(reshape(x4, 1, 4), n, 1))
    elseif ndims(x) == 2
        sz = size(x)
        # Accept (nparams, 1) (degenerate) -- broadcast across spectral.
        if sz[2] == 1
            x4 = sz[1] == 4 ? collect(x[:, 1]) : _pad4(collect(x[:, 1]), eltyp)
            return arr_type(repeat(reshape(x4, 1, 4), n, 1))
        end
        # (nparams, n_spec) -- transpose to (n_spec, nparams), pad to 4 cols.
        if sz[1] in (1, 2, 3, 4) && sz[2] == n
            xt = collect(transpose(x))   # (n, nparams)
            return sz[1] == 4 ? arr_type(xt) : arr_type(_pad4_cols(xt, eltyp))
        end
        # (n_spec, nparams) already in the assembly orientation.
        if sz[1] == n && sz[2] in (1, 2, 3, 4)
            return sz[2] == 4 ? arr_type(collect(x)) : arr_type(_pad4_cols(collect(x), eltyp))
        end
    end
    error("createAero: unexpected shape ", size(x), " for Mie parameter linearisation array; expected scalar, (nparams,), or 2-D with one dim equal to n_spec=", n)
end

_pad4(v::AbstractVector, T) = (out = zeros(T, 4); out[1:length(v)] .= v; out)
_pad4_cols(M::AbstractMatrix, T) = (out = zeros(T, size(M, 1), 4); out[:, 1:size(M, 2)] .= M; out)

function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺, 
                    τ̇Aer, lin_aerosol_optics, AerŻ⁺⁺, AerŻ⁻⁺,
                    arr_type)

    (; fᵗ, ω̃) = aerosol_optics
    (; ḟᵗ, ω̃̇) = lin_aerosol_optics

    n  = size(τAer,1)
    τ̇Aer = arr_type(collect(τ̇Aer'))  # (7, n) → (n, 7) — transpose upstream input
    #fᵗ = arr_type(fᵗ)
    #ω̃  = arr_type(ω̃)
    #ḟᵗ = arr_type(ḟᵗ)
    #ω̃̇  = arr_type(ω̃̇)

    sz = size(AerŻ⁺⁺)
    tmpŻ⁺⁺ = AerŻ⁺⁺
    tmpŻ⁻⁺ = AerŻ⁻⁺
    AerŻ⁺⁺ = arr_type(zeros(size(AerZ⁺⁺,1), size(AerZ⁺⁺,2), 7))#, n)
    AerŻ⁻⁺ = arr_type(zeros(size(AerZ⁻⁺,1), size(AerZ⁻⁺,2), 7))#, n)
    AerŻ⁺⁺[:,:,2:5] .= permutedims(tmpŻ⁺⁺, (2, 3, 1))
    AerŻ⁻⁺[:,:,2:5] .= permutedims(tmpŻ⁻⁺, (2, 3, 1))

    # Ensure arrays are in the right memory space (CPU or GPU)
    ω̃  = (ω̃ isa Number) ? arr_type(fill(ω̃,n)) : arr_type(ω̃)
    #fᵗ  = (fᵗ isa Number) ? arr_type(fill(fᵗ,n)) : arr_type(fᵗ)
    
    ω̃̇ = arr_type(ω̃̇)
    ḟᵗ = arr_type(ḟᵗ)
    
    # ω̃̇ and ḟᵗ from Mie are derivatives w.r.t. the 4 Mie
    # parameters (n_r, n_i, r_m, sigma_r). They arrive in any of:
    #   ()                 -- scalar (e.g. NoTruncation case)
    #   (nparams,)         -- spectrally flat (legacy delta-BGE pre-v0.7)
    #   (nparams, 1)       -- degenerate matrix form
    #   (nparams, n_spec)  -- fully resolved (v0.7 Mie Dual chain)
    # `_lift_mie_param_to_n_x_4` normalises all four to (n_spec, 4) so
    # we can write columns 2:5 of the 7-parameter
    # [tau_ref, n_r, n_i, r_m, sigma_r, p0, sigma_p] layout uniformly.
    # Replaces an `eachindex` loop on ḟᵗ that overflowed when the
    # matrix form contained n_spec*4 > 6 entries (BoundsError on the
    # docs build's `_jacobian_aod_plot` after the v0.7 YAML migration).
    ω̃̇_block = _lift_mie_param_to_n_x_4(ω̃̇, n, arr_type)
    ḟᵗ_block = _lift_mie_param_to_n_x_4(ḟᵗ, n, arr_type)
    ω̃̇ = arr_type(zeros(n, 7));  ω̃̇[:, 2:5] .= ω̃̇_block
    ḟᵗ = arr_type(zeros(n, 7));  ḟᵗ[:, 2:5] .= ḟᵗ_block

    # Forward modified properties
    τ_mod = (1 .- fᵗ * ω̃) .* τAer
    ϖ_mod = (1 .- fᵗ) .* ω̃ ./ (1 .- fᵗ * ω̃)

    # Allocate linearized outputs
    τ̇_mod = arr_type(zeros(n,7)) #similar(τ̇Aer)
    ϖ̇_mod = arr_type(zeros(n,7)) #similar(ω̃̇)

    # Bug 18 fix: derivatives w.r.t. 7 aerosol sub-params [τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]
    τ̇_mod[:,1] .= (1 .- fᵗ * ω̃) .* τ̇Aer[:,1]
    ϖ̇_mod[:,1] .= 0.0
    # Vectorized form over iparam dimension
    # Dimensions: iparam × spectral
    tmp = fᵗ .* ω̃̇[:,2:5] .+ ω̃ .* ḟᵗ[:,2:5]  # (nSpec, iparam)
    τ̇_mod[:,2:5] .= (1 .- fᵗ * ω̃) .* τ̇Aer[:,2:5] .- tmp .* τAer  # (nSpec, iparam)
    ϖ̇_mod[:,2:5] .= (ω̃̇[:,2:5] .* (1 - fᵗ) .- ḟᵗ[:,2:5] .* (ω̃ .* (1 .- ω̃))) ./ (1 .- fᵗ * ω̃).^2

    τ̇_mod[:,6:7] .= (1 .- fᵗ * ω̃) .* τ̇Aer[:,6:7]
    ϖ̇_mod[:,6:7] .= 0.0
    
    return CoreScatteringOpticalProperties(τ_mod, ϖ_mod, AerZ⁺⁺, AerZ⁻⁺), 
        CoreScatteringOpticalPropertiesLin(τ̇_mod, ϖ̇_mod, AerŻ⁺⁺, AerŻ⁻⁺)
end


# Extract scattering definitions and integrated absorptions for the source function!
"""
    extractEffectiveProps(lods, lods_lin)

Extract effective scattering properties from combined layer optical property arrays.

Computes the cumulative optical depth sums and determines scattering interfaces for
each layer. The scattering interface type (`ScatteringInterface_00`, `01`, `10`, `11`)
controls which adding method variant is dispatched during the interaction step.

# Returns
- `scattering_interfaces_all`: Array of `ScatteringInterface` types per layer.
- `τ_sum_all`: Cumulative optical depth `[nSpec × (Nz+1)]`.
- `τ̇_sum_all`: Cumulative τ derivative `[nSpec × Nparams × (Nz+1)]`.
"""
function extractEffectiveProps(
                lods::Array,#{CoreScatteringOpticalProperties{FT},1}
                lods_lin::Array) #where FT

    FT    = eltype(lods[1].τ)
    nSpec = length(lods[1].τ)
    nZ    = length(lods)
    nParams = size(lods_lin[1].τ̇, 2)
    # First the Scattering Interfaces:
    scattering_interface = ScatteringInterface_00()
    scattering_interfaces_all = []
    τ_sum_all = similar(lods[1].τ,(nSpec,nZ+1))
    τ_sum_all[:,1] .= 0
    τ̇_sum_all = similar(lods_lin[1].τ̇,(nSpec,nParams,nZ+1))
    τ̇_sum_all[:,:,1] .= 0
    for iz =1:nZ
        # Need to check max entries in Z matrices here as well later!
        scatter = maximum(lods[iz].τ .* lods[iz].ϖ) > 2eps(FT)
        scattering_interface = get_scattering_interface(scattering_interface, scatter, iz)
        push!(scattering_interfaces_all, scattering_interface)
        @views τ_sum_all[:,iz+1] = τ_sum_all[:,iz] + lods[iz].τ 
        for ip=1:nParams
            τ̇_sum_all[:,ip,iz+1] = τ̇_sum_all[:,ip,iz] + lods_lin[iz].τ̇[:,ip] 
        end
    end
    return scattering_interfaces_all, τ_sum_all, τ̇_sum_all
end

"""
    expandOpticalProperties(in, in_lin, arr_type; expand_Z=false)

Convert optical properties and their linearized derivatives to the target
device array type.  When `expand_Z=false` (default for RT-kernel callers)
singleton-Z arrays (`size(Z,3)==1`) are **not** replicated to `nSpec` copies —
the elemental and fused linearized kernels already branch on `size(Z,3)>1` and
use index `n2=1` when the phase matrix is spectrally flat, so no replication
is needed.  Set `expand_Z=true` only when the result will be concatenated
along the spectral dimension (e.g. the `Base.:*` multi-layer combiner).
"""
function expandOpticalProperties(in::CoreScatteringOpticalProperties,
                                  in_lin::CoreScatteringOpticalPropertiesLin,
                                  arr_type;
                                  expand_Z::Bool = false)
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = in 
    (; τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺) = in_lin 
    @assert length(τ) == length(ϖ) "τ and ϖ sizes need to match"
    @assert length(τ̇) == length(ϖ̇) "τ̇ and ϖ̇ sizes need to match"

    if size(Z⁺⁺, 3) == 1
        if expand_Z
            Z⁺⁺ = _repeat(Z⁺⁺, 1, 1, length(τ))
            Z⁻⁺ = _repeat(Z⁻⁺, 1, 1, length(τ))
            Ż⁺⁺ = _repeat(reshape(Ż⁺⁺, size(Ż⁺⁺,1), size(Ż⁺⁺,2), 1, size(Ż⁺⁺,3)), 1, 1, length(τ), 1)
            Ż⁻⁺ = _repeat(reshape(Ż⁻⁺, size(Ż⁻⁺,1), size(Ż⁻⁺,2), 1, size(Ż⁻⁺,3)), 1, 1, length(τ), 1)
            return CoreScatteringOpticalProperties(
                        _to_device(arr_type, τ), _to_device(arr_type, ϖ),
                        _to_device(arr_type, Z⁺⁺), _to_device(arr_type, Z⁻⁺)),
                   CoreScatteringOpticalPropertiesLin(
                        _to_device(arr_type, τ̇), _to_device(arr_type, ϖ̇),
                        _to_device(arr_type, Ż⁺⁺), _to_device(arr_type, Ż⁻⁺))
        else
            # Fast path: fused linearized kernels branch on size(Z,3) and
            # size(Ż,3), using n2=1 / n2_lin=1 for singleton spectral dim.
            # Ensure 3-D / 4-D shapes so Z[i,j,1] and Ż[i,j,1,p] are valid.
            Z3⁺⁺ = _to_device(arr_type, _ensure_3d(Z⁺⁺))
            Z3⁻⁺ = _to_device(arr_type, _ensure_3d(Z⁻⁺))
            Ż4⁺⁺ = _to_device(arr_type, ndims(Ż⁺⁺) == 3 ? reshape(Ż⁺⁺, size(Ż⁺⁺,1), size(Ż⁺⁺,2), 1, size(Ż⁺⁺,3)) : Ż⁺⁺)
            Ż4⁻⁺ = _to_device(arr_type, ndims(Ż⁻⁺) == 3 ? reshape(Ż⁻⁺, size(Ż⁻⁺,1), size(Ż⁻⁺,2), 1, size(Ż⁻⁺,3)) : Ż⁻⁺)
            return CoreScatteringOpticalProperties(
                        _to_device(arr_type, τ), _to_device(arr_type, ϖ), Z3⁺⁺, Z3⁻⁺),
                   CoreScatteringOpticalPropertiesLin(
                        _to_device(arr_type, τ̇), _to_device(arr_type, ϖ̇), Ż4⁺⁺, Ż4⁻⁺)
        end
    else
        @assert size(Z⁺⁺, 3) == length(τ) "Z and τ dimensions need to match"
        return CoreScatteringOpticalProperties(
                    _to_device(arr_type, τ), _to_device(arr_type, ϖ),
                    _to_device(arr_type, Z⁺⁺), _to_device(arr_type, Z⁻⁺)),
               CoreScatteringOpticalPropertiesLin(
                    _to_device(arr_type, τ̇), _to_device(arr_type, ϖ̇),
                    _to_device(arr_type, Ż⁺⁺), _to_device(arr_type, Ż⁻⁺))
    end
end
