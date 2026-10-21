#=
 
This file contains RT elemental-related functions
 
=#
"""
    elemental!(pol_type, SFI, τ_sum, τ̇_sum, dτ, F₀, computed_layer_properties,
               computed_layer_properties_lin, m, ndoubl, scatter, quad_points,
               added_layer, added_layer_lin, architecture)

Tangent-linear partner of the forward [`elemental!`](@ref) kernel: builds
the elemental-layer reflection, transmission, and source matrices **and
their derivatives** with respect to the three core layer variables
``(\\tau, \\varpi, \\mathbf{Z})`` in a single sweep.

# Forward formulas

The forward kernel uses the **exact finite-δ single-scatter** formulas of
Fell (1997, FU Berlin PhD thesis, Eqs. 1.52–1.56), restated as Eqs. (10)–(11)
of Sanghavi & Frankenberg (2023, JQSRT 311:108791) — **not** the
infinitesimal-δ linear limit of Sanghavi et al. (2014, JQSRT 133:412–433),
Eqs. (19)–(20). See [`get_elem_rt!`](@ref) and the Concepts/04 § Elemental
page for the side-by-side comparison.

```math
\\mathbf{r}^{-+}(\\mu_i, \\mu_j) = \\varpi \\, \\mathbf{Z}^{-+}(\\mu_i, \\mu_j)
  \\frac{\\mu_j}{\\mu_i + \\mu_j} \\left(1 - e^{-d\\tau(1/\\mu_i + 1/\\mu_j)}\\right) w_j
```

```math
\\mathbf{t}^{++}(\\mu_i, \\mu_j) = \\delta_{ij} e^{-d\\tau/\\mu_i} +
  \\varpi \\, \\mathbf{Z}^{++}(\\mu_i, \\mu_j)
  \\frac{\\mu_j}{\\mu_i - \\mu_j} \\left(e^{-d\\tau/\\mu_i} - e^{-d\\tau/\\mu_j}\\right) w_j
```

# Linearization (Sanghavi 2014 App. C)

Uses the product/chain-rule framework of Sanghavi 2014 App. C on the
finite-δ formulas above. Eqs. (C.8)–(C.10) give the infinitesimal-layer
and thermal-source versions; the finite-δ and direct-solar expressions
here follow by differentiating SF2023-II (10)–(11). For each matrix
element, three local partials are stored along the last axis:

- ``\\dot{\\mathbf{M}}[\\,..,1]``: ``\\partial \\mathbf{M}/\\partial(d\\tau)`` — optical-depth derivative
- ``\\dot{\\mathbf{M}}[\\,..,2]``: ``\\partial \\mathbf{M}/\\partial\\varpi`` — single-scatter albedo derivative
- ``\\dot{\\mathbf{M}}[\\,..,3]``: ``\\partial \\mathbf{M}/\\partial\\mathbf{Z}`` — phase-matrix derivative

These three "core" derivatives feed the chain rule that is fused directly
into `get_elem_rt_fused!` (non-SFI path) and `get_elem_rt_SFI_fused!`
(direct-solar path).  These fused kernels expand the core derivatives into
the supplied tangent directions (physical parameters or a compact local
optical basis) using the boundary inputs
`CoreScatteringOpticalPropertiesLin = (\\dot{\\tau}, \\dot{\\varpi}, \\dot{\\mathbf{Z}}^{++}, \\dot{\\mathbf{Z}}^{-+})`.

When `SFI=true`, source vectors ``\\mathbf{j}_0^+, \\mathbf{j}_0^-`` and
their derivatives are computed for the direct-solar contribution.

# Concepts page
See [Linearization — operator-level chain rule](../../docs/src/pages/concepts/06_linearization.md)
for the AD-boundary diagram, the parameter-strategy table, and the link
back to ParameterLayout for column ordering.

# Arguments
- `pol_type`: Polarization type (`Stokes_I`/`IQ`/`IQU`/`IQUV`).
- `SFI::Bool`: Whether to compute source-function integration terms.
- `τ_sum`: Cumulative optical depth above this layer `[nSpec]`.
- `τ̇_sum`: Derivative of cumulative τ w.r.t. parameters `[nSpec × Nparams]`.
- `dτ`: Elemental optical depth ``\\tau/2^{n_d}`` `[nSpec]`.
- `F₀`: Solar irradiance Stokes vector `[nStokes × nSpec]`.
- `computed_layer_properties`: Forward `CoreScatteringOpticalProperties`.
- `computed_layer_properties_lin`: `CoreScatteringOpticalPropertiesLin` =
  the AD-boundary handoff carrying ``(\\dot{\\tau}, \\dot{\\varpi}, \\dot{\\mathbf{Z}}^{++}, \\dot{\\mathbf{Z}}^{-+})``.
- `m::Int`: Fourier moment index.
- `ndoubl::Int`: Number of doublings.
- `scatter::Bool`: Whether the layer scatters.
- `quad_points`: Quadrature points and weights.
- `added_layer`, `added_layer_lin`: Forward + linearized RT matrices,
  written in place.
- `architecture`: `CPU`, `GPU`, or `MetalGPU`.
"""
function elemental!(pol_type, SFI::Bool,
                τ_sum::AbstractArray,#{FT2,1}, #Suniti
                τ̇_sum::AbstractArray,
                dτ::AbstractArray,
                F₀::AbstractArray,#{FT,2},    # Stokes vector of solar/stellar irradiance
                computed_layer_properties,
                computed_layer_properties_lin,
                m::Int,                     # m: fourier moment
                ndoubl::Int,                # ndoubl: number of doubling computations needed
                scatter::Bool,              # scatter: flag indicating scattering
                quad_points::QuadPoints{FT}, # struct with quadrature points, weights,
                added_layer::AddedLayer{FT},
                added_layer_lin::AddedLayerLin{FT},
                architecture;
                prepared_sources::AbstractSource =
                    SFI ? SourceSet((PreparedSolarBeam{FT,
                            typeof(array_type(architecture)(F₀))}(array_type(architecture)(F₀)),)) :
                          NoSource()
                ) where {FT<:AbstractFloat}

    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻) = added_layer
    (; ṙ⁺⁻, ṙ⁻⁺, ṫ⁻⁻, ṫ⁺⁺, J̇₀⁺, J̇₀⁻) = added_layer_lin
    (; qp_μ, iμ₀, μ₀, wt_μN, qp_μN) = quad_points
    (; τ, ϖ, Z⁺⁺, Z⁻⁺, Z₀⁺, Z₀⁻) = computed_layer_properties
    (; τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺, Ż₀⁺, Ż₀⁻) = computed_layer_properties_lin

    arr_type = array_type(architecture)
    qp_μN = arr_type(qp_μN)
    wt_μN = arr_type(wt_μN)
    τ_sum = arr_type(τ_sum)
    τ̇_sum = arr_type(τ̇_sum)
    I₀    = arr_type(pol_type.I₀)
    D = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))

    device = devi(architecture)

    # Chain-rule inputs (convert to device arrays)
    nparams = size(τ̇, 2)   # τ̇ is [nSpec, Nparams] — Nparams in last dim
    dτ̇_dev = arr_type(τ̇ ./ FT(2^ndoubl))   # elemental τ̇
    ϖ̇_dev  = arr_type(ϖ̇)
    Ż⁻⁺_dev = arr_type(Ż⁻⁺)
    Ż⁺⁺_dev = arr_type(Ż⁺⁺)

    # If in scattering mode:
    if scatter
   
        # for m==0, ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.5, while
        # for m>0,  ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.25  
        wct02 = fourier_weight(m, FT)
        wct2  = scaled_weights(m, wt_μN)
        # Zero forward, 3-core, AND ap_ arrays
        r⁻⁺ .= zero(FT)
        t⁺⁺ .= zero(FT)
        ṙ⁻⁺ .= zero(FT)
        ṫ⁺⁺ .= zero(FT)
        j₀⁺ .= zero(FT)
        j₀⁻ .= zero(FT)
        J̇₀⁺ .= zero(FT)
        J̇₀⁻ .= zero(FT)
        added_layer_lin.ap_ṫ⁺⁺ .= zero(FT)
        added_layer_lin.ap_ṫ⁻⁻ .= zero(FT)
        added_layer_lin.ap_ṙ⁻⁺ .= zero(FT)
        added_layer_lin.ap_ṙ⁺⁻ .= zero(FT)
        added_layer_lin.ap_J̇₀⁺ .= zero(FT)
        added_layer_lin.ap_J̇₀⁻ .= zero(FT)

        # Fused elemental + chain rule kernel
        kernel! = get_elem_rt_fused!(device)
        event = kernel!(r⁻⁺, t⁺⁺,
                    ṙ⁻⁺, ṫ⁺⁺,
                    added_layer_lin.ap_ṙ⁻⁺, added_layer_lin.ap_ṫ⁺⁺,
                    added_layer_lin.ap_ṙ⁺⁻, added_layer_lin.ap_ṫ⁻⁻,
                    ϖ, dτ, Z⁻⁺, Z⁺⁺,
                    dτ̇_dev, ϖ̇_dev, Ż⁻⁺_dev, Ż⁺⁺_dev,
                    qp_μN, wct2,
                    nparams, ndoubl, pol_type.n,
                    ndrange=size(r⁻⁺))
        synchronize_if_gpu()

        # Phase 3.5: source contributions dispatch via prepared_sources.
        # NoSource → no-op; PreparedSolarBeam → fused SFI+chain-rule+Bug-22
        # kernel; SourceSet → unrolled tuple iteration over each member.
        # The legacy `if SFI` block is gone; a single dispatch line handles
        # every scene. Bit-equal to the previous inline path when
        # prepared_sources is the default (SFI ? SolarBeam : NoSource).
        solar_external = Z₀⁺ !== nothing
        if solar_external && has_solar_beam(prepared_sources)
            columns = added_layer.solar_columns
            columns_lin = added_layer_lin.solar_columns
            (columns === nothing || columns_lin === nothing) && throw(ArgumentError(
                "external-solar linearization requires forward and tangent solar-column carriers"))
            Dpol = arr_type(pol_type.D)
            kernel! = get_elem_rt_solar_columns!(device)
            kernel!(columns.R₀⁻⁺, columns.R₀⁺⁻, columns.T₀⁺⁺, columns.T₀⁻⁻,
                    ϖ, dτ, Z₀⁻, Z₀⁺, qp_μN, μ₀, wct02, pol_type.n, Dpol,
                    ndrange=size(columns.R₀⁻⁺))
            synchronize_if_gpu()
            kernel! = get_elem_rt_solar_columns_lin!(device)
            kernel!(columns_lin.Ṙ₀⁻⁺, columns_lin.Ṙ₀⁺⁻,
                    columns_lin.Ṫ₀⁺⁺, columns_lin.Ṫ₀⁻⁻,
                    ϖ, dτ, Z₀⁻, Z₀⁺, dτ̇_dev, ϖ̇_dev,
                    arr_type(Ż₀⁻), arr_type(Ż₀⁺), qp_μN, μ₀,
                    wct02, pol_type.n, Dpol,
                    ndrange=(size(columns_lin.Ṙ₀⁻⁺,1), size(columns_lin.Ṙ₀⁻⁺,2),
                             size(columns_lin.Ṙ₀⁻⁺,3), nparams))
            synchronize_if_gpu()
            kernel! = apply_solar_columns_lin!(device)
            kernel!(j₀⁺, j₀⁻, added_layer_lin.ap_J̇₀⁺, added_layer_lin.ap_J̇₀⁻,
                    columns.R₀⁻⁺, columns.T₀⁺⁺,
                    columns_lin.Ṙ₀⁻⁺, columns_lin.Ṫ₀⁺⁺,
                    arr_type(F₀), τ_sum, τ̇_sum, μ₀, ndoubl, pol_type.n, Dpol,
                    ndrange=(size(added_layer_lin.ap_J̇₀⁺,1), 1,
                             size(added_layer_lin.ap_J̇₀⁺,3), nparams))
            synchronize_if_gpu()
        else
            source_tangent!(prepared_sources,
            j₀⁺, j₀⁻,
            J̇₀⁺, J̇₀⁻,
            added_layer_lin.ap_J̇₀⁺, added_layer_lin.ap_J̇₀⁻,
            ϖ, dτ,
            τ_sum, τ̇_sum,
            solar_external ? Z₀⁻ : Z⁻⁺, solar_external ? Z₀⁺ : Z⁺⁺,
            dτ̇_dev, ϖ̇_dev,
            solar_external ? arr_type(Ż₀⁻) : Ż⁻⁺_dev,
            solar_external ? arr_type(Ż₀⁺) : Ż⁺⁺_dev,
            qp_μN, ndoubl, wct02,
            pol_type.n, I₀, μ₀, solar_external ? 1 : iμ₀, D, nparams,
            architecture)
            synchronize_if_gpu()
        end

        # Apply D Matrix to forward quantities (fused kernel handles derivative D internally)
        apply_D_matrix_elemental!(ndoubl, pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

        # SFI D-matrix already applied inside fused kernel for sources;
        # this post-pass is a no-op for NoSource. Dispatch hides the gate.
        if has_solar_beam(prepared_sources)
            apply_D_matrix_elemental_SFI!(ndoubl, pol_type.n, j₀⁻)
        end
    else
        # No scattering: zero ap_ arrays, set transmission only
        added_layer_lin.ap_ṫ⁺⁺ .= zero(FT)
        added_layer_lin.ap_ṫ⁻⁻ .= zero(FT)
        added_layer_lin.ap_ṙ⁻⁺ .= zero(FT)
        added_layer_lin.ap_ṙ⁺⁻ .= zero(FT)
        added_layer_lin.ap_J̇₀⁺ .= zero(FT)
        added_layer_lin.ap_J̇₀⁻ .= zero(FT)

        t⁺⁺[:] = Diagonal{exp(-τ ./ qp_μN)}
        t⁻⁻[:] = Diagonal{exp(-τ ./ qp_μN)}
        ṫ⁺⁺[:, :, :, 1] = Diagonal{exp(-τ ./ qp_μN).*(-1 ./ qp_μN)}
        ṫ⁻⁻[:, :, :, 1] = Diagonal{exp(-τ ./ qp_μN).*(-1 ./ qp_μN)}

        # Chain rule for no-scatter: ap_ṫ = ṫ[1]*dτ̇ (only τ derivative matters)
        nspec_here = size(τ, 1)
        for iparam = 1:nparams
            for iλ = 1:nspec_here
                @views added_layer_lin.ap_ṫ⁺⁺[:,:,iλ,iparam] .= ṫ⁺⁺[:,:,iλ,1] .* dτ̇_dev[iλ,iparam]
                @views added_layer_lin.ap_ṫ⁻⁻[:,:,iλ,iparam] .= ṫ⁻⁻[:,:,iλ,1] .* dτ̇_dev[iλ,iparam]
            end
        end
    end
end

"""
    apply_D_elemental!(ndoubl, pol_n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻)

KernelAbstractions D-matrix symmetry kernel for linearized elemental R/T
operators. It applies the same Stokes parity signs to the forward matrices and
their three core derivative slots so reverse-direction operators remain
consistent with the elastic D-symmetry convention.
"""
@kernel function apply_D_elemental!(ndoubl, pol_n, 
                                r⁻⁺, @Const(t⁺⁺), r⁺⁻, t⁻⁻,
                                ṙ⁻⁺, @Const(ṫ⁺⁺), ṙ⁺⁻, ṫ⁻⁻)
    i, j, n = @index(Global, NTuple) #how best to do this for linearization? Is : okay, or should I use an iparam index?

    if ndoubl < 1
        ii = mod1(i, pol_n)
        jj = mod1(j, pol_n)
        if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2)) 
            r⁺⁻[i,j,n] = r⁻⁺[i,j,n]
            t⁻⁻[i,j,n] = t⁺⁺[i,j,n]
            ṙ⁺⁻[i,j,n,1] = ṙ⁻⁺[i,j,n,1]
            ṙ⁺⁻[i,j,n,2] = ṙ⁻⁺[i,j,n,2]
            ṙ⁺⁻[i,j,n,3] = ṙ⁻⁺[i,j,n,3]
            ṫ⁻⁻[i,j,n,1] = ṫ⁺⁺[i,j,n,1]
            ṫ⁻⁻[i,j,n,2] = ṫ⁺⁺[i,j,n,2]
            ṫ⁻⁻[i,j,n,3] = ṫ⁺⁺[i,j,n,3]
        else
            r⁺⁻[i,j,n] = -r⁻⁺[i,j,n] 
            t⁻⁻[i,j,n] = -t⁺⁺[i,j,n] 
            ṙ⁺⁻[i,j,n,1] = -ṙ⁻⁺[i,j,n,1] 
            ṙ⁺⁻[i,j,n,2] = -ṙ⁻⁺[i,j,n,2] 
            ṙ⁺⁻[i,j,n,3] = -ṙ⁻⁺[i,j,n,3] 
            ṫ⁻⁻[i,j,n,1] = -ṫ⁺⁺[i,j,n,1] 
            ṫ⁻⁻[i,j,n,2] = -ṫ⁺⁺[i,j,n,2] 
            ṫ⁻⁻[i,j,n,3] = -ṫ⁺⁺[i,j,n,3] 
        end
    else
        if mod1(i, pol_n) > 2
            r⁻⁺[i,j,n] = - r⁻⁺[i,j,n]
            ṙ⁻⁺[i,j,n,1] = - ṙ⁻⁺[i,j,n,1]
            ṙ⁻⁺[i,j,n,2] = - ṙ⁻⁺[i,j,n,2]
            ṙ⁻⁺[i,j,n,3] = - ṙ⁻⁺[i,j,n,3]
        end 
    end
    nothing
end

"""
    apply_D_elemental_SFI!(ndoubl, pol_n, J₀⁻, J̇₀⁻)

KernelAbstractions D-matrix symmetry kernel for linearized elemental source
vectors. It negates the upwelling `U/V` source components and their three core
derivative slots when source-vector D-symmetry must be applied outside the
doubling update.
"""
@kernel function apply_D_elemental_SFI!(ndoubl, pol_n, J₀⁻, J̇₀⁻)
    i, _, n = @index(Global, NTuple)
    
    if ndoubl>1
        if mod1(i, pol_n) > 2
            J₀⁻[i, 1, n] = - J₀⁻[i, 1, n]
            J̇₀⁻[i, 1, n, 1] = - J̇₀⁻[i, 1, n, 1]
            J̇₀⁻[i, 1, n, 2] = - J̇₀⁻[i, 1, n, 2]
            J̇₀⁻[i, 1, n, 3] = - J̇₀⁻[i, 1, n, 3]
        end 
    end
    nothing
end

function apply_D_matrix_elemental!(ndoubl::Int, n_stokes::Int, 
                    r⁻⁺::AbstractArray{FT,3}, 
                    t⁺⁺::AbstractArray{FT,3}, 
                    r⁺⁻::AbstractArray{FT,3}, 
                    t⁻⁻::AbstractArray{FT,3},
                    ṙ⁻⁺::AbstractArray{FT,4}, 
                    ṫ⁺⁺::AbstractArray{FT,4}, 
                    ṙ⁺⁻::AbstractArray{FT,4}, 
                    ṫ⁻⁻::AbstractArray{FT,4}) where {FT}
    device = devi(architecture(r⁻⁺))
    applyD_kernel! = apply_D_elemental!(device)
    event = applyD_kernel!(ndoubl,n_stokes, 
                        r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, 
                        ṙ⁻⁺, ṫ⁺⁺, ṙ⁺⁻, ṫ⁻⁻, 
                        ndrange=size(r⁻⁺));
    #wait(device, event);
    synchronize_if_gpu();
    return nothing
end

function apply_D_matrix_elemental_SFI!(ndoubl::Int, n_stokes::Int, 
                                J₀⁻::AbstractArray{FT,3},
                                J̇₀⁻::AbstractArray{FT,4}) where {FT}
    if ndoubl > 1
        return nothing
    else 
        device = devi(architecture(J₀⁻))
        applyD_kernel! = apply_D_elemental_SFI!(device)
        event = applyD_kernel!(ndoubl,n_stokes, J₀⁻, J̇₀⁻, ndrange=size(J₀⁻));
        #wait(device, event);
        synchronize_if_gpu();
        return nothing
    end
end
