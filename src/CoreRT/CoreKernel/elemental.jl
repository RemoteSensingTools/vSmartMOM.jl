#=
============================================================================
Elemental layer  (1/3 of the CoreKernel adding-doubling solver)
============================================================================

For a single thin homogeneous slab of optical thickness dτ this builds the
single-scattering reflection / transmission / source operators

    r⁻⁺ = (ϖ / 4μ) · Z⁻⁺ · dτ                  reflection (single-scatter)
    t⁺⁺ = I − [(I − ϖ·Z⁺⁺/4μ) · dτ]            transmission (1 − ext)
    j₀± = (ϖ · Z(μ, μ₀) / 4μ) · F₀ · e^(-τ̄/μ₀) solar source vectors (SFI)

where ϖ is the single-scattering albedo, Z± are the phase-matrix moments
assembled from the Fourier-decomposed Greek coefficients (see
`Scattering/compute_Z_matrices.jl`), μ₀ is the cosine of the solar zenith,
and F₀ is the incident solar flux.  The D-matrix symmetry relation is then
applied to recover `r⁺⁻` and `t⁻⁻` from `r⁻⁺` and `t⁺⁺`.

The slab is intentionally thin (≪ 1 mean free path); doubling repeatedly
combines two such slabs into one of double the thickness, and interaction
glues the resulting homogeneous layers together with the layer above /
below.  The trio is therefore:

    elemental!  →  doubling!  →  interaction!  (top → bottom of column)

Numerically stable for arbitrary absorption (Fell 1997, Eqs. 1.52–1.56).
Sanghavi et al. 2014, JQSRT 133:412–433, §3 gives the vector formulation
used here.  GPU kernels are dispatched through KernelAbstractions.
============================================================================
=#

"""
    elemental!(pol_type, SFI, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺, m, ndoubl, scatter, quad_points, added_layer, I_static, architecture)

Build the thin single-scattering elemental layer matrices (r, t, j) for a
homogeneous slab with optical depth `dτ`.

Uses exact single-scattering formulas (Fell 1997, Eqs. 1.52--1.56) that
remain numerically stable for arbitrary absorption.  GPU kernels are
dispatched via KernelAbstractions.  The `D`-matrix symmetry relation is
applied at the end to fill `r⁺⁻` and `t⁻⁻` from `r⁻⁺` and `t⁺⁺`.

When `SFI == true`, the solar source terms `J₀⁺` and `J₀⁻` are computed;
otherwise they are left untouched.

# Arguments
- `pol_type`: polarization type (Stokes_I, Stokes_IQ, Stokes_IQU, or Stokes_IQUV)
- `SFI::Bool`: use Source Function Integration for the solar beam
- `τ_sum`: cumulative optical depth above the current layer (per λ)
- `dτ_λ`: elemental-layer total optical depth (per λ, including gas absorption)
- `dτ`: elemental-layer scattering optical depth (scalar)
- `ϖ_λ`: single-scattering albedo per λ (with gas absorption)
- `ϖ`: single-scattering albedo (scattering only, no gas absorption)
- `Z⁺⁺`, `Z⁻⁺`: forward/backward scattering phase-matrix expansions
- `m::Int`: Fourier azimuth moment index
- `ndoubl::Int`: number of subsequent doublings (affects `D`-matrix application)
- `scatter::Bool`: whether the layer scatters
- `quad_points`: [`QuadPoints`](@ref) with quadrature angles and weights
- `added_layer`: [`AddedLayer`](@ref) or `AddedLayerRS` to fill in-place
- `I_static`: pre-allocated identity matrix (batched)
- `architecture`: CPU or GPU architecture selector
"""
function elemental!(pol_type, SFI::Bool, 
                            τ_sum::AbstractArray,       #{FT2,1}, #Suniti
                            dτ_λ::AbstractArray{FT,1},  # dτ_λ: total optical depth of elemental layer (per λ)
                            dτ::FT,                     # dτ:   scattering optical depth of elemental layer (scalar)
                            ϖ_λ::AbstractArray{FT,1},   # ϖ_λ: single scattering albedo of elemental layer (per λ, absorptions by gases included)
                            ϖ::FT,                      # ϖ: single scattering albedo of elemental layer (no trace gas absorption included)
                            Z⁺⁺::AbstractArray{FT,2},   # Z matrix
                            Z⁻⁺::AbstractArray{FT,2},   # Z matrix
                            F₀::AbstractArray,          # F₀: solar irradiance Stokes vector (nStokes, nSpec)
                            m::Int,                     # m: fourier moment
                            ndoubl::Int,                # ndoubl: number of doubling computations needed
                            scatter::Bool,              # scatter: flag indicating scattering
                            quad_points::QuadPoints{FT2}, # struct with quadrature points, weights,
                            added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}},
                            I_static,
                            architecture) where {FT<:Real,FT2}

    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻) = added_layer
    (; qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀) = quad_points
    #@unpack ϖ_Cabannes = RS_type
    arr_type = array_type(architecture)
    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM
    # @show Array(τ_sum)[1], Array(dτ_λ)[1], Array(ϖ_λ)[1], Array(Z⁺⁺)[1,1]
    # Later on, we can have Zs also vary with index, pretty easy here:
    # Z⁺⁺_ = repeat(Z⁺⁺, 1, 1, 1)
    Z⁺⁺_ = reshape(Z⁺⁺, (size(Z⁺⁺,1), size(Z⁺⁺,2),1))
    # Z⁻⁺_ = repeat(Z⁻⁺, 1, 1, 1)
    Z⁻⁺_ = reshape(Z⁻⁺, (size(Z⁺⁺,1), size(Z⁺⁺,2),1))

    D = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))
    I₀_NquadN = arr_type(zeros(FT,size(qp_μN,1))); #incident irradiation
    i_end     = pol_type.n*iμ₀
    I₀_NquadN[iμ₀Nstart:i_end] = pol_type.I₀

    device = devi(architecture)

    # If in scattering mode:
    if scatter
   
        NquadN = length(qp_μN)

        # Needs explanation still, different weights: 
        # for m==0, ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.5, while
        # for m>0,  ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.25  
        wct02 = fourier_weight(m, FT)
        wct2  = scaled_weights(m, wt_μN)
        wct0  = wct02 * ϖ * dτ
        wct   = wct02 * ϖ * wt_μN

        # Get the diagonal matrices first
        d_qp  = Diagonal(1 ./ qp_μN)
        d_wct = Diagonal(wct)

        # Calculate r⁻⁺ and t⁺⁺
        kernel! = get_elem_rt!(device)
        event = kernel!(r⁻⁺, t⁺⁺, ϖ_λ, dτ_λ, Z⁻⁺, Z⁺⁺,
            qp_μN, wct2, ndrange=size(r⁻⁺));
        synchronize_if_gpu()

        if SFI
            kernel! = get_elem_rt_SFI!(device)
            event = kernel!(J₀⁺, J₀⁻, ϖ_λ, dτ_λ, τ_sum, Z⁻⁺, Z⁺⁺, arr_type(F₀), qp_μN, ndoubl, wct02, pol_type.n, arr_type(pol_type.I₀), iμ₀, D, ndrange=size(J₀⁺))
            synchronize_if_gpu()
        end

        # Apply D Matrix
        apply_D_matrix_elemental!(ndoubl, pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

        if SFI
            apply_D_matrix_elemental_SFI!(ndoubl, pol_type.n, J₀⁻)
        end      
    else 
        # Note: τ is not defined here
        t⁺⁺[:] = Diagonal{exp(-τ ./ qp_μN)}
        t⁻⁻[:] = Diagonal{exp(-τ ./ qp_μN)}
    end    
    #@pack! added_layer = r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, J₀⁺, J₀⁻   
end

"""
    elemental!(pol_type, SFI, τ_sum, dτ, computed_layer_properties, m, ndoubl, scatter, quad_points, added_layer, architecture)

Variant of [`elemental!`](@ref) that accepts a
[`CoreScatteringOpticalProperties`](@ref) bundle instead of individual
`ϖ_λ`, `Z⁺⁺`, `Z⁻⁺` arrays.  Optical properties are unpacked internally
and the same GPU-accelerated single-scattering formulas are applied.
"""
function elemental!(pol_type, SFI::Bool,
                            τ_sum::AbstractArray,#{FT2,1}, #Suniti
                            dτ::AbstractArray,
                            F₀::AbstractArray,          # F₀: solar irradiance Stokes vector (nStokes, nSpec)
                            computed_layer_properties::CoreScatteringOpticalProperties,
                            m::Int,                     # m: fourier moment
                            ndoubl::Int,                # ndoubl: number of doubling computations needed 
                            scatter::Bool,              # scatter: flag indicating scattering
                            quad_points::QuadPoints{FT2}, # struct with quadrature points, weights, 
                            added_layer::Union{AddedLayer{FT},AddedLayerRS{FT}}, 
                            architecture) where {FT<:Real,FT2}

    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺) = added_layer
    j₀⁺ = added_layer.j₀⁺
    j₀⁻ = added_layer.j₀⁻
    (; qp_μ, iμ₀, wt_μN, qp_μN) = quad_points
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = computed_layer_properties
    #@show M
    arr_type = array_type(architecture)

    # Need to check with paper nomenclature. This is basically eqs. 19-20 in vSmartMOM
    I₀    = arr_type(pol_type.I₀)
    D     = Diagonal(arr_type(repeat(pol_type.D, size(qp_μ,1))))

    device = devi(architecture)

    # If in scattering mode:
    if scatter
        # for m==0, ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.5, while
        # for m>0,  ₀∫²ᵖⁱ cos²(mϕ)dϕ/4π = 0.25  
        wct02 = fourier_weight(m, FT)
        wct2  = scaled_weights(m, wt_μN)
 
        # More computationally intensive definition of a single scattering layer with variable (0-∞) absorption
        # with absorption in batch mode, low tau_scatt but higher tau_total, needs exact equations
        kernel! = get_elem_rt!(device)
        event = kernel!(r⁻⁺, t⁺⁺, ϖ, dτ, Z⁻⁺, Z⁺⁺, qp_μN, wct2, ndrange=size(r⁻⁺)); 
        #wait(device, event)
        synchronize_if_gpu()

        # SFI part
        kernel! = get_elem_rt_SFI!(device)
        event = kernel!(j₀⁺, j₀⁻, ϖ, dτ, arr_type(τ_sum), Z⁻⁺, Z⁺⁺, arr_type(F₀), qp_μN, ndoubl, wct02, pol_type.n, I₀, iμ₀, D, ndrange=size(j₀⁺))
        #wait(device, event)
        synchronize_if_gpu()
        
        # Apply D Matrix
        apply_D_matrix_elemental!(ndoubl, pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

        # apply D matrix for SFI
        apply_D_matrix_elemental_SFI!(ndoubl, pol_type.n, j₀⁻)   
    else
        # Note: τ is not defined here
        t⁺⁺ .= Diagonal{exp(-τ ./ qp_μN)}
        t⁻⁻ .= Diagonal{exp(-τ ./ qp_μN)}
    end    
end

"""
    get_elem_rt!(r⁻⁺, t⁺⁺, ϖ_λ, dτ_λ, Z⁻⁺, Z⁺⁺, μ, wct)

KernelAbstractions kernel that fills the elemental layer's reflection
`r⁻⁺` and transmission `t⁺⁺` matrices using the **exact finite-δ
single-scatter formulas** of Fell (1997, FU Berlin PhD thesis), Eqs.
1.52–1.56, restated in Sanghavi & Frankenberg (2023, JQSRT 311:108791),
Eqs. (10)–(11):

    r⁻⁺[i,j] = ϖ_λ · Z⁻⁺[i,j] · (μⱼ/(μᵢ+μⱼ)) · wⱼ · (1 - exp(-δτ(1/μᵢ + 1/μⱼ)))
    t⁺⁺[i,j] = ϖ_λ · Z⁺⁺[i,j] · (μⱼ/(μᵢ-μⱼ)) · wⱼ · (exp(-δτ/μᵢ) - exp(-δτ/μⱼ))     (i ≠ j)
    t⁺⁺[i,i] = exp(-δτ/μᵢ) · (1 + ϖ_λ · Z⁺⁺[i,i] · (δτ/μᵢ) · wᵢ)                    (i = j, L'Hôpital limit)

These are NOT the infinitesimal-δ linearization of Sanghavi et al. (2014,
JQSRT 133:412–433) Eqs. (19)–(20) used by many other MOM codes. Those are
linear in δ and require very thin elemental layers (large N_doubl) for
accuracy. The exact finite-δ form here lets the elemental layer be much
thicker at the same single-scatter accuracy → smaller N_doubl → less
round-off through doubling, stable in Float32 on GPU.

The numerical-stability discipline:
- `1 - exp(-x)` is computed as `-expm1(-x)`, accurate when x is small.
- `exp(-a) - exp(-b)` is computed as `expdiff_neg(a, b)`, accurate when
  a ≈ b. This case occurs when nearby μ_i, μ_j make the difference of
  exponentials catastrophically cancel in the naive evaluation.
- The `μ[i] == μ[j]` branch handles the diagonal limit where the i ≠ j
  formula has an `(μⱼ/(μᵢ-μⱼ))` apparent singularity. The L'Hôpital limit
  used here is the closed-form replacement.

Streams with quadrature weight below `eps(FT)` (e.g. zero-weight viewing
nodes appended for direct ray-tracing in the Block-Radau scheme) get
trivial Beer-law transmission and zero reflection.

# Arguments
- `r⁻⁺::AbstractArray{FT,3}`: reflection matrix `(NquadN, NquadN, nSpec)`,
  written in place.
- `t⁺⁺::AbstractArray{FT,3}`: transmission matrix `(NquadN, NquadN, nSpec)`,
  written in place.
- `ϖ_λ::AbstractArray{FT,1}`: per-wavelength SSA (= τ_scat·ϖ / τ_λ from
  the layer-optics build); shape `(nSpec,)`.
- `dτ_λ::AbstractArray{FT,1}`: per-wavelength **total** elemental optical
  depth (absorption included); shape `(nSpec,)`. **Note** that this is
  the total τ_λ at the elemental thickness — the doubling count was sized
  by the scattering-only τ_scat in `get_dtau_ndoubl` (rt_kernel.jl:245).
  See `docs/src/pages/concepts/03c_mixing.md`.
- `Z⁻⁺`, `Z⁺⁺`: backscatter and forward Fourier-moment phase matrices,
  shape `(NquadN, NquadN, 1)` or `(NquadN, NquadN, nSpec)`. The kernel
  detects which.
- `μ`: quadrature angle cosines.
- `wct`: quadrature-weight factors carrying the polarization-stream count
  (Stokes_n) for index conversion.

# Concepts page
See [The MOM Solver § Elemental layer](../../docs/src/pages/concepts/04_mom_solver.md#elemental-layer)
for the side-by-side comparison with the linear S2014 limit, the
worked-example O₂ A-band layer, and the `expm1`/`expdiff_neg` rationale.
"""
@kernel function get_elem_rt!(r⁻⁺, t⁺⁺, @Const(ϖ_λ), @Const(dτ_λ), @Const(Z⁻⁺), @Const(Z⁺⁺), @Const(μ), @Const(wct))
    FT = eltype(r⁻⁺)
    n2 = 1
    i, j, n = @index(Global, NTuple)
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    if (wct[j] > eps(FT)) 
        # 𝐑⁻⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁻⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ+μⱼ)) ̇(1 - exp{-τ ̇(1/μᵢ + 1/μⱼ)}) ̇𝑤ⱼ
        r⁻⁺[i,j,n] = 
            ϖ_λ[n] * Z⁻⁺[i,j,n2] * 
            #Z⁻⁺[i,j] * 
            (μ[j] / (μ[i] + μ[j])) * wct[j] * 
            -expm1(-dτ_λ[n] * ((1 / μ[i]) + (1 / μ[j])))
                    
        if (μ[i] == μ[j])
            # 𝐓⁺⁺(μᵢ, μᵢ) = (exp{-τ/μᵢ} + ϖ ̇𝐙⁺⁺(μᵢ, μᵢ) ̇(τ/μᵢ) ̇exp{-τ/μᵢ}) ̇𝑤ᵢ
            if i == j
                t⁺⁺[i,j,n] = 
                    exp(-dτ_λ[n] / μ[i]) *
                    (1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * (dτ_λ[n] / μ[i]) * wct[i])
                    #(1 + ϖ_λ[n] * Z⁺⁺[i,i] * (dτ_λ[n] / μ[i]) * wct[i])
            else
                t⁺⁺[i,j,n] = exp(-dτ_λ[n] / μ[j]) *
                    (ϖ_λ[n] * Z⁺⁺[i,j,n2] * (dτ_λ[n] / μ[i]) * wct[j])
            end
        else
    
            # 𝐓⁺⁺(μᵢ, μⱼ) = ϖ ̇𝐙⁺⁺(μᵢ, μⱼ) ̇(μⱼ/(μᵢ-μⱼ)) ̇(exp{-τ/μᵢ} - exp{-τ/μⱼ}) ̇𝑤ⱼ
            # (𝑖 ≠ 𝑗)
            t⁺⁺[i,j,n] = 
                ϖ_λ[n] * Z⁺⁺[i,j,n2] * 
                #Z⁺⁺[i,j] * 
                (μ[j] / (μ[i] - μ[j])) * wct[j] * 
                expdiff_neg(dτ_λ[n] / μ[i], dτ_λ[n] / μ[j])
        end
    else
        r⁻⁺[i,j,n] = zero(FT)
        if i==j
            t⁺⁺[i,j,n] = exp(-dτ_λ[n] / μ[i]) #Suniti
        else
            t⁺⁺[i,j,n] = zero(FT)
        end
    end
    nothing
end

"""
    get_elem_rt_SFI!(J₀⁺, J₀⁻, ϖ_λ, dτ_λ, τ_sum, Z⁻⁺, Z⁺⁺, F₀, μ,
                     ndoubl, wct02, nStokes, I₀, iμ0, D)

KernelAbstractions elemental source-function kernel for the direct solar beam.
Each workitem owns one stream/spectral element `(i, n)`, forms the local
`Z⁺⁺ * F₀` and `Z⁻⁺ * F₀` Stokes contractions for the incident solar stream,
then writes the exact finite-δ downwelling and upwelling source vectors. The
kernel applies cumulative beam attenuation `exp(-τ_sum / μ₀)` and the
D-matrix sign for upwelling source terms when the elemental layer will be
doubled.
"""
@kernel function get_elem_rt_SFI!(J₀⁺, J₀⁻, @Const(ϖ_λ), @Const(dτ_λ), @Const(τ_sum), @Const(Z⁻⁺), @Const(Z⁺⁺), @Const(F₀), @Const(μ), ndoubl, wct02, nStokes, @Const(I₀), iμ0, @Const(D))
    i_start  = nStokes*(iμ0-1) + 1 
    i_end    = nStokes*iμ0
    
    i, _, n = @index(Global, NTuple) ##Suniti: What are Global and Ntuple?
    FT = eltype(I₀)
    J₀⁺[i, 1, n] = zero(FT)
    J₀⁻[i, 1, n] = zero(FT)
    n2=1
    if size(Z⁻⁺,3)>1
        n2 = n
    end
    
    Z⁺⁺_I₀ = zero(FT);
    Z⁻⁺_I₀ = zero(FT);
    
    for ii = i_start:i_end
        # `n2` tracks the Z matrix spectral axis, which can be length 1 for
        # m-invariant Rayleigh/Cabannes fast paths.  The solar source is always
        # defined on the radiance grid, so keep its spectral index at `n`.
        Z⁺⁺_I₀ += Z⁺⁺[i,ii,n2] * F₀[ii-i_start+1,n]
        Z⁻⁺_I₀ += Z⁻⁺[i,ii,n2] * F₀[ii-i_start+1,n]
    end

    if (i >= i_start) & (i <= i_end)
        ctr = i-i_start+1
        # See Eq. 1.54 in Fell
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * (dτ(λ)/μ₀) * exp(-dτ(λ)/μ₀)
        J₀⁺[i, 1, n] = wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * (dτ_λ[n] / μ[i]) * exp(-dτ_λ[n] / μ[i])
    else
        # J₀⁺ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁺⁺ * I₀ * [μ₀ / (μᵢ - μ₀)] * [exp(-dτ(λ)/μᵢ) - exp(-dτ(λ)/μ₀)]
        # See Eq. 1.53 in Fell
        J₀⁺[i, 1, n] = 
        wct02 * ϖ_λ[n] * Z⁺⁺_I₀ * (μ[i_start] / (μ[i] - μ[i_start])) * 
        expdiff_neg(dτ_λ[n] / μ[i], dτ_λ[n] / μ[i_start])
    end
    #J₀⁻ = 0.25*(1+δ(m,0)) * ϖ(λ) * Z⁻⁺ * I₀ * [μ₀ / (μᵢ + μ₀)] * [1 - exp{-dτ(λ)(1/μᵢ + 1/μ₀)}]
    # See Eq. 1.52 in Fell
    J₀⁻[i, 1, n] = wct02 * ϖ_λ[n] * Z⁻⁺_I₀ * (μ[i_start] / (μ[i] + μ[i_start])) * (-expm1(-dτ_λ[n] * ((1 / μ[i]) + (1 / μ[i_start]))))

    J₀⁺[i, 1, n] *= exp(-τ_sum[n]/μ[i_start])
    J₀⁻[i, 1, n] *= exp(-τ_sum[n]/μ[i_start])

    if ndoubl >= 1
        J₀⁻[i, 1, n] = D[i,i]*J₀⁻[i, 1, n] #D = Diagonal{1,1,-1,-1,...Nquad times}
    end  
    nothing
end

"""
    apply_D_elemental!(ndoubl, pol_n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

KernelAbstractions D-matrix symmetry kernel for the elemental elastic layer.
For undoubled layers it fills the reverse-direction operators `r⁺⁻` and `t⁻⁻`
from `r⁻⁺` and `t⁺⁺` using the Stokes parity signs. For layers that will be
doubled it applies the row sign to `r⁻⁺` in place so the doubling step can use
the cheaper symmetric operator form.
"""
@kernel function apply_D_elemental!(ndoubl, pol_n, r⁻⁺, @Const(t⁺⁺), r⁺⁻, t⁻⁻)
    i, j, n = @index(Global, NTuple)

    if ndoubl < 1
        ii = mod1(i, pol_n)
        jj = mod1(j, pol_n)
        if ((ii <= 2) & (jj <= 2)) | ((ii > 2) & (jj > 2))
            r⁺⁻[i,j,n] = r⁻⁺[i,j,n]
            t⁻⁻[i,j,n] = t⁺⁺[i,j,n]
        else
            r⁺⁻[i,j,n] = -r⁻⁺[i,j,n]
            t⁻⁻[i,j,n] = -t⁺⁺[i,j,n]
        end
    else
        if mod1(i, pol_n) > 2
            r⁻⁺[i,j,n] = - r⁻⁺[i,j,n]
        end
    end
    nothing
end

"""
    apply_D_elemental_SFI!(ndoubl, pol_n, J₀⁻)

KernelAbstractions D-matrix symmetry kernel for the elemental solar source
vector. It negates the upwelling `U/V` Stokes rows of `J₀⁻` when the
source-vector symmetry has not already been handled by the doubling path.
"""
@kernel function apply_D_elemental_SFI!(ndoubl, pol_n, J₀⁻)
    i, _, n = @index(Global, NTuple)

    if ndoubl>1
        if mod1(i, pol_n) > 2
            J₀⁻[i, 1, n] = - J₀⁻[i, 1, n]
        end 
    end
    nothing
end

function apply_D_matrix_elemental!(ndoubl::Int, n_stokes::Int, r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3}) where {FT}
    device = devi(architecture(r⁻⁺))
    applyD_kernel! = apply_D_elemental!(device)
    event = applyD_kernel!(ndoubl,n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
    #wait(device, event);
    synchronize_if_gpu();
    return nothing
end

function apply_D_matrix_elemental_SFI!(ndoubl::Int, n_stokes::Int, J₀⁻::AbstractArray{FT,3}) where {FT}
    if ndoubl > 1
        return nothing
    else 
        device = devi(architecture(J₀⁻))
        applyD_kernel! = apply_D_elemental_SFI!(device)
        event = applyD_kernel!(ndoubl,n_stokes, J₀⁻, ndrange=size(J₀⁻));
        #wait(device, event);
        synchronize_if_gpu();
        return nothing
    end
end
