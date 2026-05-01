#=
 
This file contains RT doubling-related functions
 
=#

"""
    $(FUNCTIONNAME)(pol_type, SFI, expk, ndoubl::Int, added_layer::AddedLayer, I_static::AbstractArray{FT}, 
                    architecture) where {FT}

Compute homogenous layer matrices from its elemental layer using Doubling 
"""
function doubling_helper!(pol_type, 
                          SFI, 
                          expk, 
                          ndoubl::Int, 
                          added_layer::M,
                          I_static::AbstractArray{FT}, 
                          architecture) where {FT,M}

    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻, temp1, temp2, temp1_ptr, temp2_ptr,
       dbl_gp_refl, dbl_j₁⁺, dbl_j₁⁻) = added_layer
    dev = devi(architecture)

    ndoubl == 0 && return nothing
    
    @timeit "doubling_allocs" begin
    tt⁺⁺_gp_refl = dbl_gp_refl === nothing ? similar(t⁺⁺) : dbl_gp_refl
    j₁⁺ = dbl_j₁⁺ === nothing ? similar(j₀⁺) : dbl_j₁⁺
    j₁⁻ = dbl_j₁⁻ === nothing ? similar(j₀⁻) : dbl_j₁⁻
    end
    #temp = similar(t⁺⁺)
    # Pointers to avoid memory allocation in CUBLAS routines
    #@timeit "Pointers" gp_ptrs   = CUBLAS.unsafe_strided_batch(gp_refl)
    #@timeit "Pointers" temp_ptrs = CUBLAS.unsafe_strided_batch(temp)
    # Loop over number of doublings
    for n = 1:ndoubl
        @timeit "Batch Inv Doubling" compute_geometric_progression!(temp1, tt⁺⁺_gp_refl, r⁻⁺, t⁺⁺, I_static, temp2, temp1_ptr, temp2_ptr)
        @timeit "source_update" doubling_source_update!(j₀⁺, j₀⁻, j₁⁺, j₁⁻, r⁻⁺, tt⁺⁺_gp_refl, expk)
        @timeit "rt_update" doubling_rt_update!(r⁻⁺, t⁺⁺, tt⁺⁺_gp_refl, expk)
    end
    @timeit "sync_doubling" synchronize_if_gpu()

    @timeit "apply_D_matrix" begin
    apply_D_matrix!(pol_type.n, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    apply_D_matrix_SFI!(pol_type.n, j₀⁻)
    end
#    CUBLAS.unsafe_free!(temp_ptrs);
#    CUBLAS.unsafe_free!(gp_ptrs);
    return nothing 
end

"""
    doubling!(pol_type, SFI, expk, ndoubl, added_layer, I_static, architecture)

Double the elemental layer `ndoubl` times to build the full homogeneous-layer
reflectance, transmission, and source matrices stored in `added_layer`.

Delegates to [`doubling_helper!`](@ref), which iteratively applies the
adding equations (see Eqs. 8 in the Raman paper draft) and then restores
the `D`-matrix symmetry.  After completion a GPU synchronisation barrier
is issued.

# Arguments
- `pol_type`: polarization type (determines `D`-matrix structure)
- `SFI`: whether Source Function Integration is active
- `expk`: `exp(-dτ/μ₀)` attenuation factor (doubled each iteration)
- `ndoubl::Int`: number of doubling steps
- `added_layer`: [`AddedLayer`](@ref) whose `r`, `t`, `j` fields are updated in-place
- `I_static`: pre-allocated batched identity matrix
- `architecture`: CPU or GPU selector
"""
function doubling!(pol_type, SFI, 
                    expk,
                    ndoubl::Int, 
                    added_layer::AddedLayer{M},#{FT},
                    I_static::AbstractArray{FT}, 
                    architecture) where {FT,M}

    doubling_helper!(pol_type, SFI, 
                    expk, ndoubl, added_layer, I_static, architecture)
    synchronize_if_gpu()
end

"""
    apply_D!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

KernelAbstractions kernel that recovers the four homogeneous-layer operators
from a single direction's matrices using the polarization D-matrix symmetry
of Sanghavi et al. (2014), JQSRT 133:412–433.

For a homogeneous layer with `D = diag(1, 1, -1, -1)` per stream:

    T_ab = D · T_ba · D       (Eq. 29)
    R_ab = D · R_ba · D       (Eq. 30)

vSmartMOM's [`doubling_helper!`](@ref) computes only one direction during the
inner loop using the *starred* quantity `R*_10 = D · R_10` (Eq. 31), which
halves the cost. After the loop, this kernel reconstructs the four operators
following Eq. (32):

    T_ba ← T_ba                 (no change)
    R_ba ← D · R*_ba
    T_ab ← D · T_ba · D
    R_ab ← R*_ba · D

The kernel does this in two in-place passes per `(iμ, jμ, n)` index:

1. Row-multiply `r⁻⁺` by `D` (negate rows i > 2). After this, `r⁻⁺`
   holds `R*_10 = D · R_10`.
2. Write the reverse-direction operators using the (i,j)-parity sign table
   for `D[i] · D[j] = ±1`. Same-parity (both ≤ 2 or both > 2) → +1;
   mixed parity → −1.

# Arguments
- `n_stokes::Int`: number of Stokes components carried (1, 3, or 4 for
  `Stokes_I`/`IQU`/`IQUV`); chosen from `pol_type.n`.
- `r⁻⁺::AbstractArray{FT,3}`: reflection (downward → upward). On entry, the
  value computed by the doubling loop. **Modified in place** — rows with
  Stokes index `i > 2` are negated. On exit, holds `D · R_10`.
- `t⁺⁺::AbstractArray{FT,3}`: transmission (downward). Read-only.
- `r⁺⁻::AbstractArray{FT,3}`: written from `r⁻⁺` with the parity sign rule.
- `t⁻⁻::AbstractArray{FT,3}`: written from `t⁺⁺` with the parity sign rule.

# Concepts page
See [The MOM Solver § Doubling](../../docs/src/pages/concepts/04_mom_solver.md)
for the equation derivation, a stream-by-stream worked example, and the
side-by-side mapping back to the doubling inner loop.
"""
@kernel function apply_D!(n_stokes::Int,  r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
    iμ, jμ, n = @index(Global, NTuple)
    i = mod1(iμ, n_stokes)
    j = mod1(jμ, n_stokes)

    # Pass 1: row-multiply r⁻⁺ by D = diag(1,1,-1,-1) — negate rows i > 2.
    # After this, r⁻⁺ holds R*_10 = D · R_10 (Sanghavi 2014, Eq. 31).
    if (i > 2)
        r⁻⁺[iμ,jμ,n] = - r⁻⁺[iμ, jμ,n]
    end

    # Pass 2: recover the four homogeneous-layer operators via Eq. (32),
    # using the (i,j)-parity table for D[i]·D[j]:
    #   same-parity (both ≤ 2 or both > 2) → D[i]·D[j] = +1
    #   mixed parity                       → D[i]·D[j] = -1
    if ((i <= 2) & (j <= 2)) | ((i > 2) & (j > 2))
        r⁺⁻[iμ,jμ,n] = r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = t⁺⁺[iμ,jμ,n]
    else
        r⁺⁻[iμ,jμ,n] = - r⁻⁺[iμ,jμ,n]
        t⁻⁻[iμ,jμ,n] = - t⁺⁺[iμ,jμ,n]
    end

end

"""
    apply_D_SFI!(n_stokes, J₀⁻)

Companion to [`apply_D!`](@ref) for the source-function-integration vector.
Negates the Stokes-`U`/`V` components (i > 2) of `J₀⁻` in place to apply
the D-matrix symmetry to the upwelling source vector.
"""
@kernel function apply_D_SFI!(n_stokes::Int, J₀⁻)
    iμ, _, n = @index(Global, NTuple)
    i = mod1(iμ, n_stokes)
    if (i > 2)
        J₀⁻[iμ, 1, n] = - J₀⁻[iμ, 1, n]
    end
end

"""
    apply_D_matrix!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)

Host-side launcher for [`apply_D!`](@ref). Selects the
KernelAbstractions backend from `architecture(r⁻⁺)` and invokes the kernel
over the full `(NquadN, NquadN, nSpec)` index space.

For scalar runs (`n_stokes == 1`) the polarization symmetry is trivial — the
two reverse-direction matrices are simple copies — so the kernel is bypassed.
"""
@inline function apply_D_matrix!(n_stokes::Int, r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3}) where {FT}
    if n_stokes == 1
        r⁺⁻ .= r⁻⁺
        t⁻⁻ .= t⁺⁺
        return nothing
    else
        device = devi(architecture(r⁻⁺))
        applyD_kernel! = apply_D!(device)
        event = applyD_kernel!(n_stokes, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, ndrange=size(r⁻⁺));
        #wait(device, event);
        synchronize_if_gpu();
        return nothing
    end
end


@inline function apply_D_matrix_SFI!(n_stokes::Int, J₀⁻::AbstractArray{FT,3}) where {FT}
    n_stokes == 1 && return nothing
    device = devi(architecture(J₀⁻))
    applyD_kernel! = apply_D_SFI!(device)
    event = applyD_kernel!(n_stokes, J₀⁻, ndrange=size(J₀⁻));
    #wait(device, event);
    synchronize_if_gpu();
    nothing
end
