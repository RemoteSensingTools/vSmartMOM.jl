#=
============================================================================
Interaction  (3/3 of the CoreKernel adding-doubling solver)
============================================================================

Combines the *composite* layer accumulated above (R, T, J — uppercase) with
the freshly *added* homogeneous layer below (r, t, j — lowercase) using the
matrix-operator-method adding equations (Sanghavi et al. 2014, JQSRT
133:412–433, Eqs. 23–28):

    G    = (E − r⁻⁺ · R⁺⁻)⁻¹                  geometric-series factor
    R'⁻⁺ = R⁻⁺ + T⁻⁻ · G · r⁻⁺ · T⁺⁺
    T'⁺⁺ = t⁺⁺ · G · T⁺⁺
    J'₀± = J₀± + T · G · (r · J₀∓ + j₀±)      source cascade

Called once per atmospheric layer as the column is assembled top-to-bottom.
The `R⁺⁻ / T⁻⁻` half is recovered by D-matrix symmetry (vector case) or
direct copy (scalar case).

Four `ScatteringInterface_{ab}` traits dispatch to specialized kernels:
`_00` (neither layer scatters; Beer-law product, no inverse), `_01` /
`_10` (only one side scatters; one of r and R⁺⁻ is zero, no inverse), and
`_11` (both scatter; full adding equations).  This trait-dispatch pattern
keeps every branch of the algorithm in its own type-stable method instead
of an `if`-tree inside one large function.

See `elemental.jl` and `doubling.jl` for the upstream layer construction,
and `docs/src/pages/concepts/04_mom_solver.md` for the prose walkthrough.
============================================================================
=#

"""
    interaction_helper!(::ScatteringInterface_00, SFI, composite_layer, added_layer, I_static)

Specialization of the matrix-operator-method **adding** step (Sanghavi
et al. 2014, JQSRT 133:412–433, Eqs. 23–28) for the case where neither the
composite layer above nor the newly added layer below scatters.

Both R-fields are zero by construction; the geometric-series factor
`(E − r⁻⁺ R⁺⁻)⁻¹` collapses to the identity. The full adding equations
reduce to:

    J₀⁺ ← j₀⁺ + t⁺⁺ · J₀⁺                  (downward source cascade)
    J₀⁻ ← J₀⁻ + T⁻⁻ · j₀⁻                  (upward source cascade)
    T⁻⁻ ← t⁻⁻ · T⁻⁻
    T⁺⁺ ← t⁺⁺ · T⁺⁺                        (Beer-law transmission only)

No matrix inversion is required. See [`interaction_helper!(::ScatteringInterface_11, ...)`](@ref)
for the general case and `docs/src/pages/concepts/04_mom_solver.md` for the
prose walkthrough.
"""
function interaction_helper!(::ScatteringInterface_00, SFI,
                                composite_layer::CompositeLayer{FT},
                                added_layer::AddedLayer{FT},
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}
    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻, j₀_by_src) = added_layer
    (; R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻, J₀_by_src) = composite_layer
    mA, v1, v2 = _interaction_scratch(added_layer)

    # Source Function — legacy solar slot (in-place Phase 0; same GEMMs and
    # broadcasts in the same order as the allocating form — bit-identical)
    _bmm!(v1, t⁺⁺, J₀⁺); J₀⁺ .= j₀⁺ .+ v1
    _bmm!(v1, T⁻⁻, j₀⁻); J₀⁻ .= J₀⁻ .+ v1
    # Per-source slots (same formula; uses pre-mutation T⁻⁻)
    for (key, slot) in pairs(j₀_by_src)
        cslot = J₀_by_src[key]
        _bmm!(v1, t⁺⁺, cslot.J₀⁺); cslot.J₀⁺ .= slot.j₀⁺ .+ v1
        _bmm!(v1, T⁻⁻, slot.j₀⁻);  cslot.J₀⁻ .= cslot.J₀⁻ .+ v1
    end

    # Batched multiplication between added and composite
    _bmm!(mA, t⁻⁻, T⁻⁻); T⁻⁻ .= mA
    _bmm!(mA, t⁺⁺, T⁺⁺); T⁺⁺ .= mA
end

"""
    interaction_helper!(::ScatteringInterface_01, SFI, composite_layer, added_layer, I_static)

Specialization of the matrix-operator-method **adding** step (Sanghavi
et al. 2014, JQSRT 133:412–433, Eqs. 23–28) for the case where only the
newly added (lower) layer scatters; the composite layer above is
non-scattering so `R⁺⁻ = 0` and the geometric-series factor
`(E − r⁻⁺ R⁺⁻)⁻¹` collapses to the identity.

The composite reflectance now picks up the added layer's `r⁻⁺` carried
through the composite transmission both ways:

    J₀⁻ ← J₀⁻ + T⁻⁻ · (r⁻⁺ J₀⁺ + j₀⁻)
    J₀⁺ ← j₀⁺ + t⁺⁺ · J₀⁺
    R⁻⁺ ← T⁻⁻ · r⁻⁺ · T⁺⁺                   (Eq. 24 with R₀₁ = 0)
    R⁺⁻ ← r⁺⁻
    T⁺⁺ ← t⁺⁺ · T⁺⁺
    T⁻⁻ ← T⁻⁻ · t⁻⁻

No matrix inversion is required. See [`interaction_helper!(::ScatteringInterface_11, ...)`](@ref)
for the general case.
"""
function interaction_helper!(::ScatteringInterface_01, SFI,
                                composite_layer::CompositeLayer{FT},
                                added_layer::AddedLayer{FT},
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}
    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻, j₀_by_src) = added_layer
    (; R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻, J₀_by_src) = composite_layer
    mA, v1, v2 = _interaction_scratch(added_layer)

    # Source Function — legacy solar slot (in-place Phase 0)
    _bmm!(v1, r⁻⁺, J₀⁺); v1 .= v1 .+ j₀⁻
    _bmm!(v2, T⁻⁻, v1);  J₀⁻ .= J₀⁻ .+ v2
    _bmm!(v1, t⁺⁺, J₀⁺); J₀⁺ .= j₀⁺ .+ v1
    # Per-source slots (uses pre-mutation T⁻⁻ and r⁻⁺)
    for (key, slot) in pairs(j₀_by_src)
        cslot = J₀_by_src[key]
        _bmm!(v1, r⁻⁺, cslot.J₀⁺); v1 .= v1 .+ slot.j₀⁻
        _bmm!(v2, T⁻⁻, v1);        cslot.J₀⁻ .= cslot.J₀⁻ .+ v2
        _bmm!(v1, t⁺⁺, cslot.J₀⁺); cslot.J₀⁺ .= slot.j₀⁺ .+ v1
    end

    # Batched multiplication between added and composite
    _bmm!(mA, T⁻⁻, r⁻⁺); _bmm!(R⁻⁺, mA, T⁺⁺)   # R⁻⁺ is a pure output here
    R⁺⁻ .= r⁺⁻
    _bmm!(mA, t⁺⁺, T⁺⁺); T⁺⁺ .= mA
    _bmm!(mA, T⁻⁻, t⁻⁻); T⁻⁻ .= mA
end

"""
    interaction_helper!(::ScatteringInterface_10, SFI, composite_layer, added_layer, I_static)

Specialization of the matrix-operator-method **adding** step (Sanghavi
et al. 2014, JQSRT 133:412–433, Eqs. 23–28) for the case where only the
composite (upper) layer scatters; the newly added layer below is
non-scattering so `r⁻⁺ = 0` and the geometric-series factor
`(E − R⁺⁻ r⁻⁺)⁻¹` collapses to the identity.

`R⁻⁺` is unchanged because nothing below it reflects; `R⁺⁻` is the existing
composite reflectance carried through the new layer's Beer-law transmission:

    J₀⁺ ← j₀⁺ + t⁺⁺ · (J₀⁺ + R⁺⁻ j₀⁻)
    J₀⁻ ← J₀⁻ + T⁻⁻ · j₀⁻
    T⁺⁺ ← t⁺⁺ · T⁺⁺
    T⁻⁻ ← T⁻⁻ · t⁻⁻
    R⁺⁻ ← t⁺⁺ · R⁺⁻ · t⁻⁻                  (Eq. 26 with R₂₁ = 0)

No matrix inversion is required. See [`interaction_helper!(::ScatteringInterface_11, ...)`](@ref)
for the general case.
"""
function interaction_helper!(::ScatteringInterface_10, SFI,
                                composite_layer::CompositeLayer{FT},
                                added_layer::AddedLayer{FT},
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}
    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻, j₀_by_src) = added_layer
    (; R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻, J₀_by_src) = composite_layer
    mA, v1, v2 = _interaction_scratch(added_layer)

    # Source Function — legacy solar slot (in-place Phase 0)
    _bmm!(v1, R⁺⁻, j₀⁻); v1 .= J₀⁺ .+ v1
    _bmm!(v2, t⁺⁺, v1);  J₀⁺ .= j₀⁺ .+ v2
    _bmm!(v1, T⁻⁻, j₀⁻); J₀⁻ .= J₀⁻ .+ v1
    # Per-source slots (uses pre-mutation R⁺⁻ and T⁻⁻)
    for (key, slot) in pairs(j₀_by_src)
        cslot = J₀_by_src[key]
        _bmm!(v1, R⁺⁻, slot.j₀⁻); v1 .= cslot.J₀⁺ .+ v1
        _bmm!(v2, t⁺⁺, v1);       cslot.J₀⁺ .= slot.j₀⁺ .+ v2
        _bmm!(v1, T⁻⁻, slot.j₀⁻); cslot.J₀⁻ .= cslot.J₀⁻ .+ v1
    end

    # Batched multiplication between added and composite
    _bmm!(mA, t⁺⁺, T⁺⁺); T⁺⁺ .= mA
    _bmm!(mA, T⁻⁻, t⁻⁻); T⁻⁻ .= mA
    _bmm!(mA, t⁺⁺, R⁺⁻); _bmm!(R⁺⁻, mA, t⁻⁻)   # R⁺⁻ read before overwrite
end

"""
    interaction_helper!(::ScatteringInterface_11, SFI, composite_layer, added_layer, I_static)

Full matrix-operator-method **adding** step: both the composite layer above
and the newly added layer below scatter. This is the algebraic heart of MOM
and the reason the method handles arbitrarily inhomogeneous atmospheres.

Implements Sanghavi et al. (2014, JQSRT 133:412–433), Eqs. (23)–(28),
restated for the layered atmosphere of Sanghavi & Frankenberg (2023, JQSRT
311:108791) Eq. (12). With the convention that the added (lower) layer is
labelled `21` and the composite (upper) layer is labelled `01`, two batched
matrix inversions form the geometric-series factors that capture infinitely
many reflections between the two layers:

    T01_inv = T⁻⁻ · (E − r⁻⁺ R⁺⁻)⁻¹     # Eq. (24): "T₀₁(I − R₂₁R₀₁)⁻¹"
    T21_inv = t⁺⁺ · (E − R⁺⁻ r⁻⁺)⁻¹     # Eq. (23): "T₂₁(I − R₀₁R₂₁)⁻¹"

Then the composite operators are updated **in place** following the
adding equations (lowercase = added layer; uppercase = composite):

    J₀⁻ ← J₀⁻ + T01_inv · (r⁻⁺ J₀⁺ + j₀⁻)        Eq. (28): J₀₂⁻
    R⁻⁺ ← R⁻⁺ + T01_inv · r⁻⁺ · T⁺⁺              Eq. (24): R₂₀
    T⁻⁻ ← T01_inv · t⁻⁻                          Eq. (25): T₀₂

    J₀⁺ ← j₀⁺ + T21_inv · (J₀⁺ + R⁺⁻ j₀⁻)        Eq. (27): J₂₀⁺
    T⁺⁺ ← T21_inv · T⁺⁺                          Eq. (23): T₂₀
    R⁺⁻ ← r⁺⁻ + T21_inv · R⁺⁻ · t⁻⁻              Eq. (26): R₀₂

The `(E − r⁻⁺ R⁺⁻)⁻¹` factor is the matrix geometric series capturing photon
paths that bounce between the two layers any number of times before
escaping. `batch_inv!` dispatches by array type to threaded BLAS (CPU),
CUBLAS strided (CUDA), or the portable KA LU kernel (Metal) — see
`docs/src/pages/concepts/07_architecture.md`.

This kernel runs once per atmospheric layer (after the TOA layer is copied
directly into the composite). The other three `ScatteringInterface_*`
methods drop terms that are zero by construction; this one is the full case.

# Concepts page
See [The MOM Solver § Adding / Interaction](../../docs/src/pages/concepts/04_mom_solver.md#adding--interaction)
for the prose walkthrough and the dispatch table that selects between the
four `ScatteringInterface_*` cases.
"""
function interaction_helper!(::ScatteringInterface_11, SFI,
                                composite_layer::CompositeLayer{FT},
                                added_layer::AddedLayer{FT},
                                I_static::AbstractArray{FT2}) where {FT<:Real,FT2}

    (; r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺, j₀⁺, j₀⁻, j₀_by_src,
       temp1, temp2, temp1_ptr, temp2_ptr) = added_layer     #these are aliases to the respective struct elements
    (; R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻, J₀⁺, J₀⁻, J₀_by_src) = composite_layer #these are aliases to the respective struct elements
    # Phase-0 in-place scratch: the T01_inv/T21_inv products live in mA so
    # temp1/temp2 remain free GEMM scratch on BOTH the fused and fallback
    # paths (the fallback needs them for the explicit inverse).
    mA, v1, v2 = _interaction_scratch(added_layer)

    # X₂₁ refers to added layer, X₁₀ to composite layer!

    # T01_inv = T⁻⁻·(E − r⁻⁺·R⁺⁻)⁻¹ — the adding factor reused by the J₀⁻, R⁻⁺ and
    # T⁻⁻ updates below (Sanghavi et al. 2014, Eqs. 23–28). The bare inverse is only
    # an intermediate, so on GPU compute it as a single fused LU + right-solve (no
    # explicit inverse, no separate GEMM); fall back to inverse-then-multiply for
    # large N / on CPU. `temp1` is repurposed as the T01_inv output.
    if _use_fused_solve(temp1)
        @timeit "interaction inv1 fused" ka_fused_solve!(mA, r⁻⁺, R⁺⁻, T⁻⁻,
                                                         KernelAbstractions.get_backend(mA))
    else
        _bmm!(temp1, r⁻⁺, R⁺⁻)
        temp2 .= I_static .- temp1
        @timeit "interaction inv1 bla" batch_inv!(temp1, temp2, temp1_ptr, temp2_ptr)
        _bmm!(mA, T⁻⁻, temp1)
    end
    T01_inv = mA

    # J₀₂⁻ = J₀₁⁻ + T₀₁(1-R₂₁R₀₁)⁻¹(R₂₁J₁₀⁺+J₁₂⁻) — legacy solar slot
    _bmm!(v1, r⁻⁺, J₀⁺); v1 .= v1 .+ j₀⁻
    _bmm!(v2, T01_inv, v1); J₀⁻ .= J₀⁻ .+ v2
    # Per-source J₀⁻ slots (same formula, T01_inv reused; uses pre-mutation r⁻⁺)
    for (key, slot) in pairs(j₀_by_src)
        cslot = J₀_by_src[key]
        _bmm!(v1, r⁻⁺, cslot.J₀⁺); v1 .= v1 .+ slot.j₀⁻
        _bmm!(v2, T01_inv, v1);    cslot.J₀⁻ .= cslot.J₀⁻ .+ v2
    end

    # R₂₀ = R₁₀ + T₀₁(I-R₂₁R₀₁)⁻¹ R₂₁T₁₀
    _bmm!(temp1, T01_inv, r⁻⁺)
    _bmm!(temp2, temp1, T⁺⁺); R⁻⁺ .= R⁻⁺ .+ temp2

    # T₀₂ = T₀₁(1-R₂₁R₀₁)⁻¹T₁₂
    _bmm!(temp1, T01_inv, t⁻⁻); T⁻⁻ .= temp1

    # Repeating for mirror-reflected directions

    # T21_inv = t⁺⁺·(E − R⁺⁻·r⁻⁺)⁻¹ — reused by the J₀⁺, T⁺⁺ and R⁺⁻ updates below.
    # Same fused LU + right-solve as inv1 (temp1 reused as the T21_inv output).
    if _use_fused_solve(temp1)
        @timeit "interaction inv2 fused" ka_fused_solve!(mA, R⁺⁻, r⁻⁺, t⁺⁺,
                                                         KernelAbstractions.get_backend(mA))
    else
        _bmm!(temp1, R⁺⁻, r⁻⁺)
        temp2 .= I_static .- temp1
        @timeit "interaction inv2" batch_inv!(temp1, temp2, temp1_ptr, temp2_ptr)
        _bmm!(mA, t⁺⁺, temp1)
    end
    T21_inv = mA

    # J₂₀⁺ = J₂₁⁺ + T₂₁(I-R₀₁R₂₁)⁻¹(J₁₀ + R₀₁J₁₂⁻ ) — legacy solar slot
    _bmm!(v1, R⁺⁻, j₀⁻); v1 .= J₀⁺ .+ v1
    _bmm!(v2, T21_inv, v1); J₀⁺ .= j₀⁺ .+ v2
    # Per-source J₀⁺ slots (same formula, T21_inv reused; uses pre-mutation R⁺⁻)
    for (key, slot) in pairs(j₀_by_src)
        cslot = J₀_by_src[key]
        _bmm!(v1, R⁺⁻, slot.j₀⁻); v1 .= cslot.J₀⁺ .+ v1
        _bmm!(v2, T21_inv, v1);   cslot.J₀⁺ .= slot.j₀⁺ .+ v2
    end

    # T₂₀ = T₂₁(I-R₀₁R₂₁)⁻¹T₁₀
    _bmm!(temp1, T21_inv, T⁺⁺); T⁺⁺ .= temp1

    # R₀₂ = R₁₂ + T₂₁(1-R₀₁R₂₁)⁻¹R₀₁T₁₂
    _bmm!(temp1, T21_inv, R⁺⁻)
    _bmm!(temp2, temp1, t⁻⁻); R⁺⁻ .= r⁺⁻ .+ temp2
end

"""
    interaction!(scattering_interface, SFI, composite_layer, added_layer, I_static)

Combine the accumulated [`CompositeLayer`](@ref) (from above) with a newly
doubled [`AddedLayer`](@ref) (below) using the adding equations.

Dispatches to the appropriate [`interaction_helper!`](@ref) based on the
`scattering_interface` type (`ScatteringInterface_00`, `_01`, `_10`, or
`_11`), then issues a GPU synchronisation barrier.
"""
function interaction!(scattering_interface::AbstractScatteringInterface, SFI,
                        composite_layer::CompositeLayer{FT}, 
                        added_layer::AddedLayer{FT},
                        I_static::AbstractArray{FT2}) where {FT<:Real,FT2}

    interaction_helper!(scattering_interface, SFI, composite_layer, added_layer, I_static)
    synchronize_if_gpu()
end
