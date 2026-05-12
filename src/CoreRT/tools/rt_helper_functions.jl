#=

This file contains helper functions that are used throughout the vSmartMOM module

=#

"""
    get_scattering_interface(scattering_interface, scatter, iz)

Determine the `ScatteringInterface` type for layer `iz` based on previous interface and current scattering.

- **iz == 1**: TOA → `ScatteringInterface_11` (scattering) or `00` (non-scattering).
- **iz > 1**: Toggles between 00↔01 or 10↔11 depending on whether the composite layer scatters.
"""
@inline function get_scattering_interface(scattering_interface, scatter, iz)

    # First layer (TOA)
    if (iz == 1)

        # If scattering, 11. If non-scattering, 00. 
        scattering_interface = scatter ? ScatteringInterface_11() : ScatteringInterface_00()
    
    # Not the first layer (not TOA)
    else
        
        # If kn was 1, then toggle between 0/0 and 0/1 
        # Else, toggle between 1/0 and 1/1
        scattering_interface = (scattering_interface isa ScatteringInterface_00) ? 
                                    (!scatter ? ScatteringInterface_00() : ScatteringInterface_01()) : 
                                    (!scatter ? ScatteringInterface_10() : ScatteringInterface_11())
    end
    return scattering_interface
end

"""
    doubling_number(dτ_max, τ_end)

Compute elemental optical depth and doubling count for the adding-doubling method.

Finds the smallest `ndoubl` such that ``d\\tau = \\tau_{\\text{end}} / 2^{n_d} \\leq d\\tau_{\\text{max}}``.
Ensures the elemental layer is optically thin enough for accurate single-scattering.

# Returns
- `dτ`: Elemental optical depth.
- `ndoubl`: Number of doubling iterations.

See Sanghavi & Stephens (2013).
"""
@inline function doubling_number(dτ_max, τ_end)
    FT = eltype(dτ_max)
    if τ_end <= dτ_max
        return τ_end, 0
    else
        q1 = log10(FT(2))
        q2 = log10(dτ_max)
        q3 = log10(τ_end)
        tlimit = (q3 - q2) / q1
        nlimit = floor(Int, tlimit)
        diff = tlimit - nlimit
        if diff < eps(FT)
            return dτ_max, nlimit
        else
            ndoubl = nlimit + 1       
            x = q3 - q1 * ndoubl
            dτ = FT(10)^x
            return dτ, ndoubl
        end 
    end
end

"""
    nearest_point(f_array, f)

Index of `f_array` whose value is closest to `f`. Used to map VZA to quadrature points.
"""
@inline nearest_point(f_array, f) = argmin(abs.(f_array.-f))

"""
    get_indices(iμ, pol_type)

Stokes-scaled indices for quadrature point `iμ`. Returns `(st_iμ, istart, iend)` where
`istart:iend` spans the Stokes components for that direction.
"""
@inline function get_indices(iμ::Integer, pol_type::AbstractPolarizationType) 
    st_iμ = (iμ - 1) * pol_type.n
    istart = st_iμ + 1
    iend   = st_iμ + pol_type.n
    return st_iμ, istart, iend
end

"""Allocate zero matrix for RT reflection/transmission: `[nμ, nμ, nSpec]`."""
@inline default_matrix(FT, arr_type, dims, nSpec) = arr_type(zeros(FT, (dims[1], dims[2], nSpec)))
"Default matrix in ieRT calculation (zeros)"
@inline default_matrix_ie(FT, arr_type, dims, nSpec, nRaman) = arr_type(zeros(FT, (dims[1], dims[2], nSpec, nRaman)))

"Default J matrix in RT calculation (zeros)"
@inline default_J_matrix(FT, arr_type, dims, nSpec) = arr_type(zeros(FT, (dims[1], 1, nSpec)))
"Default J matrix in ieRT calculation (zeros)"
@inline default_J_matrix_ie(FT, arr_type, dims, nSpec, nRaman) = arr_type(zeros(FT, (dims[1], 1, nSpec, nRaman)))

"Default matrix in RT calculation (zeros) — multi-sensor variant"
@inline default_matrix(FT, arr_type, NSens, dims, nSpec) = [arr_type(zeros(FT, (dims[1], dims[2], nSpec))) for _ in 1:NSens]
"Default matrix in ieRT calculation (zeros) — multi-sensor variant"
@inline default_matrix_ie(FT, arr_type, NSens, dims, nSpec, nRaman) = [zeros(FT, (dims[1], dims[2], nSpec, nRaman)) for _ in 1:NSens]

"Default J matrix in RT calculation (zeros) — multi-sensor variant"
@inline default_J_matrix(FT, arr_type, NSens, dims, nSpec) = [arr_type(zeros(FT, (dims[1], 1, nSpec))) for _ in 1:NSens]
"Default J matrix in ieRT calculation (zeros) — multi-sensor variant"
@inline default_J_matrix_ie(FT, arr_type, NSens, dims, nSpec, nRaman) = [zeros(FT, (dims[1], 1, nSpec, nRaman)) for _ in 1:NSens]

##### Only for testing, random matrices:
"Default matrix in RT calculation (random)"
default_matrix_rand(FT, arr_type, dims, nSpec)   = arr_type(randn(FT, tuple(dims[1], dims[2], nSpec)))

"Default J matrix in RT calculation (random)"
default_J_matrix_rand(FT, arr_type, dims, nSpec) = arr_type(randn(FT, tuple(dims[1], 1, nSpec)))


"""
    make_added_layer(RS_type, FT, arr_type, dims, nSpec; prepared_sources=NoSource())

Construct an `AddedLayer` with zero-initialized R, T, J₀ matrices for elastic RT.

When `prepared_sources` contains any source that declares a non-`nothing`
[`source_key`](@ref) (e.g. `:thermal` for `PreparedThermalEmission`), a
matching [`SourceSlot`](@ref) is allocated under `j₀_by_src[key]`. Solar
SFI and surface SIF use the legacy `j₀⁺/j₀⁻` fields, so they don't get a
separate slot. Empty prepared-sources → empty NT, bit-equal to pre-A.2a.
"""
function make_added_layer(RS_type::Union{noRS, noRS_plus}, FT, arr_type, dims, nSpec;
                          prepared_sources::AbstractSource = NoSource())
    t1 = default_matrix(FT, arr_type, dims, nSpec)
    t2 = default_matrix(FT, arr_type, dims, nSpec)
    t1_ptr = batched_pointer_cache(t1)
    t2_ptr = batched_pointer_cache(t2)
    j₀_by_src = _make_added_source_slots(prepared_sources, FT, arr_type, dims, nSpec)
    return AddedLayer(
        r⁻⁺ = default_matrix(FT, arr_type, dims, nSpec),
        t⁺⁺ = default_matrix(FT, arr_type, dims, nSpec),
        r⁺⁻ = default_matrix(FT, arr_type, dims, nSpec),
        t⁻⁻ = default_matrix(FT, arr_type, dims, nSpec),
        j₀⁺ = default_J_matrix(FT, arr_type, dims, nSpec),
        j₀⁻ = default_J_matrix(FT, arr_type, dims, nSpec),
        temp1 = t1, temp2 = t2,
        temp1_ptr = t1_ptr, temp2_ptr = t2_ptr,
        dbl_gp_refl = default_matrix(FT, arr_type, dims, nSpec),
        dbl_j₁⁺ = default_J_matrix(FT, arr_type, dims, nSpec),
        dbl_j₁⁻ = default_J_matrix(FT, arr_type, dims, nSpec),
        j₀_by_src = j₀_by_src,
    )
end

# ============================================================================
# v0.7 Phase A.2a — per-source slot allocators for the elastic noRS path.
#
# Iterate the prepared-sources tuple, collect distinct non-nothing source_key
# values, and build a NamedTuple of SourceSlot (added-layer) /
# CompositeSourceSlot (composite-layer) carriers.
#
# Solar (`source_key === nothing`) is skipped — its j₀ lives in the legacy
# fields. NoSource is also skipped. Same for any source type that opts out by
# returning `nothing` from `source_key`.
# ============================================================================
_collect_source_keys(::NoSource) = Symbol[]
function _collect_source_keys(src::AbstractPreparedSource)
    k = source_key(src)
    return k === nothing ? Symbol[] : Symbol[k]
end
function _collect_source_keys(s::SourceSet)
    keys = Symbol[]
    for member in s.sources
        k = source_key(member)
        k === nothing && continue
        k in keys && continue
        push!(keys, k)
    end
    return keys
end
# Raw (unprepared) sources should never reach the allocator — but if they do,
# fall back to the prepared-source contract by reading the trait.
_collect_source_keys(src::AbstractSource) = Symbol[]

function _make_added_source_slots(prepared_sources::AbstractSource, FT, arr_type, dims, nSpec)
    keys = _collect_source_keys(prepared_sources)
    isempty(keys) && return NamedTuple()
    nspec_vec = ones(FT, nSpec)
    slots = map(keys) do key
        SourceSlot(
            j₀⁺     = default_J_matrix(FT, arr_type, dims, nSpec),
            j₀⁻     = default_J_matrix(FT, arr_type, dims, nSpec),
            dbl_j₁⁺ = default_J_matrix(FT, arr_type, dims, nSpec),
            dbl_j₁⁻ = default_J_matrix(FT, arr_type, dims, nSpec),
            # For thermal (and every other current per-source slot type)
            # expk = ones; squaring during the doubling loop is a no-op. We
            # default to `ones` here so the doubling kernel can treat every
            # slot uniformly; future sources can override `source_expk_init`
            # and a separate seed step before the doubling loop will reseed.
            expk    = arr_type(copy(nspec_vec)),
        )
    end
    return NamedTuple{Tuple(keys)}(Tuple(slots))
end

function _make_composite_source_slots(prepared_sources::AbstractSource, FT, arr_type, dims, nSpec)
    keys = _collect_source_keys(prepared_sources)
    isempty(keys) && return NamedTuple()
    slots = map(keys) do _key
        CompositeSourceSlot(
            J₀⁺ = default_J_matrix(FT, arr_type, dims, nSpec),
            J₀⁻ = default_J_matrix(FT, arr_type, dims, nSpec),
        )
    end
    return NamedTuple{Tuple(keys)}(Tuple(slots))
end

"""Construct an `AddedLayerRS` with inelastic (Raman) matrices for Raman scattering."""
make_added_layer(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus}, FT, arr_type, dims, nSpec)  = AddedLayerRS(
                                                default_matrix(FT, arr_type, dims, nSpec), 
                                                default_matrix(FT, arr_type, dims, nSpec), 
                                                default_matrix(FT, arr_type, dims, nSpec),
                                                default_matrix(FT, arr_type, dims, nSpec),
                                                default_J_matrix(FT, arr_type, dims, nSpec),
                                                default_J_matrix(FT, arr_type, dims, nSpec),
                                                default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman), 
                                                default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman), 
                                                default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                default_J_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                default_J_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman)
                                                )
                                                         

"Make a random added layer, supplying all random matrices"
function make_added_layer_rand(RS_type::Union{noRS, noRS_plus}, FT, arr_type, dims, nSpec)
    t1 = default_matrix_rand(FT, arr_type, dims, nSpec)
    t2 = default_matrix_rand(FT, arr_type, dims, nSpec)
    return AddedLayer(
        r⁻⁺ = default_matrix_rand(FT, arr_type, dims, nSpec),
        t⁺⁺ = default_matrix_rand(FT, arr_type, dims, nSpec),
        r⁺⁻ = default_matrix_rand(FT, arr_type, dims, nSpec),
        t⁻⁻ = default_matrix_rand(FT, arr_type, dims, nSpec),
        j₀⁺ = default_J_matrix_rand(FT, arr_type, dims, nSpec),
        j₀⁻ = default_J_matrix_rand(FT, arr_type, dims, nSpec),
        temp1 = t1,
        temp2 = t2,
        temp1_ptr = batched_pointer_cache(t1),
        temp2_ptr = batched_pointer_cache(t2),
        dbl_gp_refl = default_matrix_rand(FT, arr_type, dims, nSpec),
        dbl_j₁⁺ = default_J_matrix_rand(FT, arr_type, dims, nSpec),
        dbl_j₁⁻ = default_J_matrix_rand(FT, arr_type, dims, nSpec),
    )
end
                                                         
"""Construct a `CompositeLayer` with zero-initialized R, T, J₀ for elastic RT.

When `prepared_sources` declares any non-`nothing` [`source_key`](@ref),
matching [`CompositeSourceSlot`](@ref)s are allocated under `J₀_by_src`.
"""
make_composite_layer(RS_type::Union{noRS, noRS_plus},
    FT, arr_type, dims, nSpec;
    prepared_sources::AbstractSource = NoSource()) =
    CompositeLayer(
        R⁻⁺ = default_matrix(FT, arr_type, dims, nSpec),
        R⁺⁻ = default_matrix(FT, arr_type, dims, nSpec),
        T⁺⁺ = default_matrix(FT, arr_type, dims, nSpec),
        T⁻⁻ = default_matrix(FT, arr_type, dims, nSpec),
        J₀⁺ = default_J_matrix(FT, arr_type, dims, nSpec),
        J₀⁻ = default_J_matrix(FT, arr_type, dims, nSpec),
        J₀_by_src = _make_composite_source_slots(prepared_sources, FT, arr_type, dims, nSpec),
    )
"""Construct a `CompositeLayerRS` with inelastic matrices for Raman scattering."""
make_composite_layer(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus},
    FT, arr_type, dims, nSpec) = CompositeLayerRS(
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec), 
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec),
                                                        default_J_matrix(FT, arr_type, dims, nSpec),
                                                        default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman), 
                                                        default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman), 
                                                        default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                        default_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                        default_J_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman),
                                                        default_J_matrix_ie(FT, arr_type, dims, nSpec, RS_type.n_Raman)
                                                        )
                                                    
"Make a composite layer, supplying all default matrices"
make_composite_layer(RS_type::Union{noRS, noRS_plus}, 
    FT, arr_type, NSens, dims, nSpec) = 
    CompositeLayerMS(default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec))
"Make a composite layer, supplying all default matrices"
make_composite_layer(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus},
    FT, arr_type, NSens, dims, nSpec) = 
    CompositeLayerMSRS(default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman), 
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman), 
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_J_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_J_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec), 
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_J_matrix(FT, arr_type, NSens, dims, nSpec),
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman), 
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman), 
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_J_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman),
                    default_J_matrix_ie(FT, arr_type, NSens, dims, nSpec, RS_type.n_Raman)
                    )
"""
    get_layer_properties(computed_atmospheric_properties, iz, arr_type)

Extract layer `iz` optical properties from `ComputedAtmosphereProperties`.

Returns a `ComputedLayerProperties` with τ, ϖ, Z⁺⁺, Z⁻⁺, dτ, ndoubl, expk, scatter,
τ_sum, and scattering_interface for the specified layer.
"""
function get_layer_properties(computed_atmospheric_properties::ComputedAtmosphereProperties, iz, arr_type)
     (; τ_λ_all, ϖ_λ_all, τ_all, ϖ_all, Z⁺⁺_all, Z⁻⁺_all , dτ_max_all, dτ_all, ndoubl_all, dτ_λ_all, expk_all, scatter_all, τ_sum_all, fscattRayl_all,  scattering_interfaces_all) = computed_atmospheric_properties

    τ_λ = arr_type(τ_λ_all[:, iz])
    ϖ_λ = arr_type(ϖ_λ_all[:, iz])
    τ   = τ_all[iz]
    ϖ   = ϖ_all[iz]
    Z⁺⁺ = arr_type(Z⁺⁺_all[:,:,iz])
    Z⁻⁺ = arr_type(Z⁻⁺_all[:,:,iz])

    dτ_max = dτ_max_all[iz]
    dτ     = dτ_all[iz]
    ndoubl = ndoubl_all[iz]
    dτ_λ   = arr_type(dτ_λ_all[:, iz])
    expk   = arr_type(expk_all[:, iz])
    scatter = scatter_all[iz]
    τ_sum = arr_type(τ_sum_all[:,iz])
    scattering_interface = scattering_interfaces_all[iz]
    fscattRayl = fscattRayl_all[iz]
    #ϖ_Cabannes = ϖ_Cabannes_all[iz]
    # τ * ϖ should remain constant even though they individually change over wavelength
    # @assert all(i -> (i ≈ τ * ϖ), τ_λ .* ϖ_λ)

    return ComputedLayerProperties(τ_λ, ϖ_λ, τ, ϖ, 
        Z⁺⁺, Z⁻⁺, 
        dτ_max, dτ, 
        ndoubl, 
        dτ_λ, 
        expk, 
        scatter, 
        τ_sum, 
        fscattRayl, 
        scattering_interface)
end
