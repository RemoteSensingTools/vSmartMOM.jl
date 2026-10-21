#=
 
This file implements rt_kernel!, which performs the core RT routines (elemental, doubling, interaction)
 
=#

"""
    rt_kernel!(RS_type, pol_type, SFI, added_layer, added_layer_lin, 
               composite_layer, composite_layer_lin, 
               computed_layer_properties, computed_layer_properties_lin,
               scattering_interface, τ_sum, τ̇_sum, m, quad_points, 
               I_static, architecture, qp_μN, iz)

Core RT kernel for a single atmospheric layer `iz` in the linearized mode.

The direction representation is selected by dispatch on `jacobian`:

1. Build finite-thickness elemental operators and their supplied tangents.
2. Double these complete matrix directions to the layer's optical thickness.
3. For `LocalOpticalJacobian`, contract the local basis to retrieval columns
   and append the above-layer direct-beam attenuation derivative.
4. Seed the TOA composite, or add this layer to the composite above it.

The phase chain rule acts on complete matrix directions before doubling.
After doubling, matrix products couple angular indices: an elementwise
`∂r/∂Z .* dZ` contraction is invalid. A scalar-weighted combination of
**complete doubled tangent matrices**, however, is valid by linearity of the
fixed-forward-state tangent map; see S2014 (C.25)–(C.26).

# Dispatch
- `noRS`: Elastic scattering only (no Raman).
- Future: `RRS`, `VRS` for inelastic scattering (not yet linearized).

See Sanghavi, Davis & Eldering (2014, JQSRT 133:412–433) for the full
forward (Eqs. 19–32) and linearization (App. C) framework. The elemental
kernel uses the *exact* finite-δ formulas of Fell (1997) Eqs. 1.52–1.56,
restated as Sanghavi & Frankenberg (2023, JQSRT 311:108791) Eqs. (10)–(11),
not the linear S2014 Eqs. (19)–(20) limit. See `docs/src/pages/concepts/04_mom_solver.md`
§ Elemental and `docs/src/pages/concepts/06_linearization.md`.
"""
# One orchestration path: local layer propagation, then atmospheric adding.
function rt_kernel!(RS_type::noRS{FT}, pol_type, SFI,
        added_layer, added_layer_lin, composite_layer, composite_layer_lin,
        optics::CoreScatteringOpticalProperties,
        jacobian::AbstractOpticalPropertiesLin, scattering_interface,
        τ_sum, τ̇_sum, m, quad_points, I_static, architecture, qp_μN, iz;
        dτ_max_threshold::Union{Nothing,Real}=nothing,
        dτ_min_floor::Union{Nothing,Real}=nothing,
        local_workspace=nothing) where {FT}
    # The discrete doubling count is chosen from the forward state and held
    # fixed for tangent propagation, as in the existing analytic path.
    dτ, nd = get_dtau_ndoubl(optics, quad_points;
        dτ_max_threshold, dτ_min_floor)
    AT = array_type(architecture)
    build_doubled_layer_lin!(pol_type, SFI, _to_device(AT,τ_sum),
        _to_device(AT,τ̇_sum), dτ, _to_device(AT,RS_type.F₀), optics, jacobian,
        m, nd, quad_points, added_layer, added_layer_lin, architecture,
        I_static, local_workspace)
    # Both dispatches return complete physical-coordinate operator/source
    # tangents. No elementwise phase chain rule is applied after doubling.
    if iz == 1
        seed_composite_from_added!(composite_layer, composite_layer_lin,
            added_layer, added_layer_lin)
    else
        @timeit "interaction" interaction!(scattering_interface, SFI,
            composite_layer, composite_layer_lin, added_layer, added_layer_lin,
            I_static)
    end
    return nothing
end
