# Scalar optical derivatives shared by dense and factored phase assembly.
"Moment-independent, δ-M-scaled aerosol quantities for one layer."
struct LinAerosolInvariant{T1,T2}
    τ::T1
    ϖ::T1
    τ̇::T2
    ϖ̇::T2
end

"""
    NativeLayerSelection

Internal component-local decomposition of an active native atmospheric
Jacobian basis. `include_pressure` controls the shared leading pressure
column; `aerosol_columns[iaer]` contains indices in that aerosol's historical
seven-column block; `gas_columns` contains indices in the flattened native gas
block. Splitting the selection before mixing prevents inactive columns from
entering the combined phase-Jacobian tensors.
"""
struct NativeLayerSelection
    include_pressure::Bool
    aerosol_columns::Vector{Vector{Int}}
    gas_columns::Vector{Int}
end

function _native_layer_selection(layout::ActiveParameterLayout,
                                 n_aerosol::Int, n_gas::Int)
    columns = native_layer_columns(layout)
    selected_pressure = 1 in columns
    aerosol_columns = Vector{Vector{Int}}(undef, n_aerosol)
    reconstructed = Int[]
    selected_pressure && push!(reconstructed, 1)
    for iaer in 1:n_aerosol
        native = (2 + 7 * (iaer - 1)):(1 + 7 * iaer)
        local_columns = [column - first(native) + 1
                         for column in columns if column in native]
        aerosol_columns[iaer] = local_columns
        append!(reconstructed, first(native) .+ local_columns .- 1)
    end
    gas_native = (2 + 7 * n_aerosol):(1 + 7 * n_aerosol + n_gas)
    gas_columns = [column - first(gas_native) + 1
                   for column in columns if column in gas_native]
    append!(reconstructed, first(gas_native) .+ gas_columns .- 1)
    reconstructed == columns || throw(ArgumentError(
        "active native layer columns must follow pressure, aerosol-component, " *
        "then gas order; got $columns"))
    return NativeLayerSelection(selected_pressure, aerosol_columns, gas_columns)
end

"""
Moment-invariant cache used by the linearized Fourier loop. When `selection`
is non-`nothing`, aerosol and gas tangent arrays already carry only the active
component-local columns; all forward arrays remain complete.
"""
struct LinMInvariantCache
    rayl_τ_dev::Vector{Vector}
    aerosol::Vector{Vector{Vector}}
    gas::Vector{Vector}
    lin_gas::Vector{Vector}
    selection::Union{Nothing,NativeLayerSelection}
end

_same_greek(a, b) = all(name -> getfield(a, name) == getfield(b, name),
                         fieldnames(typeof(a)))

function _createAero_invariant(τAer, aerosol_optics, τ̇Aer,
                               lin_aerosol_optics, arr_type,
                               columns::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    (; fᵗ, ω̃) = aerosol_optics
    (; ḟᵗ, ω̃̇) = lin_aerosol_optics
    n = size(τAer, 1)
    τ̇Aer = _to_device(arr_type, collect(τ̇Aer'))
    ω̃ = ω̃ isa Number ? arr_type(fill(ω̃, n)) : _to_device(arr_type, ω̃)
    fᵗ = fᵗ isa Number ? arr_type(fill(fᵗ, n)) : _to_device(arr_type, fᵗ)
    ω̃̇_block = _lift_mie_param_to_n_x_4(ω̃̇, n, arr_type)
    ḟᵗ_block = _lift_mie_param_to_n_x_4(ḟᵗ, n, arr_type)
    # Truncation is part of the optical-property chain rule, not a missing RT
    # tangent. With f = fᵗ, ω = ω̃ and a = 1-fω, the transformed aerosol has
    # τ* = aτ and ω* = (1-f)ω/a (Sanghavi et al. 2014, Appendix A, Eq. A.3).
    # Therefore dτ* = a dτ - τ(f dω + ω df), and
    # dω* = [(1-f)dω - ω(1-ω)df]/a². These are evaluated below for every
    # microphysical direction; cf. the mixed-layer partials C.28–C.35.
    # The truncated phase derivative is supplied separately by Scattering.
    # Its normalization also includes df. The paper's truncation β_i is fᵗ
    # here, distinct from the degree-indexed Greek coefficient β_l.
    fω = fᵗ .* ω̃
    τ_mod = (1 .- fω) .* τAer
    ϖ_mod = (1 .- fᵗ) .* ω̃ ./ (1 .- fω)
    τ̇_mod = arr_type(zeros(eltype(τAer), n, 7))
    ϖ̇_mod = arr_type(zeros(eltype(τAer), n, 7))
    τ̇_mod[:,1] .= (1 .- fω) .* τ̇Aer[:,1]
    tmp = fᵗ .* ω̃̇_block .+ ω̃ .* ḟᵗ_block
    τ̇_mod[:,2:5] .= (1 .- fω) .* τ̇Aer[:,2:5] .- tmp .* τAer
    ϖ̇_mod[:,2:5] .= (ω̃̇_block .* (1 .- fᵗ) .-
        ḟᵗ_block .* (ω̃ .* (1 .- ω̃))) ./ (1 .- fω).^2
    τ̇_mod[:,6:7] .= (1 .- fω) .* τ̇Aer[:,6:7]
    if columns !== nothing
        τ̇_mod = τ̇_mod[:, columns]
        ϖ̇_mod = ϖ̇_mod[:, columns]
    end
    return LinAerosolInvariant(τ_mod, ϖ_mod, τ̇_mod, ϖ̇_mod)
end

"""
    build_m_invariant_cache_lin(iBand, model, lin_model;
                                active_layout=nothing)

Construct Fourier-independent Rayleigh, aerosol, and gas optical properties.
With an active layout, split its native atmospheric columns by component and
cache only those aerosol/gas tangents. This is the earliest safe selection
point: the forward mixing weights are already known, but no Fourier-dependent
phase tensor or MOM operator workspace has been allocated.
"""
function build_m_invariant_cache_lin(iBand, model, lin_model;
                                     active_layout::Union{Nothing,ActiveParameterLayout}=nothing)
    bands = iBand isa Integer ? (iBand,) : iBand
    (; τ_rayl, τ_aer, τ_abs, aerosol_optics) = model
    (; τ̇_aer, τ̇_abs, lin_aerosol_optics) = lin_model
    arr_type = CoreRT.array_type(model)
    nZ = size(τ_rayl[1], 2)
    nAero = size(τ_aer[first(bands)], 1)
    nGas = size(τ̇_abs[first(bands)], 1)
    selection = active_layout === nothing ? nothing :
        _native_layer_selection(active_layout, nAero, nGas)
    if (selection === nothing || selection.include_pressure) &&
            any(isnothing, (lin_model.τ̇_rayl_psurf, lin_model.τ̇_aer_psurf,
                            lin_model.τ̇_abs_psurf))
        throw(ArgumentError("Surface-pressure derivatives are unavailable; rebuild " *
            "with compute_pressure_jacobians=true or use a plan without pressure"))
    end
    rayl = Vector{Vector}(undef, length(bands))
    aeros = Vector{Vector{Vector}}(undef, length(bands))
    gas = Vector{Vector}(undef, length(bands))
    lin_gas = Vector{Vector}(undef, length(bands))
    for (iBi, iB) in enumerate(bands)
        rayl[iBi] = [_to_device(arr_type, τ_rayl[iB][:,iz]) for iz in 1:nZ]
        aeros[iBi] = [[_createAero_invariant(
            _to_device(arr_type, τ_aer[iB][iaer,:,iz]), aerosol_optics[iB][iaer],
            τ̇_aer[iB][iaer,:,:,iz], lin_aerosol_optics[iB][iaer], arr_type,
            selection === nothing ? nothing : selection.aerosol_columns[iaer])
            for iz in 1:nZ] for iaer in 1:nAero]
        gas[iBi] = [CoreAbsorptionOpticalProperties(
            _to_device(arr_type, τ_abs[iB][:,iz])) for iz in 1:nZ]
        gas_columns = selection === nothing ? axes(τ̇_abs[iB], 1) :
                      selection.gas_columns
        lin_gas[iBi] = [CoreAbsorptionOpticalPropertiesLin(
            _to_device(arr_type, collect(τ̇_abs[iB][gas_columns,:,iz]'))) for iz in 1:nZ]
    end
    return LinMInvariantCache(rayl, aeros, gas, lin_gas, selection)
end
