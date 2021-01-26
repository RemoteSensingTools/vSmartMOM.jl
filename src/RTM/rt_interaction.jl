# No scattering in either the added layer or the composite layer.
function rt_interaction_helper!(::ScatteringInterface_00, 
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                ::AbstractArray{FT}) where {FT}

    composite_layer.T⁻⁻[:] = added_layer.t⁻⁻ ⊠ composite_layer.T⁻⁻
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ ⊠ composite_layer.T⁺⁺
end

# No scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer, added to bottom of the composite layer.
# Produces a new, scattering composite layer.
function rt_interaction_helper!(::ScatteringInterface_01, 
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                ::AbstractArray{FT}) where {FT}
    composite_layer.R⁻⁺[:] = composite_layer.T⁻⁻ * added_layer.r⁻⁺ * composite_layer.T⁺⁺
    composite_layer.R⁺⁻[:] = added_layer.r⁺⁻
    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ * composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ * added_layer.t⁻⁻
end

# Scattering in inhomogeneous composite layer.
# no scattering in homogeneous layer which is 
# added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function rt_interaction_helper!(::ScatteringInterface_10, 
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT}) where {FT}

    composite_layer.T⁺⁺[:] = added_layer.t⁺⁺ * composite_layer.T⁺⁺
    composite_layer.T⁻⁻[:] = composite_layer.T⁻⁻ * added_layer.t⁻⁻
    composite_layer.R⁺⁻[:] = added_layer.t⁺⁺ * composite_layer.R⁺⁻ * added_layer.t⁻⁻
end

# Scattering in inhomogeneous composite layer.
# Scattering in homogeneous layer which is added to the bottom of the composite layer.
# Produces a new, scattering composite layer.
function rt_interaction_helper!(::ScatteringInterface_11, 
                                composite_layer::CompositeLayer{FT}, 
                                added_layer::AddedLayer{FT}, 
                                I_static::AbstractArray{FT}) where {FT}
    
    @unpack r⁺⁻, r⁻⁺, t⁻⁻, t⁺⁺ = added_layer
    @unpack R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻ = composite_layer

    # Used to store `(I - R⁺⁻ * r⁻⁺)⁻¹ * T⁺⁺`
    tmp_inv = similar(t⁺⁺)

    # Compute and store `(I - R⁺⁻ * r⁻⁺)⁻¹ * T⁺⁺`
    @timeit "interaction inv1" batch_solve!(tmp_inv, I_static .- R⁺⁻ ⊠ r⁻⁺, T⁺⁺)

    composite_layer.R⁻⁺[:] = R⁻⁺ + (T⁻⁻ ⊠ r⁻⁺ ⊠ tmp_inv)
    composite_layer.T⁺⁺[:] = t⁺⁺ ⊠ tmp_inv

    # Repeating for mirror-reflected directions

    # Compute and store `(I - r⁻⁺ * R⁺⁻)⁻¹ * t⁻⁻`
    @timeit "interaction inv2" batch_solve!(tmp_inv, I_static .- r⁻⁺ ⊠ R⁺⁻, t⁻⁻)

    composite_layer.R⁺⁻[:] = r⁺⁻ + t⁺⁺ ⊠ R⁺⁻ ⊠ tmp_inv
    composite_layer.T⁻⁻[:] = T⁺⁺ ⊠ tmp_inv
end


function rt_interaction!(scattering_interface::AbstractScatteringInterface, 
                        composite_layer::CompositeLayer{FT}, 
                        added_layer::AddedLayer{FT},
                        I_static::AbstractArray{FT}) where {FT}

    rt_interaction_helper!(scattering_interface, composite_layer, added_layer, I_static)
    synchronize()

end