"Simulates the full atmosphere from n distinct homogeneous layers"
# function rt_interaction!(kn, R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
#     aux2 = similar(R⁺⁻);# I_static = one(similar(R⁺⁻))
#     aux3 = similar(R⁺⁻); 
#     aux1 = similar(R⁺⁻)                     
#     # ToDo: Important output from this routine is R⁻⁺, R⁺⁻, T⁺⁺, T⁻⁻ (can be renamed to 𝐓⁻⁻, etc later)
#     # Need to check with paper nomenclature. This is basically eqs. 23-28 in vSmartMOM)
#     Nquadn = size(r⁻⁺, 1)
#     # kn = 1: no scattering in either the added layer or composite layer.
#     # kn = 2: composite layer has no scattering but added layer does.
#     # kn = 3: composite layer has scattering but added layer does not.
#     # kn = 4: both composite layer and added layer have scattering.

#     # Create temporary matrices
#     I_static = one(similar(R⁺⁻))
#     # aux1 = similar(R⁺⁻)
#     # aux2 = similar(R⁺⁻)
#     # aux3 = similar(R⁺⁻)

#     if kn == 1
#         # No scattering in either the added layer or the composite layer.
#         T⁻⁻ = t⁻⁻ * T⁻⁻
#         T⁺⁺ = t⁺⁺ * T⁺⁺
        
#         return nothing 
#     elseif kn == 2
#         # No scattering in inhomogeneous composite layer.
#         # scattering in homogeneous layer which is added 
#         # to the bottom of the composite layer.
#         # Produces a new, scattering composite layer.
#         M1 = T⁻⁻
#         M2 = T⁺⁺
#         R⁻⁺[:] = M1 * r⁻⁺ * M2
#         R⁺⁻[:] = r⁺⁻
#         T⁺⁺[:] = t⁺⁺ * M2
#         T⁻⁻[:] = M1 * t⁻⁻
#         return nothing 
#     elseif kn == 3
#         # Scattering in inhomogeneous composite layer.
#         # no scattering in homogeneous layer which is 
#         # added to the bottom of the composite layer.
#         # Produces a new, scattering composite layer.
#         T⁺⁺[:] = t⁺⁺ * T⁺⁺
#         T⁻⁻[:] = T⁻⁻ * t⁻⁻
#         R⁺⁻[:] = t⁺⁺ * R⁺⁻ * t⁻⁻
#         return nothing 
#     elseif kn == 4
#         # Scattering in inhomogeneous composite layer.
#         # scattering in homogeneous layer which is 
#         # added to the bottom of the composite layer.
#         # Produces a new, scattering composite layer.

#         # M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺
#         mul!(aux1, R⁺⁻, r⁻⁺)        # R⁺⁻ * r⁻⁺
#         @. aux1 = I_static - aux1   # (I - R⁺⁻ * r⁻⁺)
#         ldiv!(aux2, qr!(aux1), T⁺⁺) # M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺

#         # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
#         mul!(aux1, r⁻⁺, aux2)   # r⁻⁺ * M1
#         mul!(aux3, T⁻⁻, aux1)   # T⁻⁻ * r⁻⁺ * M1
#         @. R⁻⁺ = R⁻⁺ + aux3     # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
        
#         # t_T⁺⁺ = t⁺⁺ * M1
#         mul!(T⁺⁺, t⁺⁺, aux2)

#         # Repeating for mirror-reflected directions

#         # M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻
#         mul!(aux1, r⁻⁺, R⁺⁻)        # r⁻⁺ * R⁺⁻
#         @. aux1 = I_static - aux1   # (I - r⁻⁺ * R⁺⁻)
#         ldiv!(aux2, qr!(aux1), t⁻⁻) # M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻

#         # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1
#         mul!(aux3, R⁺⁻, aux2)   # R⁺⁻ * M1
#         mul!(aux1, t⁺⁺, aux3)   # t⁺⁺ * R⁺⁻ * M1
#         @. R⁺⁻ = r⁺⁻ + aux1     # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1

#         # t_T⁻⁻ = T⁺⁺ * M1
#         mul!(T⁻⁻, T⁺⁺, aux2)
                 
#         return nothing 
#     end
# end

# function rt_interaction!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, aux1, aux2, aux3)
#     # M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺;aux1 = similar(R⁺⁻)
#         # aux2 = similar(R⁺⁻);#I_static = one(similar(R⁺⁻))
#         # aux3 = similar(R⁺⁻); 
#         # aux1 = similar(R⁺⁻)
#         I_static = one(similar(R⁺⁻))
#         mul!(aux1, R⁺⁻, r⁻⁺)        # R⁺⁻ * r⁻⁺
#         @. aux1 = I_static - aux1   # (I - R⁺⁻ * r⁻⁺)
#         ldiv!(aux2, qr!(aux1), T⁺⁺) # M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺

#         # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
#         mul!(aux1, r⁻⁺, aux2)   # r⁻⁺ * M1
#         mul!(aux3, T⁻⁻, aux1)   # T⁻⁻ * r⁻⁺ * M1
#         @. R⁻⁺ = R⁻⁺ + aux3     # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
        
#         # t_T⁺⁺ = t⁺⁺ * M1
#         mul!(T⁺⁺, t⁺⁺, aux2)

#         # Repeating for mirror-reflected directions

#         # M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻
#         mul!(aux1, r⁻⁺, R⁺⁻)        # r⁻⁺ * R⁺⁻
#         @. aux1 = I_static - aux1   # (I - r⁻⁺ * R⁺⁻)
#         ldiv!(aux2, qr!(aux1), t⁻⁻) # M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻

#         # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1
#         mul!(aux3, R⁺⁻, aux2)   # R⁺⁻ * M1
#         mul!(aux1, t⁺⁺, aux3)   # t⁺⁺ * R⁺⁻ * M1
#         @. R⁺⁻ = r⁺⁻ + aux1     # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1

#         # t_T⁻⁻ = T⁺⁺ * M1
#         mul!(T⁻⁻, T⁺⁺, aux2)
                 
#         return nothing 
# end


function rt_interaction_helper!(R⁻⁺::AbstractArray{FT,3}, 
                         T⁺⁺::AbstractArray{FT,3}, 
                         R⁺⁻::AbstractArray{FT,3}, 
                         T⁻⁻::AbstractArray{FT,3}, 
                         r⁻⁺::AbstractArray{FT,3}, 
                         t⁺⁺::AbstractArray{FT,3}, 
                         r⁺⁻::AbstractArray{FT,3}, 
                         t⁻⁻::AbstractArray{FT,3},
                         I_static::AbstractArray) where {FT}

    
    R⁻⁺_ = R⁻⁺ # repeat(R⁻⁺, 1, 1, 1)
    T⁺⁺_ = T⁺⁺ # repeat(T⁺⁺, 1, 1, 1)
    R⁺⁻_ = R⁺⁻ # repeat(R⁺⁻, 1, 1, 1)
    T⁻⁻_ = T⁻⁻ # repeat(T⁻⁻, 1, 1, 1)

    r⁻⁺_ = r⁻⁺
    t⁺⁺_ = t⁺⁺
    r⁺⁻_ = r⁺⁻
    t⁻⁻_ = t⁻⁻

    # r⁻⁺_ = repeat(r⁻⁺, 1, 1, 1)
    # t⁺⁺_ = repeat(t⁺⁺, 1, 1, 1)
    # r⁺⁻_ = repeat(r⁺⁻, 1, 1, 1)
    # t⁻⁻_ = repeat(t⁻⁻, 1, 1, 1)
    
    aux1 = similar(R⁻⁺_)
    aux2 = similar(R⁻⁺_)
    aux3 = similar(R⁻⁺_)

    # Compute M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺
    aux1 = I_static .- R⁺⁻_ ⊠ r⁻⁺_
    batch_solve!(aux2, aux1, T⁺⁺_)

    # Compute t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
    aux1 = r⁻⁺_ ⊠ aux2  # r⁻⁺ * M1
    aux3 = T⁻⁻_ ⊠ aux1  # 
    R⁻⁺_  = R⁻⁺_ + aux3

    # T⁺⁺ = t⁺⁺ * M1
    T⁺⁺_ = t⁺⁺_ ⊠ aux2

    # Repeating for mirror-reflected directions
    # Compute M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻
    aux1 = I_static .- r⁻⁺_ ⊠ R⁺⁻_
    batch_solve!(aux2, aux1, t⁻⁻_)

    # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1
    aux3 = R⁺⁻_ ⊠ aux2
    aux1 = t⁺⁺_ ⊠ aux3
    R⁺⁻_  = r⁺⁻_ + aux1
    T⁻⁻_  = T⁺⁺_ ⊠ aux2

    @. r⁻⁺ = r⁻⁺_[:,:,1]
    @. t⁺⁺ = t⁺⁺_[:,:,1]
    @. r⁺⁻ = r⁺⁻_[:,:,1]
    @. t⁻⁻ = t⁻⁻_[:,:,1]

    @. R⁻⁺ = R⁻⁺_[:,:,1]
    @. T⁺⁺ = T⁺⁺_[:,:,1]
    @. R⁺⁻ = R⁺⁻_[:,:,1]
    @. T⁻⁻ = T⁻⁻_[:,:,1]

end

function rt_interaction!(R⁻⁺::AbstractArray{FT,3}, T⁺⁺::AbstractArray{FT,3}, 
                         R⁺⁻::AbstractArray{FT,3}, T⁻⁻::AbstractArray{FT,3}, 
                         r⁻⁺::AbstractArray{FT,3}, t⁺⁺::AbstractArray{FT,3}, 
                         r⁺⁻::AbstractArray{FT,3}, t⁻⁻::AbstractArray{FT,3},
                         I_static::AbstractArray) where {FT}

    rt_interaction_helper!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻, I_static)
    synchronize()

end