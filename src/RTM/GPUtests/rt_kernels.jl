@kernel function rt_interaction!(R⁻⁺_, T⁺⁺_, R⁺⁻_, T⁻⁻_, r⁻⁺_, t⁺⁺_, r⁺⁻_, t⁻⁻_, aux1_, aux2_, aux3_, pivot_, n)
    N = @index(Global)
    
    # Views:
    R⁻⁺   = view(R⁻⁺_, :, :, N)
    T⁺⁺   = view(T⁺⁺_, :, :, N)
    R⁺⁻   = view(R⁺⁻_, :, :, N)
    T⁻⁻   = view(T⁻⁻_, :, :, N) 
    r⁻⁺   = view(r⁻⁺_, :, :, N)
    t⁺⁺   = view(t⁺⁺_, :, :, N)
    r⁺⁻   = view(r⁺⁻_, :, :, N)
    t⁻⁻   = view(t⁻⁻_, :, :, N)
    aux1  = view(aux1_, :, :, N)
    aux2  = view(aux2_, :, :, N)
    aux3  = view(aux3_, :, :, N)
    pivot = view(pivot_, :,  N)

    # Compute M1 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺
    mat_multiply!(aux1, R⁺⁻, r⁻⁺,n)                 # aux1  =  R⁺⁻ * r⁻⁺
    eye_minus_matrix!(aux1,n)                       # aux1  = (I - R⁺⁻ * r⁻⁺)
    LU_decomposition!(aux1, pivot, n)
    LU_solve!(aux1, T⁺⁺, aux2, aux3, pivot, n)     # aux2 = (I - R⁺⁻ * r⁻⁺) \ T⁺⁺
          
    # Compute t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1
    mat_multiply!(aux1, r⁻⁺, aux2,n)   # r⁻⁺ * M1
    mat_multiply!(aux3, T⁻⁻, aux1,n)   # T⁻⁻ * r⁻⁺ * M1
    mat_add!(R⁻⁺,aux3,n)               # t_R⁻⁺ = R⁻⁺ + T⁻⁻ * r⁻⁺ * M1

    # t_T⁺⁺ = t⁺⁺ * M1
    mat_multiply!(T⁺⁺, t⁺⁺, aux2,n)

    # Repeating for mirror-reflected directions
    # Compute M1 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻
    mat_multiply!(aux1, r⁻⁺, R⁺⁻,n)              # r⁻⁺ * R⁺⁻
    eye_minus_matrix!(aux1,n)                    # (I - r⁻⁺ * R⁺⁻)
    LU_decomposition!(aux1, pivot, n)
    LU_solve!(aux1, t⁻⁻, aux2, aux3, pivot, n)  # aux2 = (I - r⁻⁺ * R⁺⁻) \ t⁻⁻

    # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1
    mat_multiply!(aux3, R⁺⁻, aux2,n)   # R⁺⁻ * M1
    mat_multiply!(aux1, t⁺⁺, aux3,n)   # t⁺⁺ * R⁺⁻ * M1
    mat_add!(R⁺⁻, r⁺⁻, aux1,n)     # t_R⁺⁻ = r⁺⁻ + t⁺⁺ * R⁺⁻ * M1

    # t_T⁻⁻ = T⁺⁺ * M1
    mat_multiply!(T⁻⁻, T⁺⁺, aux2,n)
end


@kernel function rt_interaction_rup!(R⁻⁺_, T⁺⁺_, R⁺⁻_, T⁻⁻_, r⁻⁺_, t⁺⁺_, r⁺⁻_, t⁻⁻_)
    N = @index(Global)
    # Views:
    R⁻⁺   = view(R⁻⁺_, :, :, N)
    T⁺⁺   = view(T⁺⁺_, :, :, N)
    R⁺⁻   = view(R⁺⁻_, :, :, N)
    T⁻⁻   = view(T⁻⁻_, :, :, N) 
    r⁻⁺   = view(r⁻⁺_, :, :, N)
    t⁺⁺   = view(t⁺⁺_, :, :, N)
    r⁺⁻   = view(r⁺⁻_, :, :, N)
    t⁻⁻   = view(t⁻⁻_, :, :, N)
    rt_interaction!(R⁻⁺, T⁺⁺, R⁺⁻, T⁻⁻, r⁻⁺, t⁺⁺, r⁺⁻, t⁻⁻)
end