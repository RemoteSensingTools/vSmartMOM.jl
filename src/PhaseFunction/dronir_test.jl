function compute_mie(size_param::Real, ref_idx::Number, angles::Vector{Float64})
    y = size_param * ref_idx
    x_stop = size_param + 4 * size_param^0.3333 + 2.0

    y_modulus = abs(y)
    n_max = round(Int, max(x_stop, y_modulus) + 15)

    mu_table = [cos(theta) for theta in angles]
    N_angles = length(angles)
    idx_forward_scatter = -1
    idx_backscatter = -1
    for i = 1:N_angles
        if angles[i] ≈ 0.0
            idx_forward_scatter = i
        end
        if angles[i] - π ≈ 0.0
            idx_backscatter = i
        end
    end

    d = zeros(ComplexF64, n_max+1)

    for n = n_max-1:-1:1
        rn = n+1
        d[n] = (rn/y) - (1 / (d[n+1] + rn/y))
    end

    pi0 = zeros(Float64, N_angles)
    pi1 = ones(Float64, N_angles)
    piX = zeros(Float64, N_angles)

    s1 = zeros(ComplexF64, N_angles)
    s2 = zeros(ComplexF64, N_angles)
    tau = zeros(Float64, N_angles)

    psi0 = cos(size_param)
    psi1 = sin(size_param)
    chi0 = -sin(size_param)
    chi1 = cos(size_param)

    S = zeros(Float64, N_angles, 4)

    xi1 = ComplexF64(psi1, -chi1)
    Qsca = 0.0

    for n = 1:round(Int, x_stop)
        fn = (2n+1) / (n*(n+1))
        psi = (2n-1) * psi1/size_param - psi0
        chi = (2n-1) * chi1/size_param - chi0
        xi = ComplexF64(psi, -chi)
        t_a = d[n] / ref_idx + n/size_param
        t_b = d[n] * ref_idx + n/size_param
        an = (t_a * psi - psi1) / (t_a * xi - xi1)
        bn = (t_b * psi - psi1) / (t_b * xi - xi1)
        Qsca += (2n+1) * (abs(an)^2 + abs(bn)^2)

        for j = 1:N_angles
            piX[j] = pi1[j]
            tau[j] = n * mu_table[j] * piX[j] - (n+1)*pi0[j]
            t = (-1)^n
            p = (-1)^(n-1)
            @show n, mu_table[j], piX[j],  tau[j], pi0[j]
            s1[j] = s1[j] + fn * (an*piX[j] + bn*tau[j])
            s2[j] = s2[j] + fn * (an*tau[j] + bn*piX[j])
        end

        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = ComplexF64(psi1, -chi1)

        for i = 1:N_angles
            pi1[i] = (2n+1)/n * mu_table[i] * pi1[i] - (n+1)*pi0[i] / n
            pi0[i] = piX[i]
        end
    end

    Qsca *= 2/size_param^2
    Qext = (idx_forward_scatter > 0) ? (4 / size_param^2) * real(s1[idx_forward_scatter]) : NaN
    Qback = (idx_backscatter > 0) ? (4 / size_param^2) * abs(s1[idx_backscatter])^2 : NaN

    for i = 1:N_angles
        S[i,1] = 0.5 * (abs(s1[i])^2 + abs(s2[i])^2)
        S[i,2] = -0.5 * (abs(s1[i])^2 - abs(s2[i])^2)
        S[i,3] = real(s2[i] * conj(s1[i]))
        S[i,4] = imag(s2[i] * conj(s1[i]))
    end

    return S, Qsca, Qext, Qback

end