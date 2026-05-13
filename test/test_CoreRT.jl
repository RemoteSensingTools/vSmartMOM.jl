using Test
using vSmartMOM
using vSmartMOM.CoreRT

const CORE_RT_BENCHMARK_DIR = joinpath(@__DIR__, "benchmarks")

@testset "compare against 6SV1" begin 

    include(joinpath(CORE_RT_BENCHMARK_DIR, "6SV1_R_trues.jl"))
    R_modeled_all = zeros(6, 3, 3, 16);
    R_deltas_all = zeros(6, 3, 3, 16);
    parameters = parameters_from_yaml(joinpath(CORE_RT_BENCHMARK_DIR, "6SV1_1.yaml"));

    function test_against_6SV1(case_i, azs, szas, λ, τ, ρ)
        for sza_i in 1:length(szas)
            for az_i in 1:length(azs)
                R_true = R_trues[case_i]
                R_modeled = view(R_modeled_all, case_i, :, :, :)
                R_deltas = view(R_deltas_all, case_i, :, :, :)
                parameters.spec_bands = [[1e7/λ, 1e7/λ + 1]]
                parameters.vaz = repeat([azs[az_i]], 16)
                parameters.sza = szas[sza_i]
                parameters.brdf = [vSmartMOM.CoreRT.LambertianSurfaceScalar(ρ)]
                model = model_from_parameters(parameters);
                model.τ_rayl[1] .= τ
                # rt_run returns radiance factor L = I/F₀; 6SV1 reference is reflectance R = πL/μ₀.
                R_modeled[sza_i, az_i, :] = π * vSmartMOM.CoreRT.rt_run(model, i_band=1)[1][:,1,1] / model.quad_points.μ₀
                R_deltas[sza_i, az_i, :] = abs.(R_true[sza_i][az_i] - R_modeled[sza_i, az_i, :]) ./ R_true[sza_i][az_i]
            end
        end
        return maximum(R_deltas_all[case_i,:,:,:])
    end

    ϵ = 0.006

    @test test_against_6SV1(1, [180, 90, 0], [23.0739, 53.1301, 78.4630], 530, 0.1, 0.0) < ϵ
    @test test_against_6SV1(2, [180, 90, 0], [0.0001, 36.8699, 66.4218], 530, 0.1, 0.25) < ϵ
    @test test_against_6SV1(3, [180, 90, 0], [0.0001, 36.8699, 66.4218], 440, 0.25, 0.0) < ϵ
    @test test_against_6SV1(4, [180, 90, 0], [23.0739, 53.1301, 78.4630], 440, 0.25, 0.25) < ϵ
    @test test_against_6SV1(5, [180, 90, 0], [23.0739, 53.1301, 78.4630], 360, 0.50, 0.0) < ϵ
    @test test_against_6SV1(6, [180, 90, 0], [0.0001, 36.8699, 66.4218], 360, 0.50, 0.25) < ϵ

end

@testset "rt_run_streams: Fourier+stream recovery of rt_run" begin
    # Verify that rt_run_streams + a manual Fourier+nearest-stream
    # reconstruction (mirroring postprocessing_vza.jl) reproduces
    # rt_run's R_SFI output exactly. This is the bit-exact plumbing
    # check for the per-moment streams export (Phase H) that lets
    # ExoOptics's disk integrator do one rt_run per band instead of
    # one per disk pixel.
    parameters = parameters_from_yaml(joinpath(CORE_RT_BENCHMARK_DIR, "natraj.yaml"))
    parameters.spec_bands = [[1e7/360.0, 1e7/360.0 + 1]]
    # A handful of VZA × VAZ pairs, all at the same SZA — same Coulson
    # geometry the rest of test_CoreRT.jl uses.
    parameters.vza = [11.4783, 23.0739, 50.2082, 73.7398]
    parameters.vaz = [0.0, 60.0, 120.0, 180.0]
    parameters.sza = acosd(0.2)

    # Direct rt_run (the reference output)
    model_a = model_from_parameters(parameters)
    model_a.τ_rayl[1] .= 0.5
    R_direct, _T_direct = vSmartMOM.CoreRT.rt_run(model_a)

    # Streams output (separate model — both build from the same
    # `parameters` so optical depth + Fourier-loop bounds match)
    model_b = model_from_parameters(parameters)
    model_b.τ_rayl[1] .= 0.5
    streams = vSmartMOM.CoreRT.rt_run_streams(model_b)

    # Reconstruct R from streams: per-(vza, vaz), Fourier-sum the
    # per-moment J⁻ at the nearest quadrature stream to cos(vza),
    # weighted by (Fourier weight) · (cos m·φ) for I,Q and
    # (sin m·φ) for U,V. This is exactly what postprocessing_vza.jl
    # does internally — we just do it offline from the streams.
    pol_n = streams.pol_n
    nVZA  = length(parameters.vza)
    nSpec = size(R_direct, 3)
    R_recon = zeros(eltype(R_direct), nVZA, pol_n, nSpec)
    for (i, vza_deg) in enumerate(parameters.vza)
        vaz_deg = parameters.vaz[i]
        # vSmartMOM's nearest-point logic (matches postprocessing_vza._precompute_vza_weights)
        iμ = argmin(abs.(streams.qp_μ .- cosd(vza_deg)))
        istart = (iμ - 1) * pol_n + 1
        iend   = iμ * pol_n
        for mi in eachindex(streams.J⁻_per_m)
            m = mi - 1
            cosmφ = cosd(m * vaz_deg)
            sinmφ = sind(m * vaz_deg)
            stokes_weights = pol_n == 1 ? cosmφ : [cosmφ, cosmφ, sinmφ, sinmφ][1:pol_n]
            for s in 1:nSpec
                slice = streams.J⁻_per_m[mi][istart:iend, 1, s]
                R_recon[i, :, s] .+= streams.weight[mi] .* stokes_weights .* slice
            end
        end
    end

    @test size(R_recon) == size(R_direct)
    @test isapprox(R_recon, R_direct; atol = 1e-12, rtol = 1e-10)
end

@testset "compare against natraj paper" begin

    I_modeled_all = zeros(7, 16);
    I_deltas_all  = zeros(7, 16);
    Q_modeled_all = zeros(7, 16);
    Q_deltas_all  = zeros(7, 16);
    U_modeled_all = zeros(7, 16);
    U_deltas_all  = zeros(7, 16);

    μ = [0.02, 0.06, 0.10, 0.16, 0.20, 0.28, 0.32, 0.40, 0.52, 0.64, 0.72, 0.84, 0.92, 0.96, 0.98, 1.00]
    ϕs = collect(0.0:30.0:180.0)
    τ = 0.5

    include(joinpath(CORE_RT_BENCHMARK_DIR, "natraj_trues.jl"))
    parameters = parameters_from_yaml(joinpath(CORE_RT_BENCHMARK_DIR, "natraj.yaml"));

    parameters.spec_bands = [[1e7/360.0, 1e7/360.0 + 1]]
    parameters.vza = acosd.(μ)
    parameters.sza = acosd.(0.2)

    for ϕ_i in 1:length(ϕs)
        parameters.vaz = repeat([ϕs[ϕ_i]], 16)
        model = model_from_parameters(parameters);
        model.τ_rayl[1] .= τ

        # rt_run returns radiance factor L = I/F₀; Natraj table is reflectance R = πL.
        R = π * CoreRT.rt_run(model, i_band=1)[1]
        @show size(R)

        I_modeled_all[ϕ_i,:] = R[:,1,1]
        Q_modeled_all[ϕ_i,:] = R[:,2,1]
        U_modeled_all[ϕ_i,:] = R[:,3,1]

        I_deltas_all[ϕ_i,:] = abs.(I_trues[:,ϕ_i] - I_modeled_all[ϕ_i,:]) ./ I_trues[:,ϕ_i]
        Q_deltas_all[ϕ_i,:] = abs.(Q_trues[:,ϕ_i] - Q_modeled_all[ϕ_i,:]) ./ Q_trues[:,ϕ_i]
        U_deltas_all[ϕ_i,:] = abs.(U_trues[:,ϕ_i] - U_modeled_all[ϕ_i,:]) ./ U_trues[:,ϕ_i]
    end

    # Tightened tolerances after switching natraj.yaml to GaussLegQuad +
    # NoTruncation: Float64 max rel-err on F64 measures 0.02% (I), 0.14% (Q,
    # |modeled| ≥ 0.01), and 0.01% (U, |modeled| ≥ 0.01) — see
    # docs/src/pages/benchmarks.md. Limits sit ~2× over measured for the
    # noisy-Q rows and ~5× over the very stable I/U rows.
    @test maximum(I_deltas_all) < 5e-4
    @test maximum(Q_deltas_all[findall(i -> i >= 0.01, Q_modeled_all)]) < 2.5e-3
    @test maximum(filter(!isnan, U_deltas_all[findall(i -> i >= 0.01, U_modeled_all)])) < 5e-4

end
