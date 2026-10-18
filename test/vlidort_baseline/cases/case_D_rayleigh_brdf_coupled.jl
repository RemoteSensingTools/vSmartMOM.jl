using Test
using vSmartMOM
using vSmartMOM.CoreRT

# Case D — one polarized Rayleigh layer over a BRDF surface, against VLIDORT
# 2.8.3 (reference_data/rayleigh_brdf_coupled_truth.jl). Unlike vacuum checks,
# multiple scattering between atmosphere and surface exercises every element
# of the surface Mueller matrix and its Fourier moments, off the principal
# plane. Surfaces: GISS Cox–Munk (with VLIDORT's exact direct bounce, matched
# by vSmartMOM's glint correction), Lambertian, and Ross–Li with vSmartMOM's
# 80° kernel zenith cap.
include(joinpath(REF_DIR, "rayleigh_brdf_coupled_truth.jl"))

function coupled_surface(s)
    s.kind == :coxmunk && return CoxMunkSurface(wind_speed = s.wind_speed,
        n_water = complex(s.n_water, 0.0), include_whitecaps = false,
        shadowing = s.shadowing)
    s.kind == :lambertian && return LambertianSurfaceScalar(s.albedo)
    s.kind == :rossli && return CoreRT.RossLiSurfaceScalar(s.fvol, s.fgeo, s.fiso)
    error("unknown surface $(s.kind)")
end

"TOA up and BOA down Stokes IQU on the reference grid [sza, vza, raz, stokes]."
function coupled_run(case)
    d = Dict{String,Any}(
        "radiative_transfer" => Dict{String,Any}(
            "spec_bands" => ["[12987.0]"], "surface" => ["LambertianSurfaceScalar(0.1)"],
            "nstreams" => case.nstreams, "quadrature_type" => "GaussLegQuad()",
            "polarization_type" => "Stokes_IQU()", "truncation" => "NoTruncation()",
            "depol" => case.depol, "float_type" => "Float64", "architecture" => "CPU()"),
        # All (VZA, RAZ) pairs as views of one model per SZA.
        "geometry" => Dict{String,Any}("sza" => 30.0,
            "vza" => repeat(COUPLED_VZA_DEG, outer = length(COUPLED_RAZ_DEG)),
            "vaz" => repeat(COUPLED_RAZ_DEG, inner = length(COUPLED_VZA_DEG)), "obs_alt" => [0]),
        "atmospheric_profile" => Dict{String,Any}("T" => [250.0], "p" => [500.0, 1000.0],
            "profile_reduction" => -1))
    params = read_parameters(d)
    # The surface must be set before construction: it sets the Fourier cut-off.
    params.brdf = [coupled_surface(case.surface)]
    toa = zeros(length(COUPLED_SZA_DEG), length(COUPLED_VZA_DEG), length(COUPLED_RAZ_DEG), 3)
    boa = similar(toa)
    nv = length(COUPLED_VZA_DEG)
    for (i, sza) in enumerate(COUPLED_SZA_DEG)
        params.sza = sza
        model = model_from_parameters(params)
        fill!(model.τ_rayl[1], case.tau)                    # one layer
        fill!(model.τ_abs[1], case.tau * (1 - case.omega))  # VLIDORT's ω = 1 − 1e-9
        R, T = rt_run(model)
        for k in eachindex(COUPLED_RAZ_DEG), j in 1:nv
            toa[i, j, k, :] .= R[(k - 1) * nv + j, :, 1]
            boa[i, j, k, :] .= T[(k - 1) * nv + j, :, 1]
        end
    end
    return toa, boa
end

@testset "VLIDORT baseline: Case D — Rayleigh layer over BRDF surfaces" begin
    flip = reshape([1.0, -1.0, -1.0], 1, 1, 1, 3)      # VLIDORT → vSmartMOM Q, U signs
    for case in COUPLED_CASES
        @testset "$(case.name)" begin
            toa, boa = coupled_run(case)
            for (label, model, truth) in (("TOA", toa, case.toa_up .* flip),
                                          ("BOA", boa, case.boa_dn .* flip))
                scale = maximum(abs, truth[:, :, :, 1])
                # VLIDORT's BOA output is diffuse only; vSmartMOM's downwelling at the
                # μ = μ₀ node also carries the direct solar beam, so that view is excluded.
                keep = [!(label == "BOA" && COUPLED_VZA_DEG[j] == COUPLED_SZA_DEG[i])
                        for i in eachindex(COUPLED_SZA_DEG), j in eachindex(COUPLED_VZA_DEG),
                            _ in eachindex(COUPLED_RAZ_DEG)]
                rel_I = maximum(abs.(model[:, :, :, 1] .- truth[:, :, :, 1])[keep] ./
                                truth[:, :, :, 1][keep])
                err_QU = maximum(abs.(model[:, :, :, 2:3] .- truth[:, :, :, 2:3])) / scale
                @info "Case D $(case.name) $label" rel_I err_QU
                # Both codes use 16 Gauss streams; residuals are quadrature-level
                # (Lambertian control: ~3e-5 in I).
                @test rel_I < 2e-4
                @test err_QU < 5e-5
            end
        end
    end
end
