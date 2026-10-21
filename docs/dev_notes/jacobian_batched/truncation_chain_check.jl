# From test/: julia --project=. ../docs/dev_notes/jacobian_batched/truncation_chain_check.jl
# Independent finite differences from raw Greek coefficients through the
# production δ-BGE fit, normalized Greek tangents, phase matrices and scalar
# aerosol invariant cache. This isolates the truncation chain from Mie/AOD.
using Test, JLD2, LinearAlgebra, vSmartMOM
using vSmartMOM.Scattering
const C = vSmartMOM.CoreRT
const S = vSmartMOM.Scattering
@load "test_pcw/PCW_AerosolOptics_v2.jld" aerosol_optics_PCW
const raw = AerosolOptics(greek_coefs=aerosol_optics_PCW.greek_coefs,
    ω̃=aerosol_optics_PCW.ω̃,k=aerosol_optics_PCW.k,fᵗ=aerosol_optics_PCW.fᵗ)
const fields = (:α,:β,:γ,:δ,:ϵ,:ζ)
const dfields = (:α̇,:β̇,:γ̇,:δ̇,:ϵ̇,:ζ̇)
const direction = map(fields) do f
    g = getproperty(raw.greek_coefs,f)
    d = 0.001 .* g .* cos.(0.4 .* eachindex(g))
    f == :β && (d[1]=0) # Preserve phase normalization at the raw boundary.
    d
end
const lg = linGreekCoefs(map(d->vcat(reshape(d,1,:),zeros(3,length(d))),direction)...)
const lr = linAerosolOptics(lin_greek_coefs=lg,ω̃̇=[0.02,0,0,0],
                          k̇=zeros(4),ḟᵗ=zeros(4))
const mod = δBGE(10)
const value, tangent = S.truncate_phase(mod,raw,lr)
function perturbed(δ)
    g = GreekCoefs(map((f,d)->getproperty(raw.greek_coefs,f) .+ δ.*d,fields,direction)...)
    S.truncate_phase(mod,AerosolOptics(greek_coefs=g,ω̃=raw.ω̃+δ*0.02,k=raw.k,fᵗ=raw.fᵗ))
end
const h = 1e-2 # Above the ill-conditioned normal-equation roundoff floor.
const plus, minus = perturbed(h), perturbed(-h)

@testset "Truncated phase and scalar f chain" begin
    fd_f = (plus.fᵗ-minus.fᵗ)/(2h)
    @test abs(tangent.ḟᵗ[1]) > 1e-8
    @test tangent.ḟᵗ[1] ≈ fd_f rtol=2e-3 atol=2e-7
    for (f,df) in zip(fields,dfields)
        fd = (getproperty(plus.greek_coefs,f)-getproperty(minus.greek_coefs,f))/(2h)
        analytic = getproperty(tangent.lin_greek_coefs,df)[1,:]
        println("TRUNCATION_FD family=$f absolute_error=$(maximum(abs,analytic-fd))")
        @test analytic ≈ fd rtol=2e-3 atol=2e-7
    end
    μ=[0.2,0.6,0.9]
    for pol in (Stokes_I(),Stokes_IQU()), m in (0,1,2)
        z = S.compute_Z_moments(pol,μ,value.greek_coefs,tangent.lin_greek_coefs,m)
        zp = S.compute_Z_moments(pol,μ,plus.greek_coefs,m)
        zm = S.compute_Z_moments(pol,μ,minus.greek_coefs,m)
        for k in 1:2
            @test z[k+2][1,:,:] ≈ (zp[k]-zm[k])/(2h) rtol=2e-3 atol=2e-7
        end
    end
    # Exercise the actual scalar cache with the same nonzero f tangent and
    # simultaneous τ/ω changes, rather than reimplementing its derivative.
    τ=[0.1,0.2,0.3]; dτ=zeros(7,3); dτ[2,:].=[0.03,0.02,0.01]
    cache=C._createAero_invariant(τ,value,dτ,tangent,Array)
    function forward_cache(a,δ)
        t=τ .+ δ.*dτ[2,:]
        (t.*(1-a.fᵗ*a.ω̃),fill((1-a.fᵗ)*a.ω̃/(1-a.fᵗ*a.ω̃),3))
    end
    cp,cm=forward_cache(plus,h),forward_cache(minus,-h)
    @test cache.τ̇[:,2] ≈ (cp[1]-cm[1])/(2h) rtol=2e-3 atol=2e-7
    @test cache.ϖ̇[:,2] ≈ (cp[2]-cm[2])/(2h) rtol=2e-3 atol=2e-7
end
