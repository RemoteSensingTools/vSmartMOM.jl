# Algebra experiment, not the production factored-core implementation.
# From test/: julia --project=. ../docs/dev_notes/jacobian_batched/core_basis_check.jl
# Compare early state expansion with delayed contraction after actual doubling,
# including the above-layer solar attenuation derivative. No timings are claimed.
using Test, Random, LinearAlgebra, vSmartMOM
using vSmartMOM.Scattering: Stokes_I, Stokes_IQU
using vSmartMOM.InelasticScattering: noRS
const C = vSmartMOM.CoreRT

function check_basis(FT, nstokes, ndoubl)
    rng = MersenneTwister(891)
    n, ns, nb, np = 6nstokes, 7, 8, 60
    pol = nstokes == 1 ? Stokes_I{FT}() : Stokes_IQU{FT}()
    a, da = C.make_added_layer(LinMode(), noRS(), FT, Array, nb, (n,n), ns)
    b, db = C.make_added_layer(LinMode(), noRS(), FT, Array, np, (n,n), ns)
    forward = (:r⁻⁺,:t⁺⁺,:r⁺⁻,:t⁻⁻,:j₀⁺,:j₀⁻)
    tangent = (:ap_ṙ⁻⁺,:ap_ṫ⁺⁺,:ap_ṙ⁺⁻,:ap_ṫ⁻⁻,:ap_J̇₀⁺,:ap_J̇₀⁻)
    coefficients = randn(rng,FT,nb,np,ns) .* FT(0.2)
    function contract(x)
        out = zeros(FT,size(x,1),size(x,2),ns,np)
        for s in 1:ns
            out[:,:,s,:] .= reshape(
                reshape(x[:,:,s,:],:,nb)*coefficients[:,:,s],size(x,1),size(x,2),np)
        end
        out
    end
    for (f,df) in zip(forward,tangent)
        x, dx = getproperty(a,f), getproperty(da,df)
        x .= rand(rng,FT,size(x)) .* FT(0.01)
        dx .= randn(rng,FT,size(dx)) .* FT(0.001)
        getproperty(b,f) .= x
        getproperty(db,df) .= contract(dx)
    end
    # Sources already include the fixed forward attenuation above this layer.
    # Early expansion differentiates it immediately; delayed expansion appends
    # -(dτ_above/μ₀)J after doubling, because J is linear in incident beam flux.
    μ₀ = FT(0.8)
    above = randn(rng,FT,ns,np) .* FT(0.02)
    attenuation_tangent = reshape(-above ./ μ₀,1,1,ns,np)
    db.ap_J̇₀⁺ .+= b.j₀⁺ .* attenuation_tangent
    db.ap_J̇₀⁻ .+= b.j₀⁻ .* attenuation_tangent
    dtau = randn(rng,FT,ns,nb) .* FT(0.001)
    dtau_state = zeros(FT,ns,np)
    for s in 1:ns
        dtau_state[s,:] .= transpose(coefficients[:,:,s])*dtau[s,:]
    end
    ident = repeat(Matrix{FT}(I,n,n),1,1,ns)
    e = fill(FT(0.95),ns)
    C.doubling_allparams!(pol,true,copy(e),ndoubl,a,da,ident,CPU(),dtau,μ₀)
    C.doubling_allparams!(pol,true,copy(e),ndoubl,b,db,ident,CPU(),dtau_state,μ₀)
    tol = FT === Float64 ? FT(1e-11) : FT(3e-5)
    for (f,df) in zip(forward,tangent)
        @test getproperty(a,f) ≈ getproperty(b,f) rtol=tol atol=tol
        delayed = contract(getproperty(da,df))
        if f in (:j₀⁺,:j₀⁻)
            delayed .+= getproperty(a,f) .* attenuation_tangent
        end
        @test delayed ≈ getproperty(db,df) rtol=tol atol=tol*FT(0.01)
    end
end

@testset "Eight core directions → sixty state columns after doubling" begin
    for FT in (Float32,Float64), nstokes in (1,3), ndoubl in (3,6)
        check_basis(FT,nstokes,ndoubl)
    end
end
