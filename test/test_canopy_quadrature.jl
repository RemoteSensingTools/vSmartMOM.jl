using Test, vSmartMOM, CanopyOptics, LinearAlgebra
using vSmartMOM.CoreRT, vSmartMOM.InelasticScattering

@testset "Canopy extinction uses the scattering leaf measure" begin
    for FT in (Float32,Float64), npol in (1,3), spectral in (false,true), clumped in (false,true)
        pol = npol == 1 ? CoreRT.Stokes_I{FT}() : CoreRT.Stokes_IQU{FT}()
        obs = CoreRT.ObsGeometry(FT(40),FT[0],FT[0],FT[0])
        qp = CoreRT.rt_set_streams(CoreRT.GaussLegQuad(),15,obs,pol,Array)
        n = length(qp.qp_μ)*npol
        added = CoreRT.make_added_layer(noRS{FT}(),FT,Array,(n,n),spectral ? 2 : 1)
        # Deliberately differs from both the former G default and CanopyQuadrature().
        quadrature = CanopyQuadrature(n_leaf=8)
        lad = beta_leaves(FT(0.4),FT(0.7))
        omega_clump = clumped ? FT(0.7) : one(FT)
        canopy = CoreRT.CanopySurface(; soil=CoreRT.LambertianSurfaceScalar(FT(0.1)),
            LAI=FT(3), LAD=lad, canopy_quadrature=quadrature,
            canopy_clumping=ConstantClumping(Ω=omega_clump),
            leaf_reflectance=spectral ? FT[0.4,0.5] : FT(0.4),
            leaf_transmittance=spectral ? FT[0.2,0.3] : FT(0.2),
            leaf_optics_grid=spectral ? FT[500,600] : nothing, grid_unit=:nm)
        spectral_grid = spectral ? FT(1e7)./FT[500,600] : nothing
        CoreRT._init_canopy_cache!(canopy,added,pol,qp,vSmartMOM.Architectures.CPU();
            spec_bands_wn=spectral_grid,m_max=0)

        # Independent normalized leaf-area integral of |normal dot ray|.
        # The one-sided moments analytically integrate leaf azimuth; G uses A.
        theta, weights = CanopyOptics.gauleg(quadrature.n_leaf,zero(FT),FT(pi)/FT(2))
        weights .*= FT.(CanopyOptics.pdf.(lad.LD,FT(2).*theta./FT(pi)))
        weights ./= sum(weights)
        expected = map(collect(qp.qp_μ)) do mu
            sum(zip(theta,weights)) do (angle,weight)
                P,N = CanopyOptics._one_sided_projection_moments(mu,cos(angle),0)
                weight*(P[1]+N[1])*omega_clump
            end
        end
        actual = Array(canopy._cache.G)
        @test eltype(actual) === FT
        @test isapprox(actual,repeat(expected;inner=npol);rtol=FT(100)*eps(FT))
        @test isapprox(actual[qp.iμ₀Nstart],expected[qp.iμ₀];rtol=FT(100)*eps(FT))
        @test all(isfinite,canopy._cache.Z_coarse_pp)
    end
end
