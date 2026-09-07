using Test, YAML, Logging, vSmartMOM
using vSmartMOM.CoreRT
using Distributions: LogNormal

# A strongly truncated radius grid makes the unnormalized weight sum differ
# appreciably from one. Check the quotient rule independently at fixed nodes.
@testset "Normalized size-distribution tangents" begin
    r = collect(range(0.01,0.2;length=101))
    wr = fill((last(r)-first(r))/(length(r)-1),length(r))
    wr[[1,end]] ./= 2
    μ, σ = log(0.15), log(1.6)
    w, dw = vSmartMOM.Scattering.compute_wₓ(LinMode(),LogNormal(μ,σ),wr,r,last(r))
    @test sum(w) ≈ 1
    @test all(abs.(sum(dw;dims=2)) .< 1e-12)
    for coordinate in 1:2
        h = 1e-6
        dp = coordinate == 1 ? LogNormal(μ+h,σ) : LogNormal(μ,σ+h)
        dm = coordinate == 1 ? LogNormal(μ-h,σ) : LogNormal(μ,σ-h)
        wp = vSmartMOM.Scattering.compute_wₓ(dp,wr,r,last(r))
        wm = vSmartMOM.Scattering.compute_wₓ(dm,wr,r,last(r))
        @test dw[coordinate,:] ≈ (wp-wm)/(2h) rtol=2e-8 atol=1e-10
    end
end

# Independently rebuild the nonlinear model: LinMode's own forward fields
# cannot detect a normalization mismatch between the two constructors.
@testset "Common aerosol reference extinction" begin
    with_logger(NullLogger()) do
        for explicit_reference in (false, true)
            cfg = YAML.load_file("test_parameters/JacobianTestFast.yaml")
            delete!(cfg, "absorption")
            cfg["radiative_transfer"]["polarization_type"] = "Stokes_IQU()"
            cfg["radiative_transfer"]["greek_beta_cutoff"] = nothing
            cfg["geometry"]["vaz"] = [0.0, 37.0]
            scat = cfg["scattering"]
            scat["r_max"] = 3.0
            # Converge radius quadrature before comparing size derivatives:
            # finite differences also move the distribution-dependent nodes.
            scat["nquad_radius"] = 300
            aerosol = only(scat["aerosols"])
            aerosol["μ"] = 0.15
            aerosol["σ"] = 1.4
            aerosol["nᵢ"] = 0.01
            scat["aerosols"] = [aerosol, merge(deepcopy(aerosol),
                Dict("nᵣ"=>1.5, "nᵢ"=>0.02, "μ"=>0.2, "τ_ref"=>0.06))]
            if explicit_reference
                scat["n_ref"] = "1.42 - 0.015im"
                scat["λ_ref"] = 0.7695 # strictly interior: no artificial unity anchor
            end
            params = read_parameters(cfg)
            # The denominator uses each mode's size distribution. Verify that
            # the quotient-rule size derivatives survive the fixed-index change.
            model_size, tangent_size = model_from_parameters(LinMode(),params)
            for ia in 1:2, (coordinate,slot) in ((1,4),(2,5))
                plus, minus = deepcopy(params), deepcopy(params)
                distribution = params.scattering_params.rt_aerosols[ia].aerosol.size_distribution
                # Native Mie columns use LogNormal's log-radius location
                # and log-radius width, not the YAML median/geometric width.
                values = [distribution.μ,distribution.σ]
                h = 1e-5
                vp, vm = copy(values), copy(values)
                vp[coordinate] += h
                vm[coordinate] -= h
                plus.scattering_params.rt_aerosols[ia].aerosol.size_distribution = LogNormal(vp[1],vp[2])
                minus.scattering_params.rt_aerosols[ia].aerosol.size_distribution = LogNormal(vm[1],vm[2])
                fd = (model_from_parameters(plus).τ_aer[1][ia,:,:] -
                      model_from_parameters(minus).τ_aer[1][ia,:,:])/(2h)
                @test tangent_size.τ̇_aer[1][ia,slot,:,:] ≈ fd rtol=3e-4 atol=2e-8
            end
            for external in (false, true)
                model, tangent = model_from_parameters(LinMode(), params; external_solar=external)
                forward = model_from_parameters(params; external_solar=external)
                @test model.τ_aer[1] ≈ forward.τ_aer[1] rtol=2e-11
                result = rt_run(model, tangent, 2, size(tangent.τ̇_abs[1],1), 1)
                run_forward(p) = external ?
                    (;toa=rt_run_toa(model_from_parameters(p;external_solar=true)), boa=nothing) :
                    rt_run(model_from_parameters(p))
                expected = external ? (;toa=rt_run_toa(forward),boa=nothing) : rt_run(forward)
                @test result.toa ≈ expected.toa rtol=2e-10 atol=1e-12
                external || @test result.boa ≈ expected.boa rtol=2e-10 atol=1e-12
                # Index perturbations keep parsed n_ref fixed, including mode 1.
                for ia in 1:2, (field,slot) in ((:nᵣ,2),(:nᵢ,3))
                    h = 1e-5
                    plus, minus = deepcopy(params), deepcopy(params)
                    aplus = plus.scattering_params.rt_aerosols[ia].aerosol
                    aminus = minus.scattering_params.rt_aerosols[ia].aerosol
                    setproperty!(aplus,field,getproperty(aplus,field)+h)
                    setproperty!(aminus,field,getproperty(aminus,field)-h)
                    rp, rm = run_forward(plus), run_forward(minus)
                    column = 1 + 7(ia-1) + slot
                    @test result.toa_jacobian[:,:,:,column] ≈ (rp.toa-rm.toa)/(2h) rtol=3e-4 atol=2e-8
                    external || @test result.boa_jacobian[:,:,:,column] ≈ (rp.boa-rm.boa)/(2h) rtol=3e-4 atol=2e-8
                end
            end
        end
    end
end
