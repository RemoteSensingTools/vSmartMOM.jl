# Shared small physical-optics scene for matrix and equivalent-source Jacobians.
function local_jacobian_fixture(FT, polarized, external; gpu=false, n_aerosols=1, nspec=5, nlayers=5)
    cfg = YAML.load_file("test_parameters/JacobianTestFast.yaml")
    delete!(cfg,"absorption")
    rt = cfg["radiative_transfer"]
    rt["float_type"] = string(FT)
    rt["architecture"] = gpu ? "GPU()" : "CPU()"
    rt["polarization_type"] = polarized ? "Stokes_IQU()" : "Stokes_I()"
    rt["surface"] = ["LambertianSurfaceScalar{$FT}(0.05)"]
    rt["greek_beta_cutoff"] = nothing
    cfg["geometry"]["vaz"] = [0.0,37.0]
    cfg["scattering"]["r_max"] = 3.0
    cfg["scattering"]["aerosols"][1]["μ"] = 0.15
    cfg["scattering"]["aerosols"][1]["σ"] = 1.4
    if n_aerosols == 0
        delete!(cfg,"scattering")
    elseif n_aerosols > 1
        first_aerosol = only(cfg["scattering"]["aerosols"])
        cfg["scattering"]["aerosols"] = [merge(deepcopy(first_aerosol),
            Dict("nᵣ"=>1.3+0.1(i-1),"τ_ref"=>0.04/i,"p₀"=>700.0-100(i-1))) for i in 1:n_aerosols]
    end
    cfg["atmospheric_profile"]["profile_reduction"] = nlayers
    p = read_parameters(cfg)
    p.spec_bands[1] = FT.(range(first(p.spec_bands[1]),last(p.spec_bands[1]);length=nspec))
    model, lin = model_from_parameters(LinMode(),p;external_solar=external)
    @test model.quad_points.external_solar == external
    # Supplied absorption isolates propagation from spectroscopy. Each gas
    # column changes one layer, and overlying beam derivatives are
    # nonzero below it. This is not a line-list/spectroscopy benchmark.
    ns,nz = size(model.τ_abs[1])
    model.τ_abs[1] .= [FT(0.01z*(1+s/ns)) for s in 1:ns,z in 1:nz]
    lin.τ̇_abs[1] .= 0
    for z in 1:nz
        lin.τ̇_abs[1][z,:,z] .= model.τ_abs[1][:,z] ./ FT(0.2)
    end
    return model,lin
end

local_jacobian_forward(model) = model.quad_points.external_solar ?
    (;toa=rt_run_toa(model),boa=nothing) : rt_run(model)
