using Test, YAML, LinearAlgebra, Logging, vSmartMOM
using vSmartMOM.CoreRT

include("local_jacobian_fixture.jl")

function check_local_jacobian(FT,polarized,external; gpu=false, n_aerosols=1,
                              selected_columns=nothing, expected_basis=nothing)
    model,lin = local_jacobian_fixture(FT,polarized,external;gpu,n_aerosols)
    AT = CoreRT.array_type(model)
    ng = size(lin.τ̇_abs[1],1)
    layout = if selected_columns === nothing
        nothing
    else
        names = ["native_$p" for p in selected_columns]
        keys = [ParameterKey(:atmosphere,Symbol(name)) for name in names]
        push!(names,"albedo")
        push!(keys,ParameterKey(:surface,:albedo;component=1,band=1))
        np = length(selected_columns)
        ActiveParameterLayout(keys,names,keys,names,selected_columns,np;
                              surface_columns=np+1:np+1)
    end
    cache = CoreRT.build_m_invariant_cache_lin(1,model,lin;active_layout=layout)
    local_cache = CoreRT.build_local_jacobian_cache(1,model,lin,cache)
    if expected_basis !== nothing
        @test size(local_cache.basis_tau,2) == expected_basis
        @test all(size(l.coefficients,2) == expected_basis for l in local_cache.layers)
        @test CoreRT.local_jacobian_size(n_aerosols,cache.selection) == expected_basis
        @test CoreRT.use_local_jacobian(:auto,layout,1,n_aerosols,cache.selection) ==
              (expected_basis < length(selected_columns))
    end
    rs = vSmartMOM.InelasticScattering.noRS{FT}()
    tol = FT === Float32 ? FT(3e-4) : FT(2e-10)
    for m in (0,1,2)
        physical,dp,_ = CoreRT.constructCoreOpticalProperties(rs,1,m,model,lin,cache)
        factored,df,_ = CoreRT.construct_local_optical_jacobians(rs,1,m,model,lin,local_cache)
        for z in eachindex(physical)
            fp,jp = CoreRT.expandOpticalProperties(physical[z],dp[z],AT)
            fl,jl = factored[z],df[z]
            # Choosing derivative coordinates must not change the forward
            # mixture, even by an ulp before repeated Float32 doubling.
            @test Array(fp.τ) == Array(fl.τ)
            @test Array(fp.ϖ) == Array(fl.ϖ)
            @test eltype(fp.ϖ) === eltype(fl.ϖ) === FT
            for name in (:Z⁺⁺,:Z⁻⁺,:Z₀⁺,:Z₀⁻)
                a,b = getproperty(fp,name),getproperty(fl,name)
                if a === nothing
                    @test b === nothing
                else
                    @test eltype(a) === eltype(b) === FT
                    @test Array(a) == Array(b)
                end
            end
            @test Array(jp.τ̇) ≈ Array(jl.τ̇) rtol=tol atol=eps(FT)
            @test Array(jp.ϖ̇) ≈ Array(jl.ϖ̇) rtol=tol atol=eps(FT)
            for name in (:Ż⁺⁺,:Ż⁻⁺,:Ż₀⁺,:Ż₀⁻)
                target = getproperty(jp,name)
                target === nothing && continue
                basis = getproperty(jl.basis,name)
                ns = size(jl.coefficients,1)
                size(basis,3) == 1 && ns > 1 && (basis=repeat(basis,1,1,ns,1))
                out = similar(basis,size(basis,1),size(basis,2),ns,size(jp.τ̇,2))
                CoreRT.contract_local_jacobian!(out,basis,jl.coefficients)
                a, b = Array(out), Array(target)
                # Use elementwise tolerances, not an array norm whose absolute
                # floor grows with the number of phase entries. Absorption
                # must add exactly zero phase columns in both representations.
                # Each Float32 scatterer addition rounds another quotient.
                # Allow 10 eps per addition near cancelled phase entries;
                # retain the relative bound on resolved derivatives. The
                # three-aerosol case has ~2e-6 differences at ~1e-3 entries.
                phase_atol = (FT === Float64 ? 32 : 10max(1,n_aerosols)) * eps(FT)
                @test all(isapprox.(a,b;rtol=tol,atol=phase_atol))
                gas_start = selected_columns === nothing ? 2+7n_aerosols :
                    1+count(p->p<=1+7n_aerosols,selected_columns)
                @test all(iszero, a[:,:,:,gas_start:end])
                @test all(iszero, b[:,:,:,gas_start:end])
            end
        end
    end
    run_basis(mode;kwargs...) = if layout === nothing
        rt_run(model,lin,n_aerosols,ng,1;jacobian_basis=mode,kwargs...)
    else
        plan = JacobianPlan(OCO_RRS_synth(),layout.keys,layout.parameter_names,[layout])
        rt_run(model,PlannedRTModelLin(lin,plan);jacobian_basis=mode,kwargs...)
    end
    early = run_basis(:physical)
    late = run_basis(:local)
    if layout !== nothing
        source = run_basis(:local;jacobian_adding=:source)
        auto = run_basis(:auto)
        for result in (source,auto)
            @test result.toa ≈ early.toa rtol=tol atol=10eps(FT)
            @test result.toa_jacobian ≈ early.toa_jacobian rtol=tol atol=10eps(FT)
            if !external
                @test result.boa ≈ early.boa rtol=tol atol=10eps(FT)
                @test result.boa_jacobian ≈ early.boa_jacobian rtol=tol atol=10eps(FT)
            end
        end
    end
    for (x,y) in zip(early,late)
        if x === nothing || y === nothing
            @test x === y
        else
            @test x ≈ y rtol=tol atol=10eps(FT)
        end
    end
    gas_native = 1+7n_aerosols+2
    if FT === Float64 && (selected_columns === nothing || gas_native in selected_columns)
        # Independent gas finite difference through the forward solver,
        # including the beam attenuation above every subsequent layer.
        z=2; h=1e-5
        original=copy(model.τ_abs[1][:,z])
        direction=copy(lin.τ̇_abs[1][z,:,z])
        model.τ_abs[1][:,z] .= original .+ h.*direction
        plus=local_jacobian_forward(model)
        model.τ_abs[1][:,z] .= original .- h.*direction
        minus=local_jacobian_forward(model)
        model.τ_abs[1][:,z] .= original
        p=selected_columns === nothing ? gas_native : findfirst(==(gas_native),selected_columns)
        @test late.toa_jacobian[:,:,:,p] ≈ (plus.toa-minus.toa)/(2h) rtol=2e-6 atol=1e-9
        if !external
            @test late.boa_jacobian[:,:,:,p] ≈ (plus.boa-minus.boa)/(2h) rtol=2e-6 atol=1e-9
        end
    end
end

# A zero aerosol column has Rayleigh-only forward phase moments, including
# exact zeros above m=2. Its AOT derivative can still introduce those moments.
# Use a one-sided, second-order difference to stay at nonnegative aerosol depth.
function check_zero_aerosol_jacobian(external)
    model,lin = local_jacobian_fixture(Float64,true,external)
    ng = size(lin.τ̇_abs[1],1)
    direction = copy(lin.τ̇_aer[1][1,1,:,:])
    model.τ_aer[1] .= 0
    lin.τ̇_aer[1][:,2:end,:,:] .= 0
    lin.τ̇_aer_psurf[1] .= 0
    base = local_jacobian_forward(model)
    h = 1e-6
    model.τ_aer[1][1,:,:] .= h .* direction
    one_step = local_jacobian_forward(model)
    model.τ_aer[1][1,:,:] .= 2h .* direction
    two_steps = local_jacobian_forward(model)
    model.τ_aer[1] .= 0
    for mode in (:physical,:local)
        result = rt_run(model,lin,1,ng,1;jacobian_basis=mode)
        @test result.toa ≈ base.toa rtol=1e-10 atol=1e-12
        if external
            @test result.boa === nothing
        else
            @test result.boa ≈ base.boa rtol=1e-10 atol=1e-12
        end
        for (field,jac) in ((:toa,:toa_jacobian),(:boa,:boa_jacobian))
            getproperty(base,field) === nothing && continue
            fd = (-3getproperty(base,field) .+ 4getproperty(one_step,field) .-
                  getproperty(two_steps,field)) ./ (2h)
            @test getproperty(result,jac)[:,:,:,2] ≈ fd rtol=3e-5 atol=2e-8
        end
    end
end

@testset "Local optical basis CPU" begin
    with_logger(NullLogger()) do
        for FT in (Float64,Float32), polarized in (false,true), external in (false,true)
            check_local_jacobian(FT,polarized,external)
        end
        for na in (0,2)
            check_local_jacobian(Float64,true,true;n_aerosols=na)
        end
        check_local_jacobian(Float32,true,true;n_aerosols=3)
        # Three fixed-microphysics species need only τ, ϖ and three mixture
        # directions. A mixed plan retains noncontiguous microphysics from
        # different species, exercising both phase and coefficient indexing.
        fixed = [1,2,7,9,14,16,21,24]
        mixed = [1,2,7,9,10,13,14,16,19,21,24]
        for external in (false,true)
            check_local_jacobian(Float64,true,external;n_aerosols=3,
                selected_columns=fixed,expected_basis=5)
        end
        check_local_jacobian(Float32,true,true;n_aerosols=3,
            selected_columns=fixed,expected_basis=5)
        check_local_jacobian(Float64,true,true;n_aerosols=3,
            selected_columns=mixed,expected_basis=8)
        for external in (false,true)
            check_zero_aerosol_jacobian(external)
        end
    end
end
if get(ENV,"VSMARTMOM_JACOBIAN_GPU_TEST","false") == "true"
    using CUDA
    CUDA.allowscalar(false)
    @testset "Local optical basis CUDA" begin
        with_logger(NullLogger()) do
            check_local_jacobian(Float64,true,true;gpu=true)
            check_local_jacobian(Float32,true,false;gpu=true)
            check_local_jacobian(Float32,true,true;gpu=true,n_aerosols=3)
            check_local_jacobian(Float32,true,true;gpu=true,n_aerosols=3,
                selected_columns=[1,2,7,9,14,16,21,24],expected_basis=5)
            check_local_jacobian(Float64,true,true;gpu=true,n_aerosols=3,
                selected_columns=[1,2,7,9,10,13,14,16,19,21,24],expected_basis=8)
        end
    end
end
