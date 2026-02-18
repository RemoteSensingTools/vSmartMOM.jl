"""
    test_jacobians.jl

Comprehensive perturbation test suite for validating analytic Jacobians 
against finite-difference approximations.

Tests cover:
  1. Forward R consistency (linearized path vs forward-only path)
  2. Surface albedo Jacobian
  3. Aerosol τ_ref Jacobian
  4. Aerosol Mie parameter Jacobians (nᵣ, nᵢ, rₚ, σₚ)
  5. Aerosol vertical profile Jacobians (p₀, σp)
  6. Gas VMR Jacobian (requires variable_molecules)
"""

using vSmartMOM, vSmartMOM.CoreRT
using Distributions: Normal, LogNormal, mean, std
using Test
using TimerOutputs
using Printf

# ──────────────────────────────────────────────────────────────────
# Utilities
# ──────────────────────────────────────────────────────────────────

"""
    relative_error(analytic, fd; floor=1e-12)

Compute element-wise relative error with a floor to avoid 0/0.
Returns the maximum relative error across all elements.
"""
function relative_error(analytic, fd; floor=1e-12)
    denom = max.(abs.(fd), floor)
    return maximum(abs.(analytic .- fd) ./ denom)
end

"""
    mean_relative_error(analytic, fd; floor=1e-12)

Compute mean relative error.
"""
function mean_relative_error(analytic, fd; floor=1e-12)
    denom = max.(abs.(fd), floor)
    return mean(abs.(analytic .- fd) ./ denom)
end

"""
    run_forward(params)

Build forward model from parameters and run forward RT for band 1.
Returns reflectance R with shape (n_vza, n_stokes, n_spec).
"""
function run_forward(params)
    model = model_from_parameters(params)
    R, T, ieJ₀⁺, ieJ₀⁻ = CoreRT.rt_run(model, i_band=1)
    return R
end

"""
    run_linearized(params; NAer, NGas, NSurf)

Build linearized model and run linearized RT for band 1.
Returns (R, T, Ṙ, Ṫ).
"""
function run_linearized(params; NAer=nothing, NGas=nothing, NSurf=nothing)
    lin_mode = LinMode()
    model, lin_model = model_from_parameters(lin_mode, params)
    
    # Auto-detect NAer if not provided
    if isnothing(NAer)
        NAer = isnothing(params.scattering_params) ? 0 : length(params.scattering_params.rt_aerosols)
    end
    # Auto-detect NGas: first dimension of τ̇_abs includes H2O slot
    if isnothing(NGas)
        NGas = size(lin_model.τ̇_abs[1], 1)
    end
    # Default NSurf = 1 (one surface albedo per band)
    if isnothing(NSurf)
        NSurf = 1
    end
    
    R, T, Ṙ, Ṫ = CoreRT.rt_run(model, lin_model, NAer, NGas, NSurf, i_band=1)
    return R, T, Ṙ, Ṫ, model, lin_model, NAer, NGas, NSurf
end

"""
    finite_difference_jacobian(params_base, R_base, perturb_func!, ε)

Compute central finite-difference Jacobian for a single parameter.
`perturb_func!(params, ε)` should modify the parameter in-place.
Returns (R_plus - R_minus) / (2ε).
"""
function finite_difference_jacobian(params_base, perturb_func!, ε)
    # Forward perturbation
    params_plus = deepcopy(params_base)
    perturb_func!(params_plus, +ε)
    R_plus = run_forward(params_plus)
    
    # Backward perturbation
    params_minus = deepcopy(params_base)
    perturb_func!(params_minus, -ε)
    R_minus = run_forward(params_minus)
    
    return (R_plus .- R_minus) ./ (2ε)
end

# ──────────────────────────────────────────────────────────────────
# Test Suite
# ──────────────────────────────────────────────────────────────────

@testset "Linearized RT Jacobian Tests" begin

    # ── Phase 2: End-to-end run & forward consistency ─────────────
    @testset "Phase 2: End-to-end forward consistency" begin
        println("\n" * "="^60)
        println("Phase 2: Loading JacobianTest.yaml and running linearized RT")
        println("="^60)
        
        params = parameters_from_yaml(joinpath(@__DIR__, "test_parameters", "JacobianTest.yaml"))
        
        # Run forward-only
        println("  Running forward-only RT...")
        R_fwd = run_forward(params)
        println("  Forward R shape: ", size(R_fwd))
        println("  Forward R[1,1,1] = ", R_fwd[1,1,1])
        
        # Run linearized
        println("  Running linearized RT...")
        R_lin, T_lin, Ṙ_lin, Ṫ_lin, model, lin_model, NAer, NGas, NSurf = run_linearized(params)
        println("  Linearized R shape: ", size(R_lin))
        println("  Linearized Ṙ shape: ", size(Ṙ_lin))
        println("  NAer=$NAer, NGas=$NGas, NSurf=$NSurf, Nparams=$(NAer*7+NGas+NSurf)")
        
        # Check shapes
        @test size(R_lin) == size(R_fwd)
        
        # Check forward reflectance matches
        max_diff = maximum(abs.(R_lin .- R_fwd))
        println("  Max |R_lin - R_fwd| = ", max_diff)
        @test max_diff < 1e-10  # Should be numerically identical
        
        println("  ✓ Phase 2 passed: forward R matches between linearized and forward-only paths")
    end

    # ── Phase 3a: Surface albedo Jacobian ─────────────────────────
    @testset "Phase 3a: Surface albedo Jacobian" begin
        println("\n" * "="^60)
        println("Phase 3a: Surface albedo Jacobian test")
        println("="^60)
        
        # Use the Rayleigh-only configuration (simpler, no absorption/aerosol)
        params = parameters_from_yaml(joinpath(@__DIR__, "test_parameters", "JacobianTestRayleigh.yaml"))
        
        # For Rayleigh-only: NAer=0, NGas=0, NSurf=1
        # But lin_model_from_parameters requires absorption_params to exist
        # Use the full config instead
        params = parameters_from_yaml(joinpath(@__DIR__, "test_parameters", "JacobianTest.yaml"))
        
        # Run linearized
        R_lin, T_lin, Ṙ_lin, Ṫ_lin, model, lin_model, NAer, NGas, NSurf = run_linearized(params)
        Nparams = NAer*7 + NGas + NSurf
        
        # Surface albedo is the last parameter
        iparam_surf = Nparams  # = NAer*7 + NGas + 1
        Ṙ_albedo_analytic = Ṙ_lin[iparam_surf, :, :, :]
        println("  Analytic dR/d(albedo) shape: ", size(Ṙ_albedo_analytic))
        println("  Analytic dR/d(albedo)[1,1,1] = ", Ṙ_albedo_analytic[1,1,1])
        
        # Finite difference
        ε = 1e-4
        function perturb_albedo!(p, δ)
            p.brdf[1] = LambertianSurfaceScalar{Float64}(p.brdf[1].albedo + δ)
        end
        
        println("  Computing finite-difference Jacobian (ε=$ε)...")
        Ṙ_albedo_fd = finite_difference_jacobian(params, perturb_albedo!, ε)
        println("  FD dR/d(albedo)[1,1,1] = ", Ṙ_albedo_fd[1,1,1])
        
        # Compare
        rel_err = relative_error(Ṙ_albedo_analytic, Ṙ_albedo_fd)
        mean_err = mean_relative_error(Ṙ_albedo_analytic, Ṙ_albedo_fd)
        println("  Max relative error: ", @sprintf("%.2e", rel_err))
        println("  Mean relative error: ", @sprintf("%.2e", mean_err))
        
        @test rel_err < 0.05  # 5% tolerance
        if rel_err < 0.05
            println("  ✓ Surface albedo Jacobian PASSED")
        else
            println("  ✗ Surface albedo Jacobian FAILED (rel_err = $rel_err)")
        end
    end
    
    # ── Phase 3b: Aerosol τ_ref Jacobian ──────────────────────────
    @testset "Phase 3b: Aerosol τ_ref Jacobian" begin
        println("\n" * "="^60)
        println("Phase 3b: Aerosol τ_ref Jacobian test")
        println("="^60)
        
        params = parameters_from_yaml(joinpath(@__DIR__, "test_parameters", "JacobianTest.yaml"))
        R_lin, T_lin, Ṙ_lin, Ṫ_lin, model, lin_model, NAer, NGas, NSurf = run_linearized(params)
        
        # τ_ref is aerosol parameter 1 (index 1 in Ṙ)
        iparam_τref = 1
        Ṙ_τref_analytic = Ṙ_lin[iparam_τref, :, :, :]
        println("  Analytic dR/d(τ_ref) shape: ", size(Ṙ_τref_analytic))
        println("  Analytic dR/d(τ_ref)[1,1,1] = ", Ṙ_τref_analytic[1,1,1])
        
        # Finite difference for τ_ref
        ε = 1e-4
        function perturb_τref!(p, δ)
            p.scattering_params.rt_aerosols[1].τ_ref += δ
        end
        
        println("  Computing finite-difference Jacobian (ε=$ε)...")
        Ṙ_τref_fd = finite_difference_jacobian(params, perturb_τref!, ε)
        println("  FD dR/d(τ_ref)[1,1,1] = ", Ṙ_τref_fd[1,1,1])
        
        rel_err = relative_error(Ṙ_τref_analytic, Ṙ_τref_fd)
        mean_err = mean_relative_error(Ṙ_τref_analytic, Ṙ_τref_fd)
        println("  Max relative error: ", @sprintf("%.2e", rel_err))
        println("  Mean relative error: ", @sprintf("%.2e", mean_err))
        
        @test rel_err < 0.05
        if rel_err < 0.05
            println("  ✓ Aerosol τ_ref Jacobian PASSED")
        else
            println("  ✗ Aerosol τ_ref Jacobian FAILED (rel_err = $rel_err)")
        end
    end
    
    # ── Phase 3c: Aerosol Mie parameters (nᵣ, nᵢ, rₚ, σₚ) ──────
    @testset "Phase 3c: Aerosol Mie parameter Jacobians" begin
        println("\n" * "="^60)
        println("Phase 3c: Aerosol Mie parameter Jacobian tests")
        println("="^60)
        
        params = parameters_from_yaml(joinpath(@__DIR__, "test_parameters", "JacobianTest.yaml"))
        R_lin, T_lin, Ṙ_lin, Ṫ_lin, model, lin_model, NAer, NGas, NSurf = run_linearized(params)
        
        # Aerosol parameter ordering: [τ_ref, nᵣ, nᵢ, rₚ, σₚ, p₀, σp]
        #                               1      2   3   4   5   6   7
        
        mie_params = [
            (2, "nᵣ", 1e-3, (p, δ) -> (p.scattering_params.rt_aerosols[1].aerosol.nᵣ += δ)),
            (3, "nᵢ", 1e-9, (p, δ) -> (p.scattering_params.rt_aerosols[1].aerosol.nᵢ += δ)),
            (4, "rₚ(μ)", 1e-4, (p, δ) -> begin
                dist = p.scattering_params.rt_aerosols[1].aerosol.size_distribution
                p.scattering_params.rt_aerosols[1].aerosol.size_distribution = LogNormal(dist.μ + δ, dist.σ)
            end),
            (5, "σₚ(σ)", 1e-3, (p, δ) -> begin
                dist = p.scattering_params.rt_aerosols[1].aerosol.size_distribution
                p.scattering_params.rt_aerosols[1].aerosol.size_distribution = LogNormal(dist.μ, dist.σ + δ)
            end),
        ]
        
        for (iparam, name, ε, perturb!) in mie_params
            println("\n  Testing dR/d($name)...")
            Ṙ_analytic = Ṙ_lin[iparam, :, :, :]
            println("    Analytic dR/d($name)[1,1,1] = ", Ṙ_analytic[1,1,1])
            
            println("    Computing FD (ε=$ε)...")
            Ṙ_fd = finite_difference_jacobian(params, perturb!, ε)
            println("    FD dR/d($name)[1,1,1] = ", Ṙ_fd[1,1,1])
            
            rel_err = relative_error(Ṙ_analytic, Ṙ_fd)
            mean_err = mean_relative_error(Ṙ_analytic, Ṙ_fd)
            println("    Max relative error: ", @sprintf("%.2e", rel_err))
            println("    Mean relative error: ", @sprintf("%.2e", mean_err))
            
            @test rel_err < 0.10  # 10% tolerance for Mie (complex chain)
            if rel_err < 0.10
                println("    ✓ $name Jacobian PASSED")
            else
                println("    ✗ $name Jacobian FAILED (rel_err = $rel_err)")
            end
        end
    end
    
    # ── Phase 3d: Aerosol vertical profile (p₀, σp) ──────────────
    @testset "Phase 3d: Aerosol vertical profile Jacobians" begin
        println("\n" * "="^60)
        println("Phase 3d: Aerosol vertical profile Jacobian tests")
        println("="^60)
        
        params = parameters_from_yaml(joinpath(@__DIR__, "test_parameters", "JacobianTest.yaml"))
        R_lin, T_lin, Ṙ_lin, Ṫ_lin, model, lin_model, NAer, NGas, NSurf = run_linearized(params)
        
        profile_params = [
            (6, "p₀", 1.0, (p, δ) -> begin
                dist = p.scattering_params.rt_aerosols[1].profile
                p.scattering_params.rt_aerosols[1].profile = Normal(mean(dist) + δ, std(dist))
            end),
            (7, "σp", 1.0, (p, δ) -> begin
                dist = p.scattering_params.rt_aerosols[1].profile
                p.scattering_params.rt_aerosols[1].profile = Normal(mean(dist), std(dist) + δ)
            end),
        ]
        
        for (iparam, name, ε, perturb!) in profile_params
            println("\n  Testing dR/d($name)...")
            Ṙ_analytic = Ṙ_lin[iparam, :, :, :]
            println("    Analytic dR/d($name)[1,1,1] = ", Ṙ_analytic[1,1,1])
            
            println("    Computing FD (ε=$ε)...")
            Ṙ_fd = finite_difference_jacobian(params, perturb!, ε)
            println("    FD dR/d($name)[1,1,1] = ", Ṙ_fd[1,1,1])
            
            rel_err = relative_error(Ṙ_analytic, Ṙ_fd)
            mean_err = mean_relative_error(Ṙ_analytic, Ṙ_fd)
            println("    Max relative error: ", @sprintf("%.2e", rel_err))
            println("    Mean relative error: ", @sprintf("%.2e", mean_err))
            
            @test rel_err < 0.10  # 10% tolerance
            if rel_err < 0.10
                println("    ✓ $name Jacobian PASSED")
            else
                println("    ✗ $name Jacobian FAILED (rel_err = $rel_err)")
            end
        end
    end
    
    # ── Phase 3e: Gas VMR Jacobian ────────────────────────────────
    @testset "Phase 3e: Gas VMR Jacobian" begin
        println("\n" * "="^60)
        println("Phase 3e: Gas VMR Jacobian test")
        println("="^60)
        
        # Create a version of the YAML with O2 as variable molecule
        params = parameters_from_yaml(joinpath(@__DIR__, "test_parameters", "JacobianTest.yaml"))
        
        # Check if variable_molecules is empty; if so, set O2 as variable
        if all(isempty.(params.absorption_params.variable_molecules))
            println("  NOTE: No variable molecules in YAML. Setting O2 as variable for this test.")
            # We need to reload with variable_molecules set
            # Modify the absorption params to have O2 as variable
            params.absorption_params.variable_molecules[1] = ["O2"]
            params.absorption_params.fixed_molecules[1] = String[]
        end
        
        # Run linearized with the variable molecule
        R_lin, T_lin, Ṙ_lin, Ṫ_lin, model, lin_model, NAer, NGas, NSurf = run_linearized(params)
        Nparams = NAer*7 + NGas + NSurf
        println("  NAer=$NAer, NGas=$NGas, NSurf=$NSurf, Nparams=$Nparams")
        
        # Gas VMR Jacobian is at index NAer*7 + 1 (H2O) or NAer*7 + 2 (first variable molecule)
        # The first slot (NAer*7+1) is always H2O. Variable gases start at NAer*7+2.
        # If O2 is the first (and only) variable gas, its index is NAer*7 + 2.
        # But we need to check: NGas = 1 + N_var_gas. If N_var_gas=1 (O2), NGas=2.
        iparam_gas = NAer*7 + NGas  # Last gas parameter (O2, since it's the only variable)
        println("  Gas VMR Jacobian index: $iparam_gas")
        
        if iparam_gas <= Nparams && NGas >= 2
            Ṙ_vmr_analytic = Ṙ_lin[iparam_gas, :, :, :]
            println("  Analytic dR/d(vmr_O2)[1,1,1] = ", Ṙ_vmr_analytic[1,1,1])
            
            ε = 0.21 * 1e-4  # 0.01% of O2 VMR
            function perturb_vmr_O2!(p, δ)
                p.absorption_params.vmr["O2"] += δ
            end
            
            println("  Computing FD (ε=$ε)...")
            Ṙ_vmr_fd = finite_difference_jacobian(params, perturb_vmr_O2!, ε)
            println("  FD dR/d(vmr_O2)[1,1,1] = ", Ṙ_vmr_fd[1,1,1])
            
            rel_err = relative_error(Ṙ_vmr_analytic, Ṙ_vmr_fd)
            mean_err = mean_relative_error(Ṙ_vmr_analytic, Ṙ_vmr_fd)
            println("  Max relative error: ", @sprintf("%.2e", rel_err))
            println("  Mean relative error: ", @sprintf("%.2e", mean_err))
            
            @test rel_err < 0.10
            if rel_err < 0.10
                println("  ✓ Gas VMR Jacobian PASSED")
            else
                println("  ✗ Gas VMR Jacobian FAILED (rel_err = $rel_err)")
            end
        else
            println("  SKIPPED: NGas=$NGas indicates no variable molecules in linearized model.")
            println("  (This test requires reloading parameters with variable_molecules in the YAML)")
            @test_skip true
        end
    end
end

println("\n" * "="^60)
println("All Jacobian tests complete!")
println("="^60)
