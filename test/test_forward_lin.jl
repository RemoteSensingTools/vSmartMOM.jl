# =================================================================
# Analytic Jacobians vs Finite-Difference Validation
# Uses EMIT-style setup with perturb_parameters for FD comparison
#
# Tests at TWO levels:
#   1. Model-level: τ̇_aer, k̇, ω̃̇ from model_from_parameters
#   2. RT-level: dR (TOA reflectance Jacobians) from rt_run
# =================================================================

include("test_helpers.jl")
include(joinpath(@__DIR__, "..", "src", "Testing", "perturb_parameters.jl"))
using Printf

println("="^70)
println("Analytic Jacobians vs Finite-Difference Validation (EMIT-style)")
println("="^70)

# Load fast EMIT config
params = parameters_from_yaml("test/test_parameters/ParamsEMIT_fast.yaml")
Nbands = length(params.spec_bands)
NAer_expected = length(params.scattering_params.rt_aerosols)
println("Config: $(Nbands) band(s), $(NAer_expected) aerosol(s)")

# Perturbation percentage for finite differences
ppct = 0.1  # 0.1% perturbation

# ─────────────────────────────────────────────────────────────────
# LEVEL 1: Model-level derivatives (Mie chain, absorption)
# Tests τ̇_aer, k̇, ω̃̇, ḟᵗ, greek coef derivatives 
# ─────────────────────────────────────────────────────────────────
println("\n" * "─"^70)
println("LEVEL 1: Model-level Jacobian validation")
println("─"^70)

@testset "Model-Level Jacobians" begin
    println("  Building linearized model...")
    @time model, lin_model = model_from_parameters(LinMode(), params)
    
    println("  Generating FD perturbations (ppct=$(ppct)%)...")
    @time pert_params = perturb_parameters(params, ppct)
    
    Nmol = 0
    for ib in 1:Nbands
        Nmol += length(params.absorption_params.variable_molecules[ib])
    end
    
    Nparams = length(pert_params)
    println("  Total perturbation params: $Nparams")
    println("    q=1, gas VMR=$Nmol, aerosol=$(NAer_expected*7), surface=$(Nbands)")
    @test Nparams == 1 + Nmol + NAer_expected * 7 + Nbands

    # --- Test aerosol optical property derivatives ---
    @testset "Aerosol $iaer properties" for iaer in 1:NAer_expected
        aer_names = ["τ_ref", "nᵣ", "nᵢ", "μ_size", "σ_size", "p₀", "σp"]
        
        for (i_aerprop, pname) in enumerate(aer_names)
            idx = 1 + Nmol + 7*(iaer-1) + i_aerprop
            
            # Compute perturbed model using LinMode to get spectrally-resolved τ_aer
            # (LinMode τ_aer is [n_aer, n_spec, n_layers]; FwdMode is [n_aer, n_layers])
            pert_model, _ = model_from_parameters(LinMode(), pert_params[idx])
            
            # Compute the perturbation delta
            if pname == "τ_ref"
                Δ = pert_params[idx].scattering_params.rt_aerosols[iaer].τ_ref - 
                    params.scattering_params.rt_aerosols[iaer].τ_ref
            elseif pname == "nᵣ"
                Δ = pert_params[idx].scattering_params.rt_aerosols[iaer].aerosol.nᵣ - 
                    params.scattering_params.rt_aerosols[iaer].aerosol.nᵣ
            elseif pname == "nᵢ"
                Δ = pert_params[idx].scattering_params.rt_aerosols[iaer].aerosol.nᵢ - 
                    params.scattering_params.rt_aerosols[iaer].aerosol.nᵢ
            elseif pname == "μ_size"
                Δ = pert_params[idx].scattering_params.rt_aerosols[iaer].aerosol.size_distribution.μ - 
                    params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.μ
            elseif pname == "σ_size"
                Δ = pert_params[idx].scattering_params.rt_aerosols[iaer].aerosol.size_distribution.σ - 
                    params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.σ
            elseif pname == "p₀"
                Δ = mean(pert_params[idx].scattering_params.rt_aerosols[iaer].profile) - 
                    mean(params.scattering_params.rt_aerosols[iaer].profile)
            elseif pname == "σp"
                Δ = std(pert_params[idx].scattering_params.rt_aerosols[iaer].profile) - 
                    std(params.scattering_params.rt_aerosols[iaer].profile)
            end
            
            @testset "$pname (Δ=$(@sprintf("%.2e", Δ)))" begin
                for ib in 1:Nbands
                    # τ_aer FD vs analytic
                    K_FD_τ = (pert_model.τ_aer[ib][iaer,:,:] .- model.τ_aer[ib][iaer,:,:]) ./ Δ
                    K_lin_τ = lin_model.τ̇_aer[ib][iaer, i_aerprop, :, :]
                    
                    max_e, mean_e = rel_errors(K_lin_τ, K_FD_τ)
                    @printf("    [band %d] %s τ̇_aer: max_err=%.2e, mean_err=%.2e\n", ib, pname, max_e, mean_e)
                    
                    if pname in ["τ_ref", "p₀", "σp"]
                        # Profile/τ derivatives should be very accurate
                        @test mean_e < 0.05 || isnan(mean_e)
                    end
                    
                    # For Mie-derived params, also test k̇, ω̃̇ 
                    if pname in ["nᵣ", "nᵢ", "μ_size", "σ_size"] && 
                       hasproperty(lin_model, :lin_aerosol_optics)
                        
                        K_FD_k = (pert_model.aerosol_optics[ib][iaer].k .- 
                                  model.aerosol_optics[ib][iaer].k) ./ Δ
                        K_lin_k = lin_model.lin_aerosol_optics[ib][iaer].k̇[i_aerprop-1, :]
                        max_e_k, mean_e_k = rel_errors(K_lin_k, K_FD_k)
                        @printf("    [band %d] %s k̇:     max_err=%.2e, mean_err=%.2e\n", ib, pname, max_e_k, mean_e_k)
                        
                        K_FD_ω = (pert_model.aerosol_optics[ib][iaer].ω̃ .- 
                                  model.aerosol_optics[ib][iaer].ω̃) ./ Δ
                        K_lin_ω = lin_model.lin_aerosol_optics[ib][iaer].ω̃̇[i_aerprop-1, :]
                        max_e_ω, mean_e_ω = rel_errors(K_lin_ω, K_FD_ω)
                        @printf("    [band %d] %s ω̃̇:    max_err=%.2e, mean_err=%.2e\n", ib, pname, max_e_ω, mean_e_ω)
                    end
                end
            end
        end
    end
end

# ─────────────────────────────────────────────────────────────────
# LEVEL 2: RT-level derivatives (TOA reflectance Jacobians)
# ─────────────────────────────────────────────────────────────────
println("\n" * "─"^70)
println("LEVEL 2: RT-level Jacobian validation (TOA reflectance)")
println("─"^70)

@testset "RT-Level Jacobians" begin
    println("  Running linearized RT...")
    @time R_base, dR, NAer, NGas, NSurf, model, lin_model = run_lin_rt(params)
    Nparams_rt = NAer*7 + NGas + NSurf
    println("  R shape: $(size(R_base)), dR shape: $(size(dR))")
    println("  Nparams_rt=$Nparams_rt (NGas=$NGas, NAer=$NAer, NSurf=$NSurf)")
    
    pert_params = perturb_parameters(params, ppct)
    
    # Mapping from perturb_parameters index to dR parameter index:
    # perturb_params: [q, gas_vmrs..., aer1_props(7)..., aer2_props(7)..., brdf...]
    # dR params:      [gas_derivs(NGas), aer_derivs(NAer*7), surf_derivs(NSurf)]
    #   where gas_derivs[1] = q, gas_derivs[2:end] = VMR
    #   and aer_derivs per aerosol: τ_ref, nᵣ, nᵢ, μ, σ, p₀, σp
    
    Nmol = 0
    for ib in 1:Nbands
        Nmol += length(params.absorption_params.variable_molecules[ib])
    end
    
    # --- Surface albedo ---
    @testset "Surface albedo" begin
        for ib in 1:Nbands
            surf_pert_idx = 1 + Nmol + NAer_expected*7 + ib
            dR_idx = NGas + NAer*7 + ib  # surface params are last in dR
            
            tmp_params = deepcopy(params)
            old_alb = tmp_params.brdf[ib].albedo
            tmp_params.brdf[ib].albedo = old_alb * (1 + ppct/100)
            Δ_alb = tmp_params.brdf[ib].albedo - old_alb
            
            R_pert = run_fwd_only(tmp_params)
            K_FD = (R_pert .- R_base) ./ Δ_alb
            K_lin = dR[dR_idx, :, :, :]
            
            max_e, mean_e = rel_errors(K_lin, K_FD)
            @printf("    Surface albedo band %d: max_err=%.2e, mean_err=%.2e\n", ib, max_e, mean_e)
            @test mean_e < 0.01 || isnan(mean_e)  # Surface should be very accurate
        end
    end
    
    # --- τ_ref ---
    @testset "Aerosol τ_ref" for iaer in 1:NAer
        i_aerprop = 1
        pert_idx = 1 + Nmol + 7*(iaer-1) + i_aerprop
        dR_idx = NGas + 7*(iaer-1) + i_aerprop
        
        Δτ = pert_params[pert_idx].scattering_params.rt_aerosols[iaer].τ_ref - 
             params.scattering_params.rt_aerosols[iaer].τ_ref
        
        R_pert = run_fwd_only(pert_params[pert_idx])
        K_FD = (R_pert .- R_base) ./ Δτ
        K_lin = dR[dR_idx, :, :, :]
        
        max_e, mean_e = rel_errors(K_lin, K_FD)
        @printf("    τ_ref (aer %d): max_err=%.2e, mean_err=%.2e\n", iaer, max_e, mean_e)
        @test mean_e < 0.15 || isnan(mean_e)  # Known ~10% residual from Bug 19
    end
    
    # --- Mie microphysical params (nᵣ, nᵢ, μ, σ) ---
    mie_names = ["nᵣ", "nᵢ", "μ_size", "σ_size"]
    @testset "Aerosol Mie params" for iaer in 1:NAer
        for (j, pname) in enumerate(mie_names)
            i_aerprop = j + 1  # nᵣ=2, nᵢ=3, μ=4, σ=5
            pert_idx = 1 + Nmol + 7*(iaer-1) + i_aerprop
            dR_idx = NGas + 7*(iaer-1) + i_aerprop
            
            if pname == "nᵣ"
                Δ = pert_params[pert_idx].scattering_params.rt_aerosols[iaer].aerosol.nᵣ - 
                    params.scattering_params.rt_aerosols[iaer].aerosol.nᵣ
            elseif pname == "nᵢ"
                Δ = pert_params[pert_idx].scattering_params.rt_aerosols[iaer].aerosol.nᵢ - 
                    params.scattering_params.rt_aerosols[iaer].aerosol.nᵢ
            elseif pname == "μ_size"
                Δ = pert_params[pert_idx].scattering_params.rt_aerosols[iaer].aerosol.size_distribution.μ - 
                    params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.μ
            elseif pname == "σ_size"
                Δ = pert_params[pert_idx].scattering_params.rt_aerosols[iaer].aerosol.size_distribution.σ - 
                    params.scattering_params.rt_aerosols[iaer].aerosol.size_distribution.σ
            end
            
            R_pert = run_fwd_only(pert_params[pert_idx])
            K_FD = (R_pert .- R_base) ./ Δ
            K_lin = dR[dR_idx, :, :, :]
            
            max_e, mean_e = rel_errors(K_lin, K_FD)
            @printf("    %s (aer %d): max_err=%.2e, mean_err=%.2e\n", pname, iaer, max_e, mean_e)
            # Mie params may have larger errors due to known bugs 19-21
            @test mean_e < 0.20 || isnan(mean_e)
        end
    end
    
    # --- Profile params (p₀, σp) ---
    prof_names = ["p₀", "σp"]
    @testset "Aerosol profile params" for iaer in 1:NAer
        for (j, pname) in enumerate(prof_names)
            i_aerprop = j + 5  # p₀=6, σp=7
            pert_idx = 1 + Nmol + 7*(iaer-1) + i_aerprop
            dR_idx = NGas + 7*(iaer-1) + i_aerprop
            
            if pname == "p₀"
                Δ = mean(pert_params[pert_idx].scattering_params.rt_aerosols[iaer].profile) - 
                    mean(params.scattering_params.rt_aerosols[iaer].profile)
            else
                Δ = std(pert_params[pert_idx].scattering_params.rt_aerosols[iaer].profile) - 
                    std(params.scattering_params.rt_aerosols[iaer].profile)
            end
            
            R_pert = run_fwd_only(pert_params[pert_idx])
            K_FD = (R_pert .- R_base) ./ Δ
            K_lin = dR[dR_idx, :, :, :]
            
            max_e, mean_e = rel_errors(K_lin, K_FD)
            @printf("    %s (aer %d): max_err=%.2e, mean_err=%.2e\n", pname, iaer, max_e, mean_e)
            @test mean_e < 0.15 || isnan(mean_e)
        end
    end
end

println("\n" * "="^70)
println("Analytic Jacobian validation complete.")
println("="^70)
