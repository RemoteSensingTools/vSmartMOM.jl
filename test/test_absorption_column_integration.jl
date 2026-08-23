using Test
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Absorption

function _constant_profile(p_half::Vector{Float64}; T=250.0, q=0.0,
                           vmr=Dict{String,Any}())
    temperatures = fill(T, length(p_half) - 1)
    humidity = fill(q, length(p_half) - 1)
    profile, _ = CoreRT.prepare_observer_profile(
        temperatures, p_half, humidity, vmr, 0.0, -1)
    return profile
end

function _write_cia_block(io, formula, T, ν, σ, reference_code)
    @assert length(ν) == length(σ)
    header = rpad(formula, 20) *
             lpad(string(first(ν)), 10) *
             lpad(string(last(ν)), 10) *
             lpad(string(length(ν)), 7) *
             lpad(string(T), 7) *
             lpad("1.000e-40", 10) *
             repeat(" ", 33) *
             lpad(reference_code, 3)
    println(io, header)
    for i in eachindex(ν)
        println(io, "$(ν[i]) $(σ[i])")
    end
end

@testset "Column-exact CIA and MT_CKD integration" begin
    # Adding HITRAN reference metadata must not remove the inferred public
    # constructor that existed when CIABlock had only four fields.
    inferred_block = Absorption.CIABlock(
        SubString("O2-O2", 1, 5), 250.0f0, [1.0, 2.0], Float32[3, 4])
    @test inferred_block isa Absorption.CIABlock{Float64}
    @test inferred_block.reference_code == ""

    # Missing spectral coverage at one temperature is not a measured zero.
    # Interpolate only among panels that cover each requested wavenumber.
    coverage_blocks = [
        Absorption.CIABlock{Float64}(
            "O2-O2", 200.0, [1.0, 2.0], [10.0, 20.0]),
        Absorption.CIABlock{Float64}(
            "O2-O2", 300.0, [2.0, 3.0], [30.0, 40.0]),
    ]
    coverage_table = Absorption.build_cia_table(
        coverage_blocks, [1.0, 2.0, 3.0, 4.0])
    coverage_sigma = zeros(Float64, 4)
    Absorption.cia_σ_at_T!(coverage_sigma, coverage_table, 250.0)
    @test coverage_sigma == [10.0, 25.0, 40.0, 0.0]

    # Where same-temperature blocks overlap, retain the finer source grid.
    overlap_blocks = [
        Absorption.CIABlock{Float64}(
            "O2-O2", 250.0, [1.0, 3.0], [100.0, 300.0]),
        Absorption.CIABlock{Float64}(
            "O2-O2", 250.0, [1.0, 2.0, 3.0], [100.0, 999.0, 300.0]),
    ]
    overlap_table = Absorption.build_cia_table(overlap_blocks, [2.0])
    overlap_sigma = zeros(Float64, 1)
    Absorption.cia_σ_at_T!(overlap_sigma, overlap_table, 250.0)
    @test only(overlap_sigma) == 999.0

    # Same isothermal 1-bar column represented by five versus twenty layers.
    # Midpoint n²*Δz changes with this subdivision; n(mid)*VCD does not.
    coarse = _constant_profile(collect(range(0.01, 1000.0; length=6));
                               vmr=Dict{String,Any}("O2" => 0.20946))
    fine = _constant_profile(collect(range(0.01, 1000.0; length=21));
                             vmr=Dict{String,Any}("O2" => 0.20946))

    cia = Absorption.CIATable("O2-O2", "O2", "O2",
                              [200.0, 300.0], fill(1e-43, 1, 2))
    τ_cia_coarse = zeros(Float64, 1, length(coarse.T))
    τ_cia_fine = zeros(Float64, 1, length(fine.T))
    Absorption.compute_τ_cia!(τ_cia_coarse, cia, coarse, coarse.vmr)
    Absorption.compute_τ_cia!(τ_cia_fine, cia, fine, fine.vmr)
    @test sum(τ_cia_coarse) ≈ sum(τ_cia_fine) rtol=2e-12

    # Mixed dry/H2O colliders use the profile's H2O density and exact H2O
    # column. The symmetric midpoint-column form must be invariant both to
    # layer subdivision and to the order used in the CIA pair label.
    coarse_cia_wet = _constant_profile(
        collect(range(0.01, 1000.0; length=6)); q=0.01,
        vmr=Dict{String,Any}("O2" => 0.20946))
    fine_cia_wet = _constant_profile(
        collect(range(0.01, 1000.0; length=21)); q=0.01,
        vmr=Dict{String,Any}("O2" => 0.20946))
    cia_o2_h2o = Absorption.CIATable(
        "O2-H2O", "O2", "H2O", [200.0, 300.0], fill(1e-43, 1, 2))
    cia_h2o_o2 = Absorption.CIATable(
        "H2O-O2", "H2O", "O2", [200.0, 300.0], fill(1e-43, 1, 2))
    τ_mixed_coarse = zeros(Float64, 1, length(coarse_cia_wet.T))
    τ_mixed_fine = zeros(Float64, 1, length(fine_cia_wet.T))
    τ_mixed_reverse = zeros(Float64, 1, length(coarse_cia_wet.T))
    Absorption.compute_τ_cia!(
        τ_mixed_coarse, cia_o2_h2o, coarse_cia_wet, coarse_cia_wet.vmr)
    Absorption.compute_τ_cia!(
        τ_mixed_fine, cia_o2_h2o, fine_cia_wet, fine_cia_wet.vmr)
    Absorption.compute_τ_cia!(
        τ_mixed_reverse, cia_h2o_o2, coarse_cia_wet, coarse_cia_wet.vmr)
    @test sum(τ_mixed_coarse) ≈ sum(τ_mixed_fine) rtol=2e-12
    @test τ_mixed_coarse ≈ τ_mixed_reverse rtol=2e-15

    # Synthetic continuum coefficients are sufficient to test the vertical
    # integration invariant. A constant q gives a constant H2O VMR.
    coarse_wet = _constant_profile(collect(range(0.01, 1000.0; length=6));
                                   q=0.01)
    fine_wet = _constant_profile(collect(range(0.01, 1000.0; length=21));
                                 q=0.01)
    band = Absorption.MTCKDBand([1e-27], [2e-27], [0.0], 1013.0, 296.0)
    ν = [5000.0]
    τ_cont_coarse = zeros(Float64, 1, length(coarse_wet.T))
    τ_cont_fine = zeros(Float64, 1, length(fine_wet.T))
    Absorption.compute_τ_h2o_continuum!(
        τ_cont_coarse, band, ν, coarse_wet, coarse_wet.vmr_h2o)
    Absorption.compute_τ_h2o_continuum!(
        τ_cont_fine, band, ν, fine_wet, fine_wet.vmr_h2o)
    @test sum(τ_cont_coarse) ≈ sum(τ_cont_fine) rtol=2e-12

    # AtmosphericProfile stores H2O as the dry-air molar ratio r=N_H2O/N_dry,
    # not as a moist-air mole fraction. Its hydrostatic columns must close
    # exactly on that definition.
    @test coarse_wet.vcd_h2o ./ coarse_wet.vcd_dry ≈
          coarse_wet.vmr_h2o rtol=2e-15
    expected_continuum = similar(τ_cont_coarse)
    for iz in eachindex(coarse_wet.T)
        r = coarse_wet.vmr_h2o[iz]
        P = coarse_wet.p_full[iz]
        T = coarse_wet.T[iz]
        p_dry = P / (1 + r)
        p_h2o = r * p_dry
        radterm = only(ν) * tanh(1.4388 * only(ν) / (2T))
        sigma = (only(band.C_self) * p_h2o / band.p_ref +
                 only(band.C_for) * p_dry / band.p_ref) * radterm
        expected_continuum[1, iz] = sigma * coarse_wet.vcd_h2o[iz]
    end
    @test τ_cont_coarse ≈ expected_continuum rtol=2e-15

    # At fixed temperature, humidity, and composition, both CIA and MT_CKD
    # scale as midpoint density/pressure times hydrostatic column.  The
    # surface-pressure tangent therefore contains both derivatives.  Verify
    # that factor against a centered finite difference of the actual opacity
    # helpers rather than only against its algebraic reconstruction.
    p_base = collect(range(0.01, 1000.0; length=6))
    δp = 1e-2
    function pressure_perturbed_profile(offset)
        p_half = copy(p_base)
        p_half[end] += offset
        return _constant_profile(
            p_half; q=0.01, vmr=Dict{String,Any}("O2" => 0.20946))
    end
    profile_minus = pressure_perturbed_profile(-δp)
    profile_base = pressure_perturbed_profile(0.0)
    profile_plus = pressure_perturbed_profile(δp)
    psurf_tangent = CoreRT.psurf_profile_tangents(profile_base)
    dry_ratio_dot = psurf_tangent.vcd_dry_dot[end] /
                    profile_base.vcd_dry[end]
    binary_ratio_dot = dry_ratio_dot + 0.5 / profile_base.p_full[end]

    function bottom_cia(profile)
        τ = zeros(Float64, 1, length(profile.T))
        Absorption.compute_τ_cia!(τ, cia, profile, profile.vmr)
        return τ[1, end]
    end
    cia_fd = (bottom_cia(profile_plus) - bottom_cia(profile_minus)) / (2δp)
    @test cia_fd ≈ bottom_cia(profile_base) * binary_ratio_dot rtol=2e-10

    function bottom_continuum(profile)
        τ = zeros(Float64, 1, length(profile.T))
        Absorption.compute_τ_h2o_continuum!(
            τ, band, ν, profile, profile.vmr_h2o)
        return τ[1, end]
    end
    continuum_fd =
        (bottom_continuum(profile_plus) - bottom_continuum(profile_minus)) / (2δp)
    @test continuum_fd ≈
          bottom_continuum(profile_base) * binary_ratio_dot rtol=2e-10

    # Reference selection and negative handling flow from parsed CIA entries
    # through the same helper used by forward, linearized, and batch-update
    # model construction. A legacy bare path fails closed where two reference
    # families overlap on the requested model grid.
    mktempdir() do dir
        path = joinpath(dir, "X-X.cia")
        open(path, "w") do io
            _write_cia_block(io, "X-X", 250.0, [1.0, 2.0, 3.0],
                             [1.0, 1.0, 1.0], "11")
            _write_cia_block(io, "X-X", 250.0, [1.0, 2.0, 3.0],
                             [2.0, 2.0, 2.0], "22")
            _write_cia_block(io, "X-X", 250.0, [1.0, 2.0, 3.0],
                             [-1.0, -1.0, -1.0], "33")
        end

        parsed_blocks = Absorption.parse_cia_file(path)
        @test [block.reference_code for block in parsed_blocks] ==
              ["11", "22", "33"]
        @test_throws ArgumentError Absorption.load_cia_table(path, [2.0])
        @test_throws ArgumentError Absorption.load_cia_table(
            path, [2.0]; reference_codes=Int[22])

        selected = Absorption.load_cia_table(
            path, [2.0]; reference_codes="22")
        selected_sigma = zeros(Float64, 1)
        Absorption.cia_σ_at_T!(selected_sigma, selected, 250.0)
        @test only(selected_sigma) == 2.0
        @test_throws ErrorException Absorption.load_cia_table(
            path, [2.0]; reference_codes="33")
        clamped = Absorption.load_cia_table(
            path, [2.0]; reference_codes=["33"],
            negative_policy=:clamp_zero)
        clamped_sigma = zeros(Float64, 1)
        Absorption.cia_σ_at_T!(clamped_sigma, clamped, 250.0)
        @test only(clamped_sigma) == 0.0

        absorption_dict = Dict{String,Any}(
            "fixed_molecules" => [String[]],
            "variable_molecules" => [String[]],
            "vmr" => Dict{String,Any}(),
            "broadening" => "Voigt()",
            "CEF" => "HumlicekWeidemann32SDErrorFunction()",
            "wing_cutoff" => 25,
            "cia_files" => Any[Dict{String,Any}(
                "path" => path,
                "reference_codes" => "22",
                "negative_policy" => "error")],
        )
        ap = vSmartMOM.IO._parse_absorption(
            Dict{String,Any}("absorption" => absorption_dict), Float64)
        @test ap.cia_files == [path]
        @test ap.cia_reference_codes == [["22"]]
        @test ap.cia_negative_policies == [:error]

        configured = CoreRT._load_configured_cia_table(
            ap, 1, [2.0], Float64)
        configured_sigma = zeros(Float64, 1)
        Absorption.cia_σ_at_T!(configured_sigma, configured, 250.0)
        @test only(configured_sigma) == 2.0

        legacy_dict = deepcopy(absorption_dict)
        legacy_dict["cia_files"] = [path]
        legacy_ap = vSmartMOM.IO._parse_absorption(
            Dict{String,Any}("absorption" => legacy_dict), Float64)
        @test legacy_ap.cia_reference_codes == [nothing]
        @test legacy_ap.cia_negative_policies == [:error]
        @test_throws ArgumentError CoreRT._load_configured_cia_table(
            legacy_ap, 1, [2.0], Float64)

        bad_policy = deepcopy(absorption_dict)
        bad_policy["cia_files"][1]["negative_policy"] = "discard"
        @test_throws ArgumentError vSmartMOM.IO._parse_absorption(
            Dict{String,Any}("absorption" => bad_policy), Float64)
        bad_codes = deepcopy(absorption_dict)
        bad_codes["cia_files"][1]["reference_codes"] = String[]
        @test_throws ArgumentError vSmartMOM.IO._parse_absorption(
            Dict{String,Any}("absorption" => bad_codes), Float64)
    end

    # Production O2-O2 visible CIA explicitly uses HITRAN reference 54
    # (Finkenzeller & Volkamer 2022). Reference 32 overlaps but represents a
    # different source family and contains negative measurements; never stitch
    # the two implicitly.
    cia_path = get(ENV, "VSMARTMOM_TEST_O2_O2_CIA",
        "/home/sanghavi/data/HITRAN_CIA/O2-O2_2024.cia")
    if isfile(cia_path)
        selected = Absorption.load_cia_table(
            cia_path, 1e7 ./ [350.0, 400.0, 450.0];
            reference_codes=["54"])
        @test selected.reference_codes == ["54"]
        @test all(selected.coverage)
        sigma_293 = zeros(Float64, 3)
        Absorption.cia_σ_at_T!(sigma_293, selected, 293.0)
        @test all(>=(0), sigma_293)
        @test sigma_293[2] ≈ 7.3726e-50 rtol=1e-12

        ultraviolet = Absorption.load_cia_table(
            cia_path, 1e7 ./ [250.0, 296.0, 298.0, 350.0];
            reference_codes="54")
        @test ultraviolet.coverage == Bool[false, false, true, true]
        @test_throws ArgumentError Absorption.load_cia_table(
            cia_path, [1e7 / 368.48])
    else
        @test_skip isfile(cia_path)
    end

    mtckd_path = get(ENV, "VSMARTMOM_TEST_MTCKD",
        "/home/sanghavi/data/MT_CKD_4.2/absco-ref_wv-mt-ckd.nc")
    if isfile(mtckd_path)
        mtckd = Absorption.load_mtckd(mtckd_path)
        @test occursin("4.3", mtckd.title)
        @test occursin("MT_CKD_4.3", mtckd.version_description)
    else
        @test_skip isfile(mtckd_path)
    end
end
