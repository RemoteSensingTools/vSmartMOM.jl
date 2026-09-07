# Synthetic inputs for portable workflow tests; no campaign data are used.
const TEST_SURFACES = (:urban, :rural, :desert, :forest)

function corrected_sif_provenance()
    radiance = 0.5 / (2π)
    return Dict{String,Any}(
        "sif_definition_version" => Int32(2),
        "sif_definition" =>
            "isotropic BOA radiance normalized by 2pi*L_lambda(760 nm)=0.5",
        "sif_case_on_label" => "angular_integral760_0p5",
        "sif_reference_wavelength_nm" => 760.0,
        "sif_upwelling_solid_angle_sr" => 2π,
        "sif_angular_integral_760_mW_m-2_nm-1" => 0.5,
        "sif_radiance_760_mW_m-2_sr-1_nm-1" => radiance,
        "sif_cosine_weighted_irradiance_760_mW_m-2_nm-1" => π * radiance,
        "sif_SIF760_mW_m-2_sr-1_per_cm-1" =>
            radiance * 760.0^2 / 1e7,
        "sif_mSIF_mW_m-2_sr-1_per_cm-2" => 1.2291230681458325e-5,
        "sif_template_wavelength_integral_mW_m-2_sr-1" =>
            15.368806005166872,
    )
end

function write_test_measurement(path, truth; provenance=Dict{String,Any}())
    mkpath(dirname(path))
    NCDataset(path, "c") do dataset
        dataset.attrib["instrument_processing_complete"] = 1
        dataset.attrib["state_index"] = truth.state_index
        dataset.attrib["sif_case"] = String(truth.sif_case)
        for (key, value) in provenance
            dataset.attrib[key] = value
        end
    end
    return path
end

function write_test_noise(path, truth; provenance=Dict{String,Any}())
    mkpath(dirname(path))
    NCDataset(path, "c") do dataset
        defDim(dataset, "measurement", 3)
        defDim(dataset, "band", 3)
        for class in ("corrected", "uncorrected")
            defVar(dataset, "measurement_$class", Float64,
                   ("measurement",))[:] = [1.0, 2.0, 3.0]
            defVar(dataset, "noise_std_$class", Float64,
                   ("measurement",))[:] = fill(0.1, 3)
            defVar(dataset, "Se_diagonal_$class", Float64,
                   ("measurement",))[:] = fill(0.01, 3)
        end
        defVar(dataset, "wavelength", Float64,
               ("measurement",))[:] = [760.0, 1600.0, 2050.0]
        defVar(dataset, "band_start_index", Int32,
               ("band",))[:] = Int32[1, 2, 3]
        defVar(dataset, "band_end_index", Int32,
               ("band",))[:] = Int32[1, 2, 3]
        dataset.attrib["noise_covariance_complete"] = 1
        dataset.attrib["state_index"] = truth.state_index
        dataset.attrib["sif_case"] = String(truth.sif_case)
        dataset.attrib["campaign"] = String(truth.campaign)
        dataset.attrib["background_co2_ppm"] = truth.background_co2_ppm
        dataset.attrib["bottom_co2_layer_index"] = truth.bottom_layer_index
        dataset.attrib["bottom_co2_ppm"] = truth.bottom_co2_ppm
        dataset.attrib["xco2_ppm"] = truth.xco2_ppm
        for (key, value) in provenance
            dataset.attrib[key] = value
        end
    end
    return path
end

function write_test_truth_table(path; bottom_layer=false)
    mkpath(dirname(path))
    open(path, "w") do io
        if bottom_layer
            println(io, "# index surface_index surface aerosol_index aerosol_case " *
                        "sif_index sif_case bottom_co2_index background_co2_ppm " *
                        "bottom_layer_index bottom_co2_ppm xco2_ppm")
        else
            println(io, "# index surface_index surface aerosol_index aerosol_case " *
                        "sif_index sif_case xco2_index xco2_ppm")
        end
        index = 0
        co2_values = bottom_layer ? (360, 380, 400, 420, 440) :
                                    (380, 400, 420, 440)
        for (surface_index, surface) in enumerate(TEST_SURFACES),
                (aerosol_index, aerosol) in enumerate((:none, :aod760_0p28)),
                (sif_index, sif) in enumerate((
                    :off, :angular_integral760_0p5)),
                (co2_index, co2) in enumerate(co2_values)
            index += 1
            if bottom_layer
                xco2 = 400 + 0.06229780735737046 * (co2 - 400)
                println(io, "$index $surface_index $surface $aerosol_index " *
                            "$aerosol $sif_index $sif $co2_index 400 16 " *
                            "$co2 $xco2")
            else
                println(io, "$index $surface_index $surface $aerosol_index " *
                            "$aerosol $sif_index $sif $co2_index $co2")
            end
        end
    end
    return path
end

