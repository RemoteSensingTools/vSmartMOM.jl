

function readHITRAN(filename, molecule, isotope, ν_min, ν_max)

    varNames = ["molec_id", "local_iso_id", "nu", "sw", "a", "gamma_air",
                "gamma_self", "elower", "n_air", "delta_air"] #, "global_upper_quanta",
                # "global_lower_quanta", "local_upper_quanta", "local_lower_quanta",
                # "ierr", "iref", "line_mixing_flag", "gp", "gpp"]
    varLengths = [2, 1, 12, 10, 10, 5, 5, 10, 4, 8] #, 15, 15, 15, 15, 6, 12, 1, 7, 7]
    idxRanges = append!([0], [sum(varLengths[1:i]) for i in 1:length(varLengths)])

    open(filename) do file
        allDataRows = []
        for ln in eachline(file)
            parts = [tryparse(Float64, replace(ln[(idxRanges[i]+1):idxRanges[i+1]], " " => "")) for i in 1:(length(varLengths))]
            if(parts[1] == molecule && (parts[2] == isotope || isotope == -1) && (ν_min <= parts[3] <= ν_max))
                allDataRows = append!(allDataRows, [parts])
            end
        end

        parts = [[allDataRows[i][j] for i in 1:length(allDataRows)] for j in 1:length(varNames)]

        if(length(parts[1]) == 0)
            throw(HitranEmptyError())
        end

        finalStruct = HitranTable(mol=parts[1], iso=parts[2], νᵢ=parts[3], Sᵢ=parts[4], Aᵢ=parts[5], γ_air=parts[6], γ_self=parts[7], E″=parts[8], n_air=parts[9], δ_air=parts[10])

        return finalStruct

    end

end
