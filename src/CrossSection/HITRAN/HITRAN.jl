
module HITRAN

using Pandas

export readHITRAN

varNames = ["molec_id", "local_iso_id", "nu", "sw", "a", "gamma_air",
            "gamma_self", "elower", "n_air", "delta_air", "global_upper_quanta",
            "global_lower_quanta", "local_upper_quanta", "local_lower_quanta",
            "ierr", "iref", "line_mixing_flag", "gp", "gpp"]
varLengths = [2, 1, 12, 10, 10, 5, 5, 10, 4, 8, 15, 15, 15, 15, 6, 12, 1, 7, 7]
idxRanges = append!([0], [sum(varLengths[1:i]) for i in 1:length(varLengths)])


function readHITRAN(filename, molecule, isotope, ν_min, ν_max)

    open(filename) do file
        allDataRows = []
        for ln in eachline(file)
            parts = [tryparse(Float64, replace(ln[(idxRanges[i]+1):idxRanges[i+1]], " " => "")) for i in 1:(length(varLengths))]
            if(parts[1] == molecule && (parts[2] == isotope || isotope == -1) && (ν_min <= parts[3] <= ν_max))
                allDataRows = append!(allDataRows, [parts])
            end
        end

        allDataDict = Dict()

        for j in 1:length(varNames)
            parts = [allDataRows[i][j] for i in 1:length(allDataRows)]
            allDataDict[varNames[j]] = parts
        end

        allDataDF = DataFrame(allDataDict)

        return allDataDF
    end

end

end
