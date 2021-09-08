#####
##### Function to parse hitran par-file data
#####

using ..Architectures: GPU

"""
    read_hitran(filepath::String, mol::Int=-1, iso::Int=-1, ν_min::Real=0, ν_max::Real=Inf)

Read/parse a HITRAN data file and return the data in [`HitranTable`](@ref) format

"""
function read_hitran(filepath::String; mol::Int=-1, iso::Int=-1, 
                     ν_min::Real=0, ν_max::Real=Inf, 
                     min_strength::Real=0)

    # Infer type from input ν_min
    FT = eltype(AbstractFloat(ν_min))

    # Declare some constant properties of the hitran line parameters, such as
    # names, lengths, and types
    varNames = ["molec_id", "local_iso_id", "nu", "sw", "a", "gamma_air",
                "gamma_self", "elower", "n_air", "delta_air", "global_upper_quanta",
                "global_lower_quanta", "local_upper_quanta", "local_lower_quanta",
                "ierr", "iref", "line_mixing_flag", "gp", "gpp"]
    varLengths = [2, 1, 12, 10, 10, 5, 5, 10, 4, 8, 15, 15, 15, 15, 6, 12, 1, 7, 7]
    varTypes = [Int64, Int64, FT, FT, FT, FT, FT, FT, FT, FT, String, String, String, String, String, String, String, FT, FT]

    # Take the line parameters' lengths and assemble an index-ranges list
    # Each value is the starting index, minus 1, of var i and ending index of var i-1
    # i.e.
    #   [0, 2, 3, 15, 25, 35, 40, 45, 55, 59, 67, 82, 97, 112, 127, 133, 145, 146, 153, 160]
    # Starting at 0 correctly produces (1,2), (3,3), (4,15), ..., (154, 160)
    idxRanges = append!([0], [sum(varLengths[1:i]) for i in 1:length(varLengths)])

    # Open the file specified
    open(filepath, "r") do file

        # Where to hold the matching rows
        rows = []

        # Loop through every line
        for ln in eachline(file)

            # Go from a line to a list of values
            values = [varTypes[i] == String ? ln[(idxRanges[i]+1):idxRanges[i+1]] : something(tryparse(varTypes[i], ln[(idxRanges[i]+1):idxRanges[i+1]]), varTypes[i](0)) for i in 1:(length(varLengths))]

            # Check that the search criteria are met (molecule, isotope, wavenumber range, and min. line-strength)
            if((values[1] == mol || mol == -1) && (values[2] == iso || iso == -1) 
                && (ν_min <= values[3] <= ν_max) && values[4] >= min_strength)

                # Add this row to the list of rows
                rows = append!(rows, [values])
            end
        end

        # Convert the row-major list of lists to a column-major list of lists
        cols = [[rows[i][j] for i in 1:length(rows)] for j in 1:length(varNames)]

        # Check that there is at least 1 matching record in the HITRAN file
        length(cols[1]) == 0 ? throw(HitranEmptyError()) : nothing

        # Return a HitranTable struct populated with the obtained columns
        return HitranTable(mol=cols[1], iso=cols[2], νᵢ=cols[3], Sᵢ=cols[4], Aᵢ=cols[5], γ_air=cols[6], γ_self=cols[7], E″=cols[8], n_air=cols[9], δ_air=cols[10], global_upper_quanta=cols[11], global_lower_quanta=cols[12], local_upper_quanta=cols[13], local_lower_quanta=cols[14], ierr=cols[15], iref=cols[16], line_mixing_flag=cols[17], g′=cols[18], g″=cols[19])

    end

end

"""
    make_hitran_model(hitran::HitranTable, 
                      broadening::AbstractBroadeningFunction; 
                      wing_cutoff::Real=40, 
                      vmr::Real=0, 
                      CEF::AbstractComplexErrorFunction=HumlicekWeidemann32SDErrorFunction(), 
                      architecture = default_architecture)

Convenience function to make a HitranModel out of the parameters (Matches make_interpolation_model)

"""
function make_hitran_model(hitran::HitranTable, 
                           broadening::AbstractBroadeningFunction; 
                           wing_cutoff::Integer=40, 
                           vmr::Union{Real, Vector}=0, 
                           CEF::AbstractComplexErrorFunction=HumlicekWeidemann32SDErrorFunction(), 
                           architecture = default_architecture)

    if architecture isa GPU && !(CEF isa HumlicekWeidemann32SDErrorFunction)
        @warn "Cross-section calculations on GPU may or may not work with this CEF (use HumlicekWeidemann32SDErrorFunction if you encounter issues)"
    end

    return HitranModel(hitran=hitran, broadening=broadening , wing_cutoff=wing_cutoff , vmr=vmr, CEF=CEF, architecture=architecture)
end