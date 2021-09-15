#=
 
This file contains helper functions to load and store partition-sum information relevant 
to absorption cross section calculation
 
=#

const tips_file_path = String(@__DIR__) * "/TIPS_2017.nc"
const TIPS_2017_ISOT_HASH_CONST = ncread(tips_file_path, "TIPS_2017_T")
const TIPS_2017_ISOQ_HASH_CONST = ncread(tips_file_path, "TIPS_2017_Q")

T_end_idxs = zeros(Int64, size(TIPS_2017_ISOT_HASH_CONST)[1], size(TIPS_2017_ISOT_HASH_CONST)[2])
Q_end_idxs = zeros(Int64, size(TIPS_2017_ISOQ_HASH_CONST)[1], size(TIPS_2017_ISOQ_HASH_CONST)[2])

for M in 1:size(TIPS_2017_ISOT_HASH_CONST)[1]
    for I in 1:size(TIPS_2017_ISOT_HASH_CONST)[2]
        t_end_idx = findfirst(x->x==-1, TIPS_2017_ISOT_HASH_CONST[M,I,:])
        q_end_idx = findfirst(x->x==-1, TIPS_2017_ISOQ_HASH_CONST[M,I,:])
        T_end_idxs[M,I] = isnothing(t_end_idx) ? length(TIPS_2017_ISOT_HASH_CONST[M,I,:]) : t_end_idx-1
        Q_end_idxs[M,I] = isnothing(q_end_idx) ? length(TIPS_2017_ISOQ_HASH_CONST[M,I,:]) : q_end_idx-1
    end
end

get_TT(M, I) = TIPS_2017_ISOT_HASH_CONST[M, I, 1:T_end_idxs[M,I]]
get_TQ(M, I) = TIPS_2017_ISOQ_HASH_CONST[M, I, 1:Q_end_idxs[M,I]]
