using Plots
using DelimitedFiles

TSZ = [168.0, 178.0]
AOT = collect(0.0001:0.1:0.7)
H2O = collect(0.05:0.5:2.5)
GNDALT = collect(0.0:5.0)

const MOD_DIR = expanduser("~/data/EMIT_MODTRANcomp/MODTRAN_out")
const MOM_DIR = expanduser("~/data/EMIT_MODTRANcomp/vSmartMOM_out")

# Round to n decimals then right-pad with '0's.
function fmt(x, n)
    s = string(round(x; digits=n))
    i = findfirst('.', s)
    i === nothing ? s * "." * "0"^n : rpad(s, i + n, '0')
end

filename1(h2o, aot, gndalt, tsz) =
    "MODTRAN_H2O$(fmt(h2o,4))_AOT$(fmt(aot,4))_GNDALT$(fmt(gndalt,3))_TSZ$(fmt(tsz,1)).dat"
filename2(h2o, aot, gndalt, tsz) =
    "vSmartMOM_H2O$(fmt(h2o,4))_AOT$(fmt(aot,4))_GNDALT$(fmt(gndalt,3))_TSZ$(fmt(tsz,1)).dat"

for h2o in H2O, aot in AOT, gndalt in GNDALT, tsz in TSZ
    path1 = joinpath(MOD_DIR, filename1(h2o, aot, gndalt, tsz))
    path2 = joinpath(MOM_DIR, filename2(h2o, aot, gndalt, tsz))
    #@show filename1(h2o, aot, gndalt, tsz)
    #@show filename2(h2o, aot, gndalt, tsz)
    (isfile(path1) && isfile(path2)) || continue
    #isfile(path2) || continue
    data1 = readdlm(path1, skipstart=2)
    data2 = readdlm(path2, skipstart=2)
    # ...
end