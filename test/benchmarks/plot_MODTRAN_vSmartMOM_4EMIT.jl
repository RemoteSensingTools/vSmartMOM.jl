using Plots
using DelimitedFiles

function fmt(x, n)
    s = string(round(x; digits=n))
    i = findfirst('.', s)
    i === nothing ? s * "." * "0"^n : rpad(s, i + n, '0')
end

H2O = collect(0.05:0.5:2.5)
AOT = collect(0.0001:0.1:0.7)
GNDALT = collect(0.0:5.0)
TSZ = [168.0, 178.0]
h2o = H2O[1]
aot = AOT[1]
gndalt = GNDALT[1]
tsz = TSZ[1]


MOD_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/MODTRAN_out/"
filename1(h2o, aot, gndalt, tsz) =
           "MODTRAN_H2O$(fmt(h2o,4))_AOT$(fmt(aot,4))_GNDALT$(fmt(gndalt,3))_TSZ$(fmt(tsz,1)).dat"
path1 = joinpath(MOD_DIR, filename1(h2o, aot, gndalt, tsz))
data1 = readdlm(path1, skipstart=2)

MOM_DIR = "/home/sanghavi/data/EMIT_MODTRANcomp/vSmartMOM_out/"
filename2(h2o, aot, gndalt, tsz) =
           "vSmartMOM_modtran_equiv_H2O$(fmt(h2o,4))_AOT$(fmt(aot,4))_GNDALT$(fmt(gndalt,3))_TSZ$(fmt(tsz,1)).dat"
path2 = joinpath(MOM_DIR, filename2(h2o, aot, gndalt, tsz))
# The vSmartMOM file has two '#'-prefixed header lines; skipping only one would
# leave the tokenised column-name line in data2 and splay sphalb into a phantom
# 9th column of empty strings.
data2 = readdlm(path2, comments=true, comment_char='#')