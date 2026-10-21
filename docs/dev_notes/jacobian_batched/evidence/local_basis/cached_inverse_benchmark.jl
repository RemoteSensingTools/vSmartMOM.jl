# Reproduce the rejected LU-cache experiment in an isolated checkout after
# applying cached_inverse_experiment.patch. The production branch omits it.
# Run from test/ with AUDIT_BACKEND=cuda AUDIT_NSPEC=512 AUDIT_STREAMS=16.
using vSmartMOM
isdefined(vSmartMOM.CoreRT,:_CACHED_JACOBIAN_INVERSE_ENABLED) ||
    error("Apply cached_inverse_experiment.patch in an isolated checkout first")
for key in ("AUDIT_COMPARE_VENDOR","AUDIT_COMPARE_MEDIUM","AUDIT_COMPARE_TILES")
    ENV[key]="false"
end
output=get(ENV,"AUDIT_OUTPUT","/tmp/cached-inverse")
ENV["AUDIT_OUTPUT"]=output*"-reference.toml"
include("../../local_basis_benchmark.jl")
vSmartMOM.CoreRT._CACHED_JACOBIAN_INVERSE_ENABLED[]=true
ENV["AUDIT_OUTPUT"]=output*"-cached.toml"
main()
