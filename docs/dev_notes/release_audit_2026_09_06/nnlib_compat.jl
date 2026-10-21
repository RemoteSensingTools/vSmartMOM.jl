using NNlib
A=ones(2,2,3); B=ones(2,2,1)
println("BEFORE_VSMARTMOM batched_mul shape=",size(NNlib.batched_mul(A,B)))
using vSmartMOM
try
    println("AFTER_VSMARTMOM batched_mul shape=",size(NNlib.batched_mul(A,B)))
catch err
    println("AFTER_VSMARTMOM ERROR=",sprint(showerror,err,catch_backtrace()))
end
println("DISPATCH=",which(NNlib.batched_mul,(typeof(A),typeof(B))))
