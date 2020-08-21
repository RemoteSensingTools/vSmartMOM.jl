using WignerSymbols
using JLD2

Nmax = 250;
A_lnm = zeros(Float32,2Nmax+1,Nmax,2Nmax+1);
B_lnm = zeros(Float32,2Nmax+1,Nmax,2Nmax+1);
for l = 0:2Nmax
    for n = 1:Nmax
        mᵒ = max(abs(l-n),n+1)
        for m = mᵒ:n+l
            fac = sqrt(2(2l+1)*(2n+1)*(2m+1));
            #println(m," ",l, " ", n)
            @inbounds A_lnm[l+1,n,m-n+1] = fac*wigner3j(m,n,l,-1,1,0);
            if l>1
                @inbounds B_lnm[l+1,n,m-n+1] = fac*wigner3j(m,n,l,-1,-1,2);
            end
        end
        println(l," ",n)
    end
end

@save "wignerTable.jld2" A_lnm B_lnm


