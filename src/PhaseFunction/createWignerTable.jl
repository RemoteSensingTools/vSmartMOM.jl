using WignerSymbols
using JLD2

Nmax = 250;
A_lnm = zeros(Float32,2Nmax+1,Nmax,2Nmax+1);
B_lnm = zeros(Float32,2Nmax+1,Nmax,2Nmax+1);

#arr_to_fill = Array{Any}(nothing, l+n+1, n+1);

function f(Nmax)
    cN = 0;
    cM = 0;
    for l = 0:2Nmax
        for n = 1:Nmax
            mᵒ    = max(l-n,n+1)
            m_max = min(l+n,Nmax)
            #println("l,n: ", l, " ", n, " " , m_max-mᵒ, " " , (n+l)-mᵒ)
            cN += max(0,m_max-mᵒ)
            cM += max(0,(n+l)-mᵒ)
            #= for m = mᵒ:n+l
                println()
                fac = sqrt(2(2l+1)*(2n+1)*(2m+1));
                #println(m," ",l, " ", n)
                @inbounds A_lnm[l+1,n,m-n+1] = fac*wigner3j(m,n,l,-1,1,0);
                if l>1
                    @inbounds B_lnm[l+1,n,m-n+1] = fac*wigner3j(m,n,l,-1,-1,2);
                end
            end
            println(l," ",n) =#
        end
    end
    return cN,cM
end

global 
global 


@save "wignerTable.jld2" A_lnm B_lnm



