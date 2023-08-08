function energy_level(mol_idx, v, J)
    E = 0.0
    for k = 0:4
        tmpE1 = pow(v+0.5, k)
        for l =0:4
            tmpE2 = tmpE1 * Y[k,l] * pow(J*(J+1),l)
            E = E + tmpE2
        end
    end
    return E; #in cm^{-1}
end