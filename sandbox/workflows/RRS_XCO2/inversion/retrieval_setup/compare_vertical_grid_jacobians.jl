#!/usr/bin/env julia
# Compare serialized outputs of vertical_grid_jacobian_states.jl bit for bit.
using Serialization, LinearAlgebra, Printf
# usage: julia compare_vertical_grid_jacobians.jl BASELINE.jls CANDIDATE.jls...
base = deserialize(ARGS[1])
for name in ARGS[2:end]
    new = deserialize(name)
    println("== $name vs baseline (7acab570)")
    for key in sort(collect(keys(base)))
        b, n = base[key], new[key]
        meas = b.measurement == n.measurement
        jac = b.jacobian == n.jacobian
        jac_rest = b.jacobian[:, 2:end] == n.jacobian[:, 2:end]
        k1 = b.jacobian[:, 1]; k1n = n.jacobian[:, 1]
        @printf("  %-20s measurement_identical=%s jacobian_identical=%s cols2:30_identical=%s  |K1_new-K1_base|/|K1_base|=%.4g\n",
                key, meas, jac, jac_rest, norm(k1n - k1) / norm(k1))
        if hasproperty(b, :xco2)
            @printf("  %-20s xco2 base=%.9f new=%.9f\n", "", b.xco2, n.xco2)
        end
    end
end
