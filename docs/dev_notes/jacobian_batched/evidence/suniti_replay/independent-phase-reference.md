# Independent Float32 phase-derivative check

The five CPU entries failing the original three-aerosol tolerance were checked
against a Float64 mixture evaluated from the same supplied Float32 component
optics and component optical-depth derivatives. This avoids differences from
recomputing Mie optics at another precision. These entries differentiate aerosol
loading, which does not change an individual aerosol phase matrix.

For each aerosol i, evaluate S_i = tau_i * omega_i in Float64 and form
Z = (S_ray * Z_ray + sum(S_i * Z_i)) / (S_ray + sum(S_i)). The reference loading
derivative is dS_i * (Z_i - Z) / (S_ray + sum(S_i)).

| Phase / layer / spectral point / parameter | Local absolute error | Physical absolute error |
|---|---:|---:|
| Z++ / 3 / 1 / 16 | 1.7894e-8 | 1.2763e-6 |
| Z++ / 3 / 4 / 16 | 6.5545e-8 | 1.6306e-6 |
| Z0+ / 4 / 1 / 2 | 1.5488e-10 | 1.2383e-6 |
| Z0+ / 4 / 4 / 2 | 5.4276e-10 | 1.9698e-6 |
| Z0+ / 4 / 5 / 2 | 2.6275e-9 | 2.2154e-6 |

The local representation is closer in all five entries. This supports the
interpretation of the remaining mismatch as Float32 cancellation in successive
physical-column quotient derivatives. It is not evidence of full retrieval
convergence or validation of every derivative coordinate.

Reproduction: independent-phase-reference.jl, run from the optimization worktree test/
directory with Julia 1.12.6, JULIA_NUM_THREADS=4, OPENBLAS_NUM_THREADS=1.
Output: independent-phase-reference.log. The diagnostic disables the cross-representation
phase closeness assertion to print the failing entries, retains the exact
forward/zero-gas checks, and skips full RT propagation. Its 316 passing checks
must not be described as the full local-Jacobian suite.

All 154 stored source hashes in Sanghavi's study checkout matched on inspection.
