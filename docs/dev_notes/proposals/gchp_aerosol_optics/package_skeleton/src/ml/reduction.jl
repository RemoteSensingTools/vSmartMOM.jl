# Dimensionality reduction for ML inputs (brief §8b). Per-layer profiles over 50–72 levels
# are too high-dimensional to feed raw. Compress each profile field onto a shared EOF/SVD
# basis and store only the leading scores.

"""
    EOFBasis

A frozen reduction basis for one profile field (e.g. extinction profile, effective-radius
profile, per-species-fraction profile, RH/AW profile): the field mean, the leading
eigenvectors (columns), and the explained-variance for choosing rank `k`. Persist this with
the dataset — inputs are only meaningful relative to the basis that produced them.
"""
struct EOFBasis
    # TODO: mean::Vector, components::Matrix (nlev × k), explained_variance::Vector
end

"""
    fit_eof(profiles::AbstractMatrix; k=nothing, var_threshold=0.99) -> EOFBasis

`profiles` is (nlev × nscenes). Center, `svd`, keep `k` components (or enough for
`var_threshold` of variance). Standard EOF/PCA.
"""
function fit_eof end   # TODO (new)

"""
    eof_reduce(profile::AbstractVector, basis::EOFBasis) -> Vector   # forward: profile → scores
    eof_expand(scores::AbstractVector, basis::EOFBasis) -> Vector    # inverse: scores → profile
"""
function eof_reduce end   # TODO (new)
function eof_expand end   # TODO (new)
