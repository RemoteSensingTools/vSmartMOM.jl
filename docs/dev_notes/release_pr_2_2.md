# Integrate Sanghavi's linearization and Raman work for v2.2

This release candidate consolidates the multisensor/surface work, retrieval-selected
analytic Jacobians, precision/convolution investigations, and Suniti Sanghavi's
updated OCO retrieval workflows. Scientific authorship is preserved in the commit
history and stated in the README, package metadata, and release documentation.

The integration also fixes Cox–Munk absolute radiance and tangent consistency,
NNlib coexistence, GPU test failure handling, and numerical dispatch ambiguities.
It removes misleading placeholder aerosol optics and adds pinned Taplo, scene
schema/parser checks, export/docstring inventories, and executable README checks.

Validation details and limits: [release validation](release_validation_2026-09-07.md).
The exact CPU/CUDA totals in that record must be final before opening this PR.
The hosted platform matrix should then run on the PR tip.

Merge with a merge commit to preserve the scientific-author records. Version
2.2.0 preserves the existing 2.1.0 tag; registration requires the package-author
version-sequence override or registry review. Verify documentation deployment
credentials before publishing. This PR does not register, tag, or deploy a release.

The separate GCHP sectional-AOD integration, archived CUDA graph experiment,
and private campaign datasets remain outside this candidate's release scope.
