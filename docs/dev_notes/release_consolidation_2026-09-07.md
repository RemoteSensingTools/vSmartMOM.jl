# Release consolidation — working record, 2026-09-07

Status: preparation in progress; this is not a release-readiness sign-off.
The isolated branch `integration/release-candidate` was created from `0540d37a`,
including the integration base `b524a0d8` and subsequent Jacobian optimization,
correctness, convergence, convolution, and precision investigations.

## Attribution contract

User-confirmed credit policy: Suniti Sanghavi receives scientific authorship
credit for the linearization and Raman work, including formulations, analytic
derivatives, and retrieval methodology. Software plumbing and speed work have
separate implementation credit. Apply this distinction to release summaries,
the manual, and the scientific references rather than deriving scientific
authorship from who ran a commit command.

Preserve the original author identity
and author date when integrating her commits. Prefer an ordinary merge when
the source history is appropriate; a selected cherry-pick preserves the author
while recording the integrator as committer. Do not squash mixed-author work
into a single cfranken-authored commit. Do not change repository-wide Git
identity or add co-author trailers (repository convention).

The user clarified: "It is her code, I am just the helper here." On that basis,
14 author records in the unpublished Jacobian work were corrected to
`sunitisanghavi <suniti.sanghavi@gmail.com>` on this local candidate. These cover
linearization formulation/correctness, equivalent-source and local optical
Jacobians, and scientific precision/convergence investigations. The other 14
commits in that segment retain cfranken's engineering/helper authorship. All
committer identities and original dates were retained. Earlier shared history,
the performance branch, and remote refs were not rewritten.

The [correction map](release_author_corrections_2026-09-07.json) records every
old/new ID in the 28-commit segment. Every commit's tree and message were checked
for exact equality; only author metadata and dependent parent IDs changed.
The corrected segment ends at `5534c7b1`. This is an authorship correction,
not a new scientific or numerical implementation.

Scientific credit is explicit in the public release notes. `Project.toml`
continues to list Suniti first, now with separate array entries for each
existing credited author. Published papers in `CITATION.bib` retain their
publication author order. Release publication by Suniti can be arranged at
handoff, but the publishing account does not replace commit authorship.

Examples requiring care: original `8bc60941` combines an equivalent-source Jacobian
implementation with an aerosol-normalization correction; `fe6b865e` combines
local optical Jacobian factorization with a performance change. Both now credit
Suniti as author with cfranken retained as committer. Do not classify an entire
commit from the word "Jacobian" in its subject.
Likewise, batching or fusing Raman kernels is performance engineering even
though the underlying Raman formulation is hers.

GitHub contribution graphs require a commit email associated with the author's
account and qualifying commits on the default branch. Association of the
existing email with her account has not been verified. See
[GitHub's contribution reference](https://docs.github.com/en/account-and-profile/reference/profile-contributions-reference).

## Snapshot before integration

The range `origin/main..0540d37a` contains 109 cfranken-authored and 14
sunitisanghavi-authored commits (including merges). The separate
`origin/main..origin/suniti_multi_sensor` range contains 25 commits, all
sunitisanghavi-authored. These ranges overlap; their counts must not be added
or interpreted as shares of scientific contribution.

Sanghavi-authored feature anchors already reachable from the candidate:

| Contribution | Commit |
| --- | --- |
| Height-aware observer radiances | `4c0758e7` |
| Height-resolved observer Jacobians | `0cae2db1` |
| Aerosol profile tangent consistency | `ee21b857` |
| Lambertian Legendre surface Jacobians | `372eaab9` |
| Surface pressure and profile gas Jacobians | `5ba80c51` |
| Wavelength-ordered opacity interpolation | `092ff2c1` |
| Moist-air opacity and CIA handling | `4774afd4` |
| External-solar TOA SFI | `41c788e3` |
| Retrieval-selected Jacobian plans | `55a0a73f` |
| Fourier convergence and OCO retrieval workflow | `b524a0d8` |

The adjacent `release_attribution_2026-09-07.json` records full commit IDs,
parents, authors, committers and dates for both ranges. Recheck it before the
final integration: unchanged source history should remain reachable, and
ported or author-corrected commits should have an explicit provenance map.

## Branch consolidation findings

- Local `feat/surface-split` (`51be2b13`) and remote
  `origin/feat/surface-split` (`46b33640`) have identical tracked trees under
  `src`, `ext`, `test`, `config`, `schemas`, and `.github`. Their apparent
  ahead/behind counts largely reflect rewritten history. Whole-tree differences
  concern privacy hooks and research products/notes. Use the sanitized remote
  ancestry already present in the candidate, rather than merging the old
  history back in.
- The original checkout has an unresolved `.gitignore` conflict and staged
  privacy hooks; it is left intact.
- Later `origin/suniti_multi_sensor` work includes SIF normalization/restart
  corrections and retrieval campaign tooling. Review and port relevant public
  source/test changes with their authorship intact; workflow relocation makes
  a blind whole-branch merge inappropriate. This integration is outstanding.
- Other feature branches still need an ancestry/patch-equivalence inventory.

## Outstanding release gates

Revalidate the earlier [release audit](release_readiness_2026-09-06.md) against
this candidate. Its aerosol-reference normalization, strict-docs attachment,
and angular-cache dimension defects have fixes on the Jacobian branch.
Cox–Munk source scaling and forward/Jacobian consistency, NNlib coexistence,
and the shipped GPU runner still require work. The candidate is not yet ready
to register.

Remaining work includes complete public-API/docstring inventory, executable
documentation examples, schema/parser agreement, pinned Taplo lint/format
checks, version/migration reconciliation, release automation, fresh dependency
resolution, and CPU/CUDA plus available platform validation on the final
candidate. The existing `v2.1.0` tag must not be reused. No release version has
been selected and no tag, push, registration, or publication has been made.
