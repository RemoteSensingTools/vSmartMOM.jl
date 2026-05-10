# Commit grouping draft — docs rewrite

Suggested 5-commit split (or 1 big commit, your call).

---

## 1. Pass 1 skeleton & nav

Restructure docs around a single RT/MOM narrative thread.

- New `docs/src/pages/concepts/` directory with 11 numbered pages (01 →
  08, plus 03a/03b/03c) following the spine: problem → MOM thesis →
  vector RTE & discretization → layer optical properties (with gas
  absorption + Mie/Rayleigh + mixing/δ-M as branches) → MOM solver
  (elemental → doubling → adding) → surfaces → linearization →
  architecture-agnostic kernels → inelastic (brief).
- `docs/make.jl` Concepts nav rebuilt to expose the arc; old per-module
  Overview pages dropped from nav (kept on disk as a safety backstop —
  to be deleted in a follow-up commit once migration is verified).

---

## 2. Pass 2 Concepts prose

Write each Concepts page in narrative order, with paper citations and
file:LINE anchors throughout.

Two load-bearing points threaded through the prose:
- **Exact finite-δ elemental, not the linear S2014 (19)–(20) limit.**
  The forward kernel uses Fell 1997 / SF2023-II Eqs (10)–(11) with
  `-expm1(-x)` and `expdiff_neg(a, b)` for numerical stability.
- **Scattering-only `N_doubl` sizing with total `τ_λ` for transmission.**
  SF2023-II Eqs (8)–(9). Lets line-by-line on hyperspectral grids
  run as one batched matmul per layer per Fourier moment on GPU.

Plus a "Why this is fast" section in Concepts/06 with the
`Ġ = -G·d(RR)·G` analysis: combined forward + linearized run costs
**less than 2× a forward-only run** because the dominant batched LU
is paid once and reused. ForwardDiff would pay `(1+N_params)×`,
finite differences `(1+N_state)×`.

---

## 3. Code docstrings

Paper-cited docstrings on the math-heavy kernels Concepts/04 and
Concepts/06 reference. Forward path:
`apply_D!`, `apply_D_SFI!`, `apply_D_matrix!`, `compute_geometric_progression!`,
`doubling_source_update!`, `doubling_rt_update!`, `get_elem_rt!`, all
four `interaction_helper!(::ScatteringInterface_*)` cases. Linearized path:
`elemental_lin!::elemental!`, `doubling_lin::doubling_helper!`,
`interaction_lin` module header, `lin_added_layer_all_params_helper!`,
`rt_kernel_lin::rt_kernel!`. Fixes several "Sanghavi & Stephens 2013"
miscitations (no such paper) → S2014 / Fell 1997 / SF2023-II.

The opaque `apply_D!` row-negation-then-sign-table sequence is now
fully explained with S2014 Eqs (29)–(32) and the
`R*_10 = D · R_10` rationale.

---

## 4. AGENTS.md + CLAUDE.md

Add `AGENTS.md` at repo root: 3-minute orientation for AI coding
agents. Includes the narrative spine, the two non-obvious tricks that
compound (exact finite-δ + constant-`N_doubl`), the 8-point
differentiator list with file:LINE anchors, the code map, conventions,
and a maintenance contract.

Update `CLAUDE.md` to point readers at AGENTS.md as the first read.

---

## 5. theory_references.md crib sheet

Overwrite `docs/dev_notes/theory_references.md` with the verified
equation↔code map covering all seven canonical papers in
`docs/papers/`. Sections A–J map every paper equation that the
Concepts arc quotes to a `file.jl:LINE` reference. Seven explicit
"design choice" passages identified for the prose to echo verbatim.
CITATION.bib entries drafted.

---

## Commit message style notes

- Imperative mood (per existing user preference).
- **No `Co-Authored-By:` lines** (per user feedback memory).
- The CLAUDE.md change in commit 4 is intentional even though it's a
  one-line edit — it pairs with AGENTS.md.
