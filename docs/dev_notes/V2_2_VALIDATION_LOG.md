# vSmartMOM v2.2 consolidation — validation log

One entry per commit on `release/v2.2-consolidation`. A commit cannot carry
its own SHA, so each entry names the change; the SHA appears in `git log`.
Environment unless stated: host curry.gps.caltech.edu, Julia 1.12.5,
AtmosphericAbsorption 0.1.2, CPU (`CUDA_VISIBLE_DEVICES=''`),
`JULIA_NUM_THREADS=4`, `OPENBLAS_NUM_THREADS=1`, tests run from `test/`.

Review gate for every commit: focused tests → independent Claude reviewer →
Codex (`gpt-6.1-sol`, reasoning effort xhigh, read-only) → fixes → re-review
by both.

## 1. docs: decision log, ExoOptics audit, handoff copy

Docs only (`docs/dev_notes/`). Handoff copied verbatim from Christian
Frankenberg's checkout (SHA-256 identical to the source file).
- Claude reviewer: factual errors in §1.1, §1.3, §1.5, §2.2 → fixed; re-review
  confirmed, final wording nits applied as suggested.
- Codex: APPROVE-WITH-FIXES (four documentation findings) → fixed; re-review
  APPROVE-WITH-FIXES (two wording nits) → applied as suggested.
