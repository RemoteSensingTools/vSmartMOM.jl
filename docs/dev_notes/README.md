# docs/dev_notes/ — Living Developer References

This directory holds **long-lived technical references** that agents and
contributors need beyond the public docs. Session handoffs, phase plans, and
commit drafts are deleted after the work lands; git history preserves them.

## Kept files

| File | What it is |
|---|---|
| `theory_references.md` | Verified paper-equation ↔ code-line crib sheet. Source of truth for all citations in docs/src. |
| `fourier_stream_resolution_plan.md` | Design specification for the v0.7 `nstreams` / Fourier-order API; referenced from Schema.md and release notes. |
| `LINEARIZATION_BUGS.md` | Catalog of 25 bugs found and fixed during the sanghavi → unified merge. Reference for understanding the linearized RT chain rule. |
| `OCO_RETRIEVAL_LINEARIZATION_HANDOFF.md` | Handoff and acceptance plan for a 33-state, three-band OCO-like retrieval with four aerosols, CO₂/H₂O profile scales, and bandwise Lambertian surfaces. |
| `raman_gpu_optimization.md` | GPU performance audit for the Raman (inelastic) path. Parked for S. Sanghavi; ~19,000 allocations/run identified. |
| `RAMAN_CODE_HANDOFF.md` | Orientation for the Raman code author: forward vs linearized paths, Raman dispatch map, what to check. |
| `MieGPUSpeedup.md` | GPU Mie implementation plan (DoubleSingle, Neumaier, five kernels). Prototype on CPU backend; full GPU pipeline unlanded. |
| `batched_kernel_benchmarks.md` | Benchmark results (A100, L40S) for custom KA kernels vs cuBLAS for batched matmul and inversion. |
| `ocean.md` | Coupled atmosphere-ocean RT design (Fresnel interface, quadrature grid mismatch, OceanOptics.jl spinoff). Phase 1 only; remaining phases deferred. |
| `vsmartmom_vs_vlidort_audit.md` | Technical comparison of vSmartMOM vs vLIDORT: feature gaps, community trust, and hardening roadmap. |
| `refactoring_roadmap.md` | Post-v2.1.0 cleanup items: kernel duplication, surface `@autolinearize`, and other deferred refactors. |
| `internal_api_cleanup.md` | Triage of the 48 internal-export symbols into four action buckets (de-export / promote / keep). Post-v2.0.0 follow-up. |
| `warning_inventory.md` | Record of the 48 missing-docstring warnings resolved in the warnonly retirement. Input set for `internal_api_cleanup.md`. |
