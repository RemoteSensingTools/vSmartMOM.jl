# RAMI4ATM data files

The two large reference datasets used by the prototype `rami_*.jl`
scripts in this directory (`rami.jl`, `rami_v2.jl`, `rami_runAll.jl`,
`rami_runPolarized.jl`) are **no longer tracked in this repository**
to keep the clone size small. Fetch them locally before running the
scripts:

| File | Size | Source |
|------|------|--------|
| `RAMI4ATM_experiments_v1.0.json` | ~46 MB | <https://romc.jrc.ec.europa.eu/_www/php/RAMI4ATM/RAMI4ATM.php> — RAMI4ATM experiment manifest |
| `ref_kurucz_2006.nc`            | ~40 MB | Kurucz 2006 solar reference spectrum (published with the RAMI4ATM dataset; same project page) |

Place both into `sandbox/rami/` next to the prototype scripts and
they'll resolve via the existing relative-path lookups.

`createAbscoRAMI.jl` writes its own absorption-cross-section LUTs to
`/net/fluo/data2/data/Rami_CS_database/` (Caltech filesystem) and is
not affected by this change.

---

These files are listed in the repository `.gitignore` so they will
not be re-tracked accidentally on commit.
