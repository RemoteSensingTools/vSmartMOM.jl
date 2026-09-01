# Local privacy guard

The `pre-commit` hook rejects newly staged files under these paths at the
repository root or beneath a relocated workflow directory:

- `RRS_XCO2/resources/`
- `RRS_XCO2/manuscript/`
- RRS/XCO2 generated `*.png`, `*.dat`, `*.nc`, and `*.csv` products
- `TODO/PRIVATE_*`
- `species_rayleigh_raman_todo.md`

It also rejects common credential and private-key filenames. This complements
`.gitignore`: unlike an ignore rule, the hook still catches `git add -f`.

The `pre-push` hook examines the complete history reachable from every pushed
ref and rejects any history containing the protected research paths or
generated RRS/XCO2 products. This also guards against accidentally publishing
local recovery refs with `git push --mirror`.

The checks match any path component named `RRS_XCO2`, not only a root-level
directory. Repository reorganizations therefore cannot bypass the privacy
boundary by moving the workflow beneath `sandbox/` or another prefix.

Enable the repository hooks once per clone/worktree with:

```bash
git config core.hooksPath .githooks
```

An intentional exception requires an explicit one-commit override:

```bash
VSMARTMOM_ALLOW_PRIVATE_PATHS=1 git commit
```

Do not use the override merely to silence the guard. Review every listed file
and confirm that it is intended for the public repository first.
