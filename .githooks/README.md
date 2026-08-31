# Local privacy guard

The `pre-commit` hook rejects newly staged files under:

- `RRS_XCO2/resources/`
- `RRS_XCO2/manuscript/`
- `TODO/PRIVATE_*`

It also rejects common credential and private-key filenames. This complements
`.gitignore`: unlike an ignore rule, the hook still catches `git add -f`.

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
