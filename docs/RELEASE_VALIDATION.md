# Release validation — 2026-10-28

Compared with GitHub `main` at `55213c0008b14ecec868c080a7bc5bda5060ee6a` (32 files). Checks ran from a separate copy of the prepared release on Linux, Python 3.12.14, NumPy 2.3.5, with BLAS threads limited to one.

Passed checks:

- Syntax compilation of all 46 Python files; all 32 GitHub snapshot files checked against Git blob hashes.
- `python tools/prepare_revision.py --check-only` and actual preparation/extraction: frozen protocols, 16 source hashes, three validation-code hashes, frozen selection-data hash, and all 1339 baseline archive file hashes.
- CLI `--help`, `validate-scenes`, and `check-formulas`.
- `run_revision.py --regression`: exact default-plan reproduction for IP-AD and IP-F20-H at scene/seed pairs (2101, 3101) and (2106, 3105), using 48 complete-plan evaluations. Two targets per plan were replayed with 200 realizations each; truth hashes, position RMSEs, and full-interval pass counts matched the archive.
- `j2_revision.py --check`: finite-difference Jacobian, independent potential derivative, exact zero-J2 reduction, and step-halving propagation checks.
- `check_semantics.py`: extension-multiplier/cap/visibility behavior, ranking-weight interface, and retention of unperturbed calendars.
- `check_j2_replay.py`: exact equivalence of the dynamics adapter to the original two-body evaluator for its representative target.
- `check_figure3.py`: exact original numerical duration-diagnostic curves.

These are packaging, numerical-consistency, and representative regression checks. **The complete 960-schedule matrix, 160 J2 replays, and 120 new-seed validation runs were not rerun for this release.** Historical result summaries and audit records remain in `revision/experiments/reference-results/` and `revision/experiments/provenance/archived-checks/`; they are not installed as fresh-run PASS records. Plot/report scripts are included, but plot generation was not rerun in this packaging validation. This check does not establish physical or flight validation.
