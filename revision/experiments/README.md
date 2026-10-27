# Frozen October revision experiments

This directory contains a separate, frozen copy of the numerical source used for the October revision. The main protocol fingerprint is `897fc155a0a40c679acec45d26910c621b77295215516c8d012c1cbfc0a2d975`. The warning-variant validation fingerprint is `a7c65e81362b54812bc5411f24e5c875691ab46993f07bc152a678358e85590d`.

Run the commands below from the **repository root**, using Python 3.12 and the root `requirements-plot.txt`. The original `requirements.txt` in this directory is retained as an environment record; the root plotting requirements also pin the transitive `six` dependency. No local `vendor/` runtime is required.

## Prepare and check

```bash
python tools/prepare_revision.py --check-only
python tools/prepare_revision.py
python revision/experiments/run_revision.py --regression
python revision/experiments/j2_revision.py --check
python revision/experiments/check_semantics.py
python revision/experiments/check_j2_replay.py
python revision/experiments/check_figure3.py
```

The first command verifies protocol fingerprints, all 16 frozen source hashes, the three validation-code hashes, the selection summary, and all 1339 baseline archive entries. The second extracts `reference-archive.zip` and creates output directories. Existing differing baseline files or selection summaries cause an error; the helper does not overwrite them. `--check-only` does not run scientific experiments or create outputs.

The regression command performs four full 48-evaluation searches and representative 200-realization target replays, comparing with archived plans and metrics. It is a meaningful regression, not an instant smoke test. The remaining checks cover J2 derivatives/propagation, duration/weight interfaces and unperturbed execution masks, equivalence of the two-body replay adapter, and the original numerical duration diagnostic. Run them sequentially.

## Full reproduction: command order matters

```bash
# 960 unique schedules; some reference evaluations are imported after regression.
python revision/experiments/run_revision.py --workers 4
# Reproduce all 80 default plans and remeasure planning times.
python revision/experiments/finalize_checks.py --baseline
# 160 fixed-calendar perturbed-truth replays.
python revision/experiments/j2_revision.py --workers 4
python revision/experiments/finalize_checks.py
python revision/experiments/audit_population_admission.py

# 120 runs on search seeds 3106–3110; frozen warning-0.7 candidate and two controls.
python revision/experiments/parameter-update/run_validation.py
python revision/experiments/parameter-update/check_and_summarize.py

# Generate main summary, plots, and TeX inserts LAST.
python revision/experiments/summarize_revision.py
```

The validation checker deliberately requires `results/summary.json` to have the exact hash of the historical data used to select warning 0.7. Preparation initializes this file from immutable `reference-results/summary.json`; it is a **selection record**, not a summary of a new run. Run the validation checker **before** the final main summary command, which regenerates timing-dependent fields and changes that file's hash. To audit validation again after regenerating the main summary, use a fresh work copy with the frozen selection record and copies of the generated validation outputs. Do not remove the hash check or silently replace a newly generated summary.

`resume_matrix.py` can resume the main matrix after an interruption, skipping already saved jobs before geometry setup. It uses four workers and the same frozen runner. Check inputs with `tools/prepare_revision.py --check-only` first. Do not run separate groups concurrently against the same output directory. The full protocol needs substantial computation and storage; no portable runtime guarantee is available. Planning CPU/wall times are recorded per run and vary by hardware. Monte Carlo evaluation times include valid cache reuse and are not standalone estimator benchmarks.

## Included inputs and historical evidence

- `src/`, `config/`, `data/`: frozen numerical source, original eight scenes, diagnostics, and 32 nested-population inputs (families 5101–5108, sizes 24/48/72/96).
- `protocol.json`: hold/heuristic/population/J2 settings, five paired search seeds 3101–3105, budget 48, and 200 realizations per target.
- `parameter-update/`: frozen candidate, validation protocol, runner, and audit/report generator. Warning 0.8 remains the original reference; warning 0.7 is a separate candidate.
- `reference-archive.zip`: 80 exact baseline plans, 80 run records, and 1179 referenced chains. Their 1339 SHA-256 values are in `provenance/reference-archive-hashes.json`. The runner uses the extracted `reference-archive/` by default; `REVISION_REFERENCE_ARCHIVE` can point to an identical expanded copy.
- `reference-results/`: compact **historical** summaries, paired differences, seed means, and admission audit. These are supplied for comparison; the large full raw output caches are excluded.
- `provenance/archived-checks/`: historical PASS records, separated from active `checks/`. They are never installed as current-run PASS records. `provenance/runtime.json` describes the prior execution environment with local machine paths removed.

One-time preparation scripts that regenerated the frozen protocol and local wait/continuation helpers are excluded. The plotting runner and its required mechanism builder are included. Root source and frozen source both expose the same parameter interfaces; future root edits do not automatically change the frozen experiment copy.

## Outputs

Fresh runs write `results/plans`, `results/runs`, `results/chains`, `results/j2`, and `parameter-update/results`. Successful numerical and completeness checks write active `checks/` and `parameter-update/checks.json`. Full aggregation requires all 960 schedules, 160 J2 replays, and all 80 default-plan reproductions; partial runs cannot produce a complete summary. The validation checker separately requires all 120 runs.

Plotting/report scripts write generated TeX inserts and figures into `revision/paper/`, and an analysis report into `revision/audit/`. These generated fragments require the manuscript's macros to compile; they are not a standalone manuscript or a separate supplementary-material submission. No paper source or editor correspondence is included in the repository.

## Interpretation

Each scene RMSE pools squared errors across targets, realizations, and evaluation times. Eight scene RMSEs are averaged equally for each search seed; mean and sample standard deviation use the five seed means. Compliance counts target–realization pairs meeting both error limits at every evaluation grid point, retaining all targets, including unobserved ones. Search-seed standard deviations are not confidence intervals, and 200 realizations are not 200 independent orbital scenes.

The main matrix contains 320 hold schedules, 360 heuristic schedules (40 overlap hold-200/IP-AD), and 320 population schedules, hence 960 unique schedules. J2 replays keep the original timestamps, cancel whole infeasible visits without replacement, and use the same executed calendar for both estimator modes. They examine model mismatch and matched J2 estimation; they do not create a J2-aware optimal schedule or model drag, observer-state uncertainty, or realized laser returns. Copied two-body `hold_risk` labels are not new J2 hold guarantees.

The 96-target study is an overload boundary with all targets retained in the denominator, not a proven constellation capacity. Validation seeds 3106–3110 reuse known geometries and Monte Carlo inputs; they test optimization-seed variation rather than generalization to unseen scenes. Covariance feasibility never implies guaranteed actual-error compliance.
