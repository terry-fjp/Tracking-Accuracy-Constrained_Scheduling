# Changes prepared on 2026-10-28

Compared with GitHub main at `55213c0008b14ecec868c080a7bc5bda5060ee6a`:

- Expose duration-extension multipliers and ranking weights through the existing configuration. Defaults remain `(1.0, 1.5)` and `(0.5, 0.3, 0.2)`. Reject invalid configurations. The numerical evaluator is unchanged.
- Add the frozen October revision experiments: hold/heuristic sensitivity, nested target populations, fixed-calendar J2 replay, and the warning-0.7 variant evaluated on new optimization seeds.
- Preserve the original warning-0.8 reference. The warning-0.7 candidate is a separate named variant and does not replace defaults.
- Include exact frozen source/protocol snapshots, synthetic population inputs, compact historical summaries, archived audit records, and the 80-run baseline archive needed by the revision runner.
- Add pinned dependencies, repository hygiene files, a standard-library preparation/integrity helper, and runnable documentation.

Historical summaries and archived checks describe prior completed experiments. They are not evidence of a fresh run of this release. Full raw revision output caches, local dependencies, manuscript sources/proofs, and editor correspondence are not included. No original GitHub file is deleted.
