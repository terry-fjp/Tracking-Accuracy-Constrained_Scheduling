# Tracking-Accuracy-Constrained Scheduling

Research code and synthetic scene inputs for **Tracking-Accuracy-Constrained Scheduling for Active–Passive Observation Constellations of Noncooperative Space-Object Clusters**.

**Authors:** Jiapeng Feng, Yuqing Li, Sensen Guo, Wei Zhao, Xuhai Wang, and Jun Zhou.  
**Corresponding author:** Jiapeng Feng — fengjp@microsate.com  
**Repository:** https://github.com/terry-fjp/Tracking-Accuracy-Constrained_Scheduling

This repository contains the scheduling and estimation implementation, experiment configurations, eight synthetic orbital scenes, and selected diagnostic data. The manuscript, LaTeX template, publication figures, plotting scripts, and full experiment-result archives are not included in this directory. Experiment outputs are generated locally by the commands below.

## Research scope

The implementation schedules joint passive optical angle measurements and active laser ranging for noncooperative space-object clusters. It combines:

- Tracking-accuracy constraints based on position–velocity covariance and a 200 s measurement-free hold interval.
- Adaptive observation durations and an optional initialization-priority rule.
- Time-dependent common-field acquisition and serial laser service.
- Independent target posteriors and covariance propagation between observations.
- Hierarchical self-adaptive differential evolution (H-SADE) for scheduling search.
- Nonlinear truth propagation and extended Kalman filter (EKF) evaluation using paired random realizations across methods.

Planning covariance feasibility and actual estimation-error compliance are evaluated separately. The experiments investigate conditional performance improvements; they do not establish universal superiority or flight readiness.

## Repository structure

```text
Tracking-Accuracy-Constrained_Scheduling/
├── README.md
├── config/
│   ├── evaluation.json           # Actual-error acceptance thresholds
│   └── experiments.json          # Planning and measurement-noise conditions
├── data/
│   ├── scene.schema.json         # Scene data schema
│   ├── scenes/
│   │   ├── 2101.json ... 2108.json
│   │   └── manifest.json         # Scene counts and numerical-content hashes
│   └── diagnostics/
│       ├── mechanism.json
│       ├── mechanism-curves.json
│       └── scene-description.json
└── src/
    └── space_cluster_scheduling/
        ├── __init__.py
        ├── __main__.py           # Main command-line interface
        ├── paths.py              # Repository-relative data/output locations
        ├── scenes.py             # Scene loading and validation
        ├── model.py              # Orbital geometry and covariance model
        ├── scheduler.py          # Scheduling, optimization, and plan replay
        ├── nonlinear_math.py     # Nonlinear dynamics and measurements
        ├── evaluation.py         # Actual-error acceptance calculations
        ├── run_accuracy.py       # Three-method scheduling and EKF experiments
        ├── run_sweep.py          # Multiple experiment conditions
        ├── run_controls.py       # Matched controls and search-seed experiments
        ├── validate.py           # Scene integrity checks
        ├── check_formulas.py     # Local numerical consistency checks
        ├── audit_accuracy.py     # Saved-plan and result audits
        ├── rescore.py            # Reclassify saved estimation errors
        └── smoke.py              # Replay against existing P1 results
```

## Environment setup

Use Python 3.12 and NumPy. The commands below target Linux/Bash; `run_controls.py` reads `/proc/cpuinfo` to record processor information. NumPy 2.3.5 is recommended for consistency with the manuscript preparation environment. No package installer configuration is included, so expose the source directory through `PYTHONPATH`.

Run these commands from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install numpy==2.3.5
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

Limiting numerical-library threads avoids CPU oversubscription when using multiple experiment workers. Choose the number of workers according to available memory and CPU resources.

## Input and numerical checks

Display available commands and validate the supplied scenes:

```bash
python -m space_cluster_scheduling --help
python -m space_cluster_scheduling validate-scenes
```

Scene validation checks the supported data fields, target and observer counts, and canonical numerical-content hashes against the manifest. The supplied scenes contain 138 targets in total.

Run local formula checks:

```bash
python -m space_cluster_scheduling check-formulas
```

This writes `results/formula-checks.json` and checks selected dynamics and measurement Jacobians, covariance updates, and hold propagation. These are numerical consistency checks, not external physical validation.

## Main experiments

### Run the reference condition

```bash
python -m space_cluster_scheduling simulate --condition P1 --workers 1
```

This runs all eight scenes with three methods, search seed 3101, 48 complete-plan evaluations per method, and 200 Monte Carlo realizations per target. It is a full experiment, not a short smoke test.

The output method identifiers map to the manuscript labels as follows:

| Output identifier | Manuscript label | Observation duration | Initialization priority | Hold screening |
|---|---|---|---|---|
| `Fixed20` | F20 | Fixed 20 s | No | No hard rejection |
| `Hold` | AD | Adaptive | No | Yes |
| `Deadline` | IP-AD | Adaptive | Yes | Yes |

By default, outputs are written to `results/experiments/P1/`. A different directory can be selected with `--output`:

```bash
python -m space_cluster_scheduling simulate --condition P1 --workers 2 --output results/custom/P1
```

### Audit and rescore completed results

After the reference simulation has completed:

```bash
python -m space_cluster_scheduling audit --condition P1
python -m space_cluster_scheduling rescore --condition P1
```

For a custom result directory, pass the same `--output` path to both commands. The audit replays plans, verifies resource constraints and applicable hold constraints, checks saved hashes and paired truth, and reaggregates accuracy metrics.

Rescoring updates saved acceptance rates and produces `threshold-summary.csv`. It does not optimize new schedules or rerun the EKF. Changing planning thresholds requires new simulations; changing actual-error acceptance thresholds can be handled through rescoring.

### Run multiple conditions

Run a selected subset:

```bash
python -m space_cluster_scheduling.run_sweep --conditions P1 P2 M1 M1P4 --workers 2
```

Run every configured condition:

```bash
python -m space_cluster_scheduling.run_sweep --workers 2
```

The configuration contains 20 conditions: four planning levels at the reference measurement noise, plus a full four-by-four measurement-noise/planning-level grid. The sweep automatically audits completed conditions.

## Matched controls and repeated search seeds

Run the control experiments:

```bash
python -m space_cluster_scheduling.run_controls --workers 2
```

The protocol includes:

- Six matched methods: F20, IP-F20, F20-H, IP-F20-H, AD, and IP-AD.
- Five paired search seeds: 3101–3105.
- Eight scenes under the P1 and P4M1 conditions, with 48 plan evaluations per method.
- Fixed-duration scans of 10, 20, 30, and 40 s with initialization priority and hold screening, under P1 with seed 3101.
- Search budgets of 24, 48, and 96 evaluations on scenes 2101 and 2106 under P1.
- Recorded optimization wall time and process CPU time, separate from preprocessing and estimation evaluation.

`IP` enables initialization priority. `H` adds hold screening to fixed-duration scheduling without shortening or extending the prescribed visit. F20 and IP-F20 do not reject visits solely for failing the hold check; adaptive methods apply hold screening.

The controls label `P4M1` corresponds to `M1P4` in `config/experiments.json`: planning thresholds of 2000 m and 100 m/s, with measurement standard deviations of 50 microradians and 30 m.

For a smaller control run:

```bash
python -m space_cluster_scheduling.run_controls --smoke --workers 2
```

This still performs six optimizations for scene 2101 and 200 Monte Carlo realizations per target. Outputs go to the fixed directory `results/controls/`. Keep the same worker setting when continuing a control run: the saved protocol fingerprint includes the worker count and runtime metadata. For a changed protocol, use a separate checkout or preserve the existing results elsewhere before starting again.

## Experiment settings and interpretation

| Setting | Value |
|---|---|
| Synthetic scenes | 2101–2108 |
| Targets per scene | 12, 18, or 24 |
| Total targets across scenes | 138 |
| Planning horizon | 0–900 s |
| Initialization deadline | 120 s |
| Measurement and evaluation time step | 0.2 s |
| Measurement-free hold interval | 200 s |
| Actual-error evaluation interval | 120–900 s |
| Monte Carlo realizations | 200 per target and plan |
| Baseline actual-error acceptance | Position ≤ 500 m and velocity ≤ 25 m/s |

`config/evaluation.json` defines four actual-error threshold pairs: 500 m/25 m/s, 1000 m/50 m/s, 1500 m/75 m/s, and 2000 m/100 m/s. Full-interval compliance requires both errors to remain within their thresholds at every evaluation grid point. All targets contribute to the evaluation, including targets without scheduled observations.

Repeated search seeds measure scheduling-search variation. Monte Carlo realizations measure estimation variation for a fixed plan; these are different sources of variability. Repeated seeds and parameter settings do not constitute additional independent orbital scenes.

## Generated outputs and replay

Main experiments generate a directory for each condition containing:

```text
results/experiments/<condition>/
├── protocol.json
├── plans/
├── chains/
├── scene-summaries/
├── scene-results.json
└── summary.json
```

Auditing adds `final-audit.json`; rescoring adds `threshold-summary.csv`. Control experiments generate their own protocol, plans, cached target results, per-run records, and `all-runs.json` under `results/controls/`.

Existing results are reused subject to the checks implemented by each runner. Use a fresh output directory or separate checkout after changing code, inputs, or experimental settings. Full simulations may require substantial computation; no runtime guarantee is provided.

The separate command `python -m space_cluster_scheduling smoke` requires completed reference plans and chains in `results/experiments/P1/`, which are not bundled here. It also requires a `docs/` output directory. After generating P1 results, it can be run with:

```bash
mkdir -p docs
python -m space_cluster_scheduling smoke
```

This compares recomputed target-level estimation metrics against existing results and writes `docs/smoke-results.json`. The command uses the fixed P1 location; its input path is not redirected by `--output`.

## Modeling limitations

All orbital scenes are synthetic. The model assumes known observer ephemerides, target association, and availability of scheduled measurements. It does not model illumination, reflectance, laser-return link budgets, acquisition failure, or observer ephemeris errors. Truth and estimation use two-body dynamics with the specified process noise; perturbations such as J2 and atmospheric drag are not included.

Feasibility and compliance are checked on a discrete time grid. Covariance feasibility does not guarantee actual-error compliance, and numerical consistency checks do not establish continuous-time or flight performance.

## Funding and disclosure

This work was supported by the Science and Technology Commission of Shanghai Municipality, grant **25DZ3100700**, for the project “Research on Models for Onboard Intelligent Mission Decision-Making and Health Management” (English translation of the supplied project title).

OpenAI Codex, version **26.901.51231**, assisted with translation and code development. The authors reviewed the code, simulation results, and entire manuscript and take responsibility for the content.

The authors declare no conflicts of interest.

## Citation and license

The associated manuscript title and author list are provided above. Publication details and a DOI should be added when available.

No license file is included in this directory. No open-source license is asserted by this README; contact the corresponding author regarding reuse permissions.
