"""Physical-error acceptance thresholds, separate from planning covariance."""
import json
import numpy as np
from .paths import ROOT
CONFIG = json.loads((ROOT / "config/evaluation.json").read_text())
THRESHOLDS = tuple(tuple(pair) for pair in CONFIG["position_velocity_thresholds"])
BASELINE = tuple(CONFIG["baseline_pair"])
def key(pair):
    return f"{pair[0]:g}/{pair[1]:g}"
def score_errors(position, velocity):
    """Joint inclusive acceptance on per-realization full-interval maxima."""
    p, v = np.asarray(position), np.asarray(velocity)
    if p.ndim != 1 or p.shape != v.shape or not len(p):
        raise ValueError("Error arrays must be nonempty paired vectors")
    if not np.isfinite(p).all() or not np.isfinite(v).all() or (p < 0).any() or (v < 0).any():
        raise ValueError("Errors must be finite and nonnegative")
    counts = {key(pair): int(np.count_nonzero((p <= pair[0]) & (v <= pair[1]))) for pair in THRESHOLDS}
    return {k: 100 * n / len(p) for k, n in counts.items()}, counts[key(BASELINE)]
def rescore_record(record):
    rates, count = score_errors(record["max_position_error_per_realization_m"], record["max_velocity_error_per_realization_m_s"])
    record.update(threshold_pass_pct=rates, full_interval_pass_count=count,
                  full_interval_pass_fraction=count / record["repetitions"],
                  evaluation_baseline=list(BASELINE))
    return record
