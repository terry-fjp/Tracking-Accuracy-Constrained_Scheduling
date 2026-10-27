"""Versioned public scene records and lossless numerical-model adapter."""
import json
import math
from pathlib import Path
ELEMENTS = ("radius_m", "inclination_deg", "raan_deg", "phase_deg")
def to_model(scene):
    if scene.get("schema_version") != "1.0":
        raise ValueError("Unsupported scene schema_version")
    from .model import Config
    config = scene["config"]
    Config(**config)  # Reject unknown planning parameters.
    for field in ("observers", "targets"):
        records = scene[field]
        if not records or [r["id"] for r in records] != list(range(len(records))):
            raise ValueError(f"{field}: IDs must equal ordered zero-based indices")
        for record in records:
            orbit = record["orbit"]
            if set(orbit) != set(ELEMENTS) or not all(math.isfinite(orbit[k]) for k in ELEMENTS):
                raise ValueError("Invalid orbital elements")
            if orbit["radius_m"] <= 6378137.0:
                raise ValueError("Orbit radius must exceed Earth radius")
    for t in scene["targets"]:
        for k in ("initial_sigma_position_m", "initial_sigma_velocity_m_s"):
            if not math.isfinite(t[k]) or t[k] <= 0:
                raise ValueError("Prior standard deviations must be positive")
    if config["dt"] <= 0 or config["init_deadline"] > config["horizon"]:
        raise ValueError("Invalid time configuration")
    return {"seed":scene["seed"], "config":config,
            "observers":[[r["orbit"][k] for k in ELEMENTS] for r in scene["observers"]],
            "targets":[[r["orbit"][k] for k in ELEMENTS] for r in scene["targets"]],
            "initial_sigma_r":[r["initial_sigma_position_m"] for r in scene["targets"]],
            "initial_sigma_v":[r["initial_sigma_velocity_m_s"] for r in scene["targets"]],
            **{k:scene[k] for k in ("epoch_utc","frame","target_source")}}
def load_scene(path):
    return to_model(json.loads(Path(path).read_text()))
