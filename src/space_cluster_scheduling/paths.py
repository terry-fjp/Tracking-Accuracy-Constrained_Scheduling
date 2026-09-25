"""Paths resolved from the checkout, never the current working directory."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
REFERENCE = ROOT / "results/experiments/P1"
