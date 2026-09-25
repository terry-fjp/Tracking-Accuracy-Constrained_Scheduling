"""One command for validation, rescoring, simulation and audits."""
import argparse
from pathlib import Path
from .paths import REFERENCE, ROOT
import json

def main():
    parser=argparse.ArgumentParser(prog="python -m space_cluster_scheduling")
    parser.add_argument("command", choices=["smoke","validate-scenes","rescore","audit","simulate","check-formulas"])
    parser.add_argument("--output", type=Path, help="Result directory (default: bundled reference results)")
    parser.add_argument("--workers", type=int, default=1, help="Independent scene processes for simulate (default: 1)")
    parser.add_argument("--condition",default="P1",help="Condition ID from config/experiments.json")
    args=parser.parse_args()
    if args.workers<1: parser.error("--workers must be positive")
    spec=json.loads((ROOT/"config/experiments.json").read_text())
    condition=next((c for c in spec["conditions"] if c["id"]==args.condition),None)
    if condition is None:parser.error("Unknown condition ID")
    output=args.output.resolve() if args.output else ROOT/"results/experiments"/args.condition
    if args.command=="smoke":
        from .smoke import main as run
        run()
    elif args.command=="validate-scenes":
        from .validate import main as run
        run()
    elif args.command=="rescore":
        from .rescore import main as run
        run(output)
    elif args.command in ("audit","simulate"):
        from . import audit_accuracy, run_accuracy
        module=audit_accuracy if args.command=="audit" else run_accuracy
        module.OUT=output
        if args.command=="simulate": module.main(workers=args.workers,condition=condition)
        else: module.main()
    else:
        from .check_formulas import main as run
        run()

if __name__=="__main__":
    main()
