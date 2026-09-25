"""Rerun one observed target in all three methods against frozen results."""
import json
import numpy as np
from .paths import DATA, REFERENCE
from .scenes import load_scene
from .model import Physics
from .run_accuracy import chain
from .evaluation import rescore_record

def main():
    scene=2102
    p=Physics(load_scene(DATA/"scenes"/f"{scene}.json"))
    methods=["Fixed20","Hold","Deadline"]
    plans={m:json.loads((REFERENCE/"plans"/f"{scene}-{m}.json").read_text())["plan"] for m in methods}
    candidates=set.intersection(*[{v["target"] for v in plan} for plan in plans.values()])
    target=min(candidates)
    checks=[]
    for method in methods:
        actual=rescore_record(chain(p,target,[v for v in plans[method] if v["target"]==target],scene,200))
        expected=json.loads((REFERENCE/"chains"/f"{scene}-{method}-{target}.json").read_text())
        for key in ["time_position_rmse_m","time_velocity_rmse_m_s","final_position_rmse_m","final_velocity_rmse_m_s",
                    "max_position_error_per_realization_m","max_velocity_error_per_realization_m_s","curves_1s"]:
            np.testing.assert_allclose(actual[key],expected[key],rtol=1e-8,atol=1e-8,err_msg=key)
        assert actual["threshold_pass_pct"]==expected["threshold_pass_pct"]
        checks.append({"method":method,"target":target,"position_rmse_m":actual["time_position_rmse_m"],
                       "velocity_rmse_m_s":actual["time_velocity_rmse_m_s"],"truth_hash_matches":actual["truth_sha256"]==expected["truth_sha256"]})
        print(method,"PASS",flush=True)
    report={"status":"PASS","scene":scene,"realizations_per_method":200,"checks":checks,
            "scope":"three full nonlinear chains; other chains verified by reaggregation and plan replay"}
    (DATA.parent/"docs/smoke-results.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
