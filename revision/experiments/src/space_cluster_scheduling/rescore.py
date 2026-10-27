"""Recompute acceptance tables from saved physical-error maxima."""
import json
import csv
from .paths import REFERENCE
from .evaluation import rescore_record, BASELINE, THRESHOLDS

def main(output=REFERENCE):
    def read(name): return json.loads((output/name).read_text())
    def write(name, data): (output/name).write_text(json.dumps(data, indent=2)+"\n")
    scenes=read("scene-results.json")
    for scene in scenes:
        records=[]
        for target in range(scene["target_count"]):
            name=f"chains/{scene['scene']}-{scene['method']}-{target}.json"
            record=rescore_record(read(name));write(name,record);records.append(record)
        scene["actual_full_interval_pass_pct"]=100*sum(r["full_interval_pass_count"] for r in records)/sum(r["repetitions"] for r in records)
        scene["threshold_pass_pct"]={k:sum(r["threshold_pass_pct"][k] for r in records)/len(records) for k in records[0]["threshold_pass_pct"]}
        scene["evaluation_baseline"]=list(BASELINE)
    write("scene-results.json",scenes)
    totals=read("summary.json")
    for total in totals:
        rows=[r for r in scenes if r["method"]==total["method"]]
        total["actual_full_interval_pass_pct"]=sum(r["actual_full_interval_pass_pct"] for r in rows)/len(rows)
        total["threshold_pass_pct"]={k:sum(r["threshold_pass_pct"][k] for r in rows)/len(rows) for k in rows[0]["threshold_pass_pct"]}
        total["evaluation_baseline"]=list(BASELINE)
    write("summary.json",totals)
    protocol=read("protocol.json")
    protocol.update(evaluation_threshold_pairs=[list(p) for p in THRESHOLDS], evaluation_baseline=list(BASELINE),
                    secondary_endpoint=f"all-grid actual position <={BASELINE[0]:g} m and velocity <={BASELINE[1]:g} m/s; all targets included",
                    threshold_sweep="Evaluator-only reclassification of this condition; plans and RMSE are unchanged by rescoring; planning settings are recorded in protocol.config_overrides")
    write("protocol.json",protocol)
    with (output/"threshold-summary.csv").open("w",newline="") as f:
        writer=csv.writer(f);writer.writerow(["method","position_threshold_m","velocity_threshold_m_s","full_interval_pass_pct"])
        for r in totals:
            for k,v in r["threshold_pass_pct"].items():writer.writerow([r["method"],*k.split("/"),v])
    print(json.dumps([{k:r[k] for k in ["method","threshold_pass_pct"]} for r in totals],indent=2))
