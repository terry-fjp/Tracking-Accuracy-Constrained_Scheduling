"""Reaggregate physical-error endpoints and audit every saved observation plan."""
import json,hashlib,ast,os
from pathlib import Path
import numpy as np
from .evaluation import THRESHOLDS, BASELINE, key, rescore_record
from .model import Physics
from .scheduler import Visit,replay
from .paths import ROOT, DATA, REFERENCE
from .scenes import load_scene
OUT=Path(os.environ.get('AEROSPACE_ACCURACY_OUTPUT',REFERENCE))

def main():
    protocol=json.loads((OUT/'protocol.json').read_text())
    summaries=json.loads((OUT/'scene-results.json').read_text());checks=[];chains=0;trajectories=0
    for scene in range(2101,2109):
        model_scene=load_scene(DATA/'scenes'/f'{scene}.json')
        model_scene['config'].update(protocol.get('config_overrides',{}))
        scene_hash=hashlib.sha256(json.dumps(model_scene,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        p=Physics(model_scene);truth={}
        for method in ['Fixed20','Hold','Deadline']:
            record=next(x for x in summaries if x['scene']==scene and x['method']==method)
            plan=json.loads((OUT/'plans'/f'{scene}-{method}.json').read_text())
            assert plan['scene_sha256']==scene_hash and plan['hold_s']==200.0
            assert plan['search']['evaluations']==48 and plan['search']['seed']==3101
            plan_hash=hashlib.sha256(json.dumps(plan['plan'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
            audit=replay(p,[Visit(**v) for v in plan['plan']])
            assert np.allclose(audit['objective'],plan['metrics']['objective'],rtol=1e-9,atol=1e-8)
            assert audit['max_label_error']<1e-8
            assert audit['resource_violations']==0
            if method!='Fixed20':assert audit['hold_violations']==0
            assert abs(audit['laser_s']-record['laser_s'])<1e-9
            rows=[]
            for i in range(p.N):
                r=json.loads((OUT/'chains'/f'{scene}-{method}-{i}.json').read_text());rows.append(r)
                assert r['scene_sha256']==scene_hash and r['plan_sha256']==plan_hash
                assert r['hold_s']==200.0 and r['repetitions']==200
                pr=np.array(r['max_position_error_per_realization_m']);pv=np.array(r['max_velocity_error_per_realization_m_s'])
                assert len(pr)==len(pv)==200 and np.isfinite(pr).all() and np.isfinite(pv).all()
                assert int(((pr<=BASELINE[0])&(pv<=BASELINE[1])).sum())==r['full_interval_pass_count']
                if i in truth:assert truth[i]==r['truth_sha256']
                truth[i]=r['truth_sha256'];chains+=1;trajectories+=200
                rates=[]
                for pair in [key(pair) for pair in THRESHOLDS]:
                    a,b=map(float,pair.split('/'));rate=100*float(np.mean((pr<=a)&(pv<=b)))
                    assert abs(rate-r['threshold_pass_pct'][pair])<1e-9;rates.append(rate)
                assert all(a<=b for a,b in zip(rates,rates[1:]))
            for field,target in [('time_position_rmse_m','position_rmse_m'),('time_velocity_rmse_m_s','velocity_rmse_m_s')]:
                actual=float(np.sqrt(np.mean([r[field]**2 for r in rows])))
                assert abs(actual-record[target])<1e-8
            actual_pass=100*sum(r['full_interval_pass_count'] for r in rows)/(p.N*200)
            assert abs(actual_pass-record['actual_full_interval_pass_pct'])<1e-9
            for pair in rows[0]['threshold_pass_pct']:
                assert abs(np.mean([r['threshold_pass_pct'][pair] for r in rows])-record['threshold_pass_pct'][pair])<1e-9
            checks.append({'scene':scene,'method':method,'resource_violations':audit['resource_violations'],
                           'hold_violations':audit['hold_violations'],'metrics_reaggregated':True,'paired_truth_identical':True})
        del p
    totals=json.loads((OUT/'summary.json').read_text())
    for total in totals:
        rows=[r for r in summaries if r['method']==total['method']]
        for field in ['position_rmse_m','velocity_rmse_m_s','final_position_rmse_m','final_velocity_rmse_m_s','actual_full_interval_pass_pct','laser_s']:
            assert abs(np.mean([r[field] for r in rows])-total[field])<1e-8
        for pair in total['threshold_pass_pct']:
            assert abs(np.mean([r['threshold_pass_pct'][pair] for r in rows])-total['threshold_pass_pct'][pair])<1e-9
    report={'status':'PASS','plans_checked':len(checks),'target_chains':chains,'target_method_realizations':trajectories,
            'paired_unique_truth_realizations':trajectories//3,'targets':chains//3,
            'thresholds_per_chain':len(THRESHOLDS),'checks':checks,
            'hold_s':200.0,'input_and_plan_hashes':'PASS',
            'scope':'reaggregation and plan replay; not external physical validation'}
    (OUT/'final-audit.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))

if __name__=='__main__':main()
