"""Prospectively fixed new-search-seed assessment of a selected warning policy."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,json,time,hashlib
from pathlib import Path
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
D=Path(__file__).resolve().parent;sys.path.insert(0,str(D.parent))
import run_revision as r
from revised_algorithm import revised_scene
P=json.loads((D/'protocol.json').read_text());OUT=D/'results'
def run(job):
 sid,seed,variant=job;file=OUT/'runs'/f'{variant}-{sid}-s{seed}.json'
 if file.exists():
  row=json.loads(file.read_text());assert row['protocol']==P['fingerprint'];return row['run_id']
 sc=r.scene_for(sid,dict(group='validation',overrides={}))
 if variant=='IP-AD-W07':sc=revised_scene(sc)
 else:sc['config']['warning']=0.8
 t=time.perf_counter();p=r.Physics(sc);pre=time.perf_counter()-t
 t=time.perf_counter();cpu=time.process_time()
 ans=r.optimize(p,seed=seed,budget=48,mode='fixed_hold' if variant=='IP-F20-H' else 'hold',deadline_priority=True)
 wall=time.perf_counter()-t;cpu=time.process_time()-cpu
 metrics=r.replay(p,ans['plan']);assert metrics['resource_violations']==metrics['hold_violations']==0
 plan=[asdict(v) for v in ans['plan']];name=f'{variant}-{sid}-s{seed}'
 r.write(OUT/'plans'/f'{name}.json',dict(plan=plan,metrics=metrics,genes=ans['genes'],search=ans['search'],protocol=P['fingerprint']))
 t=time.perf_counter();rows=[]
 for i in range(p.N):
  old=r.evaluate(p,i,plan,sid)
  # Parent cache is keyed by full frozen physical inputs and actual calendar.
  # This independent result key additionally binds the new validation protocol.
  sig=r.digest(dict(physical_calendar_signature=old['signature'],validation_protocol=P['fingerprint']))
  c=dict(old);c.update(signature=sig,physical_calendar_signature=old['signature'],validation_protocol=P['fingerprint'])
  r.write(OUT/'chains'/f'{sig}.json',c);rows.append(c)
 first={i:min(v['end'] for v in plan if v['target']==i) for i in set(v['target'] for v in plan)}
 row=dict(run_id=name,variant=variant,scene=sid,seed=seed,target_count=p.N,configuration=asdict(p.c),protocol=P['fingerprint'],physical_protocol=r.PROTOCOL['fingerprint'],position_rmse_m=float(np.sqrt(np.mean([x['time_position_rmse_m']**2 for x in rows]))),velocity_rmse_m_s=float(np.sqrt(np.mean([x['time_velocity_rmse_m_s']**2 for x in rows]))),compliance_pct=100*sum(x['full_interval_pass_count'] for x in rows)/(p.N*200),laser_s=metrics['laser_s'],occupied_s=metrics['occupied_s'],observations=len(plan),timely_targets=sum(k<=p.tick(120) for k in first.values()),late_targets=sum(k>p.tick(120) for k in first.values()),unobserved_targets=p.N-len(first),planning_wall_s=wall,planning_cpu_s=cpu,preprocessing_s=pre,evaluation_s=time.perf_counter()-t,target_results=[x['signature'] for x in rows],truth_hashes=[x['truth_sha256'] for x in rows])
 r.write(file,row);print(name,'RMSE',round(row['position_rmse_m'],1),'compliance',round(row['compliance_pct'],2),flush=True);return name
if __name__=='__main__':
 for name,h in P['code_hashes'].items():assert hashlib.sha256((D/name).read_bytes()).hexdigest()==h,name
 jobs=[(s,k,v) for s in P['scenes'] for k in P['search_seeds'] for v in P['variants']]
 with ProcessPoolExecutor(max_workers=4) as pool:
  fs=[pool.submit(run,j) for j in jobs]
  for i,f in enumerate(as_completed(fs),1):f.result();print('VALIDATION',i,'/',len(jobs),flush=True)
