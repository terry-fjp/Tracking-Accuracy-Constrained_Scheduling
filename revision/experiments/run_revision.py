"""Frozen paired revision experiments. No tuning or outcome-based filtering."""
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import sys,json,hashlib,time,copy,argparse
from pathlib import Path
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'src'))
from space_cluster_scheduling.model import Physics,make_scene
from space_cluster_scheduling.scenes import load_scene
from space_cluster_scheduling.scheduler import optimize,replay
from space_cluster_scheduling.run_accuracy import chain
WORKSPACE=ROOT.parents[2];ARCHIVE=Path(os.environ.get('REVISION_REFERENCE_ARCHIVE',str(ROOT/'reference-archive')))
OUT=ROOT/'results';PROTOCOL=json.loads((ROOT/'protocol.json').read_text())
CORE=ROOT/'src/space_cluster_scheduling'
EVAL_HASH=hashlib.sha256(b''.join((CORE/n).read_bytes() for n in ['run_accuracy.py','nonlinear_math.py','evaluation.py'])).hexdigest()

def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def write(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+f'.{os.getpid()}.tmp');tmp.write_text(json.dumps(x,ensure_ascii=False,default=lambda x:x.tolist() if hasattr(x,'tolist') else str(x)));tmp.replace(p)
def scene_for(sid,setting):
 sc=load_scene(ROOT/'data/scenes'/f'{sid}.json')
 if setting['group']=='population':
  sc=make_scene(sid,96,**sc['config']);n=setting['population']
  for k in ['targets','initial_sigma_r','initial_sigma_v']:sc[k]=sc[k][:n]
  sc['seed']=sid+3000
 sc['config'].update(setting['overrides']);return sc

def physical_target(p,i,sid):
 c=p.c
 return dict(scene=sid,target=i,orbit=p.scene['targets'][i],prior=[p.scene['initial_sigma_r'][i],p.scene['initial_sigma_v'][i]],observers=p.scene['observers'],dt=c.dt,horizon=c.horizon,initialization=c.init_deadline,q=c.q,noise=[c.sigma_angle,c.sigma_range],repetitions=200,engine=EVAL_HASH,protocol=PROTOCOL['fingerprint'])
def evaluate(p,i,plan,sid,old=None):
 schedule=[{k:v[k] for k in ['start','end','satellite']} for v in plan if v['target']==i]
 sig=digest(dict(physical=physical_target(p,i,sid),measurements=schedule));file=OUT/'chains'/f'{sig}.json'
 if file.exists():return json.loads(file.read_text())
 if old is not None:
  r=json.loads((ARCHIVE/'chains'/f'{old}.json').read_text());r['imported_source']=str(ARCHIVE/'chains'/f'{old}.json')
 else:r=chain(p,i,[v for v in plan if v['target']==i],sid,200)
 r.pop('curves_1s',None);r['signature']=sig;r['measurement_schedule']=schedule
 write(file,r);return r

def run(job):
 sid,seed,setting=job;sc=scene_for(sid,setting);t=time.perf_counter();p=Physics(sc);pre=time.perf_counter()-t;rows=[]
 for method in setting['methods']:
  name=f"{setting['id']}-{sid}-s{seed}-{method}";out=OUT/'runs'/f'{name}.json'
  if out.exists():
   row=json.loads(out.read_text());assert row['protocol']==PROTOCOL['fingerprint'];rows.append(row);continue
  oldid=f'P1-{sid}-s{seed}-b48-{method}';baseline=setting['id']=='hold-200'
  if baseline:
   oldplan=json.loads((ARCHIVE/'plans'/f'{oldid}.json').read_text());oldrow=json.loads((ARCHIVE/'runs'/f'{oldid}.json').read_text());plan=oldplan['plan'];metrics=oldplan['metrics'];wall=oldrow['planning_wall_s'];cpu=oldrow['planning_cpu_s'];ans={'search':oldplan['search'],'genes':oldplan['genes']}
  else:
   start=time.perf_counter();cput=time.process_time();ans=optimize(p,seed=seed,budget=48,mode='hold' if method=='IP-AD' else 'fixed_hold',deadline_priority=True)
   wall=time.perf_counter()-start;cpu=time.process_time()-cput;metrics=replay(p,ans['plan']);plan=[asdict(v) for v in ans['plan']]
   assert metrics['resource_violations']==metrics['hold_violations']==0
  write(OUT/'plans'/f'{name}.json',dict(run_id=name,plan=plan,metrics=metrics,genes=ans['genes'],search=ans['search'],source=str(ARCHIVE/'plans'/f'{oldid}.json') if baseline else None))
  start=time.perf_counter();rr=[evaluate(p,i,plan,sc['seed'],oldrow['target_results'][i] if baseline else None) for i in range(p.N)]
  first={i:min(v['end'] for v in plan if v['target']==i) for i in set(v['target'] for v in plan)}
  row=dict(run_id=name,setting=setting['id'],group=setting['group'],scene=sid,truth_scene=sc['seed'],seed=seed,method=method,target_count=p.N,configuration=sc['config'],protocol=PROTOCOL['fingerprint'],position_rmse_m=float(np.sqrt(np.mean([r['time_position_rmse_m']**2 for r in rr]))),velocity_rmse_m_s=float(np.sqrt(np.mean([r['time_velocity_rmse_m_s']**2 for r in rr]))),compliance_pct=100*sum(r['full_interval_pass_count'] for r in rr)/(p.N*200),laser_s=metrics['laser_s'],occupied_s=metrics['occupied_s'],observations=len(plan),timely_targets=sum(k<=p.tick(120) for k in first.values()),late_targets=sum(k>p.tick(120) for k in first.values()),unobserved_targets=p.N-len(first),planning_wall_s=wall,planning_cpu_s=cpu,preprocessing_s=pre,evaluation_s=time.perf_counter()-start,archived_baseline=baseline,target_results=[r['signature'] for r in rr],truth_hashes=[r['truth_sha256'] for r in rr])
  if baseline:
   assert abs(row['position_rmse_m']-oldrow['position_rmse_m'])<1e-9
   assert abs(row['compliance_pct']-oldrow['compliance_pct'])<1e-10
  write(out,row);rows.append(row);print(name,'Er',round(row['position_rmse_m'],2),'pass',round(row['compliance_pct'],2),'plan_s',round(wall,1),'eval_s',round(row['evaluation_s'],1),flush=True)
 return len(rows)

def regression():
 results=[]
 for sid,seed in [(2101,3101),(2106,3105)]:
  p=Physics(scene_for(sid,PROTOCOL['settings'][1]))
  for method in ['IP-AD','IP-F20-H']:
   oldid=f'P1-{sid}-s{seed}-b48-{method}';old=json.loads((ARCHIVE/'plans'/f'{oldid}.json').read_text())
   a=optimize(p,seed=seed,budget=48,mode='hold' if method=='IP-AD' else 'fixed_hold',deadline_priority=True);plan=[asdict(v) for v in a['plan']]
   assert plan==old['plan'],oldid
   oldrow=json.loads((ARCHIVE/'runs'/f'{oldid}.json').read_text())
   for i in [0,p.N-1]:
    r=chain(p,i,[v for v in plan if v['target']==i],sid,200);oldr=json.loads((ARCHIVE/'chains'/f"{oldrow['target_results'][i]}.json").read_text())
    assert r['truth_sha256']==oldr['truth_sha256']
    assert np.isclose(r['time_position_rmse_m'],oldr['time_position_rmse_m'],rtol=1e-12,atol=1e-8)
    assert r['full_interval_pass_count']==oldr['full_interval_pass_count']
   results.append({'run':oldid,'plan_exact':True,'nonlinear_replay_targets':[0,p.N-1]});print('REGRESSION PASS',oldid,flush=True)
 write(ROOT/'checks/default-regression.json',{'status':'PASS','cases':results,'archived_baseline_policy':'Reuse archived exact default-input plans and target outputs after representative regression and unchanged original evaluation-source verification.'})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--regression',action='store_true');ap.add_argument('--workers',type=int,default=4);ap.add_argument('--group');args=ap.parse_args()
 if args.regression:regression()
 else:
  assert (ROOT/'checks/default-regression.json').exists()
  jobs=[(s,k,c) for s in PROTOCOL['scenes'] for k in PROTOCOL['search_seeds'] for c in PROTOCOL['settings'] if args.group is None or c['group']==args.group]
  write(OUT/'execution.json',dict(protocol=PROTOCOL['fingerprint'],runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),jobs=len(jobs),workers=args.workers))
  with ProcessPoolExecutor(max_workers=args.workers) as pool:
   futures=[pool.submit(run,j) for j in jobs]
   for i,f in enumerate(as_completed(futures),1):f.result();print('JOBS',i,'/',len(jobs),flush=True)
