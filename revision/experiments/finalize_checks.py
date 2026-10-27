"""Final numerical, paired-input, cache and archive checks; complete outputs required."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,hashlib,time
from pathlib import Path
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor
from run_revision import *
from types import SimpleNamespace
from space_cluster_scheduling.model import Config

def baseline_check(pair):
 from j2_revision import mask_calendar
 sid,seed=pair;p=Physics(scene_for(sid,PROTOCOL['settings'][1]));out=[]
 for method in ['IP-AD','IP-F20-H']:
  name=f'hold-200-{sid}-s{seed}-{method}';file=OUT/'runs'/f'{name}.json';row=json.loads(file.read_text())
  if row.get('full_default_plan_reproduction') and row.get('unperturbed_execution_mask_exact'):continue
  t=time.perf_counter();cpu=time.process_time();a=optimize(p,seed=seed,budget=48,mode='hold' if method=='IP-AD' else 'fixed_hold',deadline_priority=True);elapsed=time.perf_counter()-t;cput=time.process_time()-cpu
  old=json.loads((ARCHIVE/'plans'/f'P1-{sid}-s{seed}-b48-{method}.json').read_text());assert [asdict(v) for v in a['plan']]==old['plan']
  kept,rejected=mask_calendar(p,old['plan'])
  assert not rejected and kept==sorted(old['plan'],key=lambda v:(v['end'],v['satellite'])),name
  planpath=OUT/'plans'/f'{name}.json';pd=json.loads(planpath.read_text());pd.update(archived_search=pd.get('archived_search',pd['search']),search=a['search'],full_default_plan_reproduction=True);write(planpath,pd)
  row.update(archived_planning_wall_s=row.get('archived_planning_wall_s',row['planning_wall_s']),archived_planning_cpu_s=row.get('archived_planning_cpu_s',row['planning_cpu_s']),planning_wall_s=elapsed,planning_cpu_s=cput,full_default_plan_reproduction=True,unperturbed_execution_mask_exact=True)
  write(file,row);out.append(name)
 return out

def verify():
 files=list((OUT/'runs').glob('*.json'));assert len(files)==960,len(files)
 rows=[json.loads(f.read_text()) for f in files];truths={};chains={};settings={c['id']:c for c in PROTOCOL['settings']}
 expected={(c['id'],sid,seed,m) for c in PROTOCOL['settings'] for sid in PROTOCOL['scenes'] for seed in PROTOCOL['search_seeds'] for m in c['methods']}
 assert {(r['setting'],r['scene'],r['seed'],r['method']) for r in rows}==expected
 for r in rows:
  assert r['protocol']==PROTOCOL['fingerprint'] and r['target_count']==r['timely_targets']+r['late_targets']+r['unobserved_targets']
  for i,h in enumerate(r['truth_hashes']):
   key=(r['truth_scene'],i)
   if key in truths:assert truths[key]==h
   truths[key]=h
  sc=scene_for(r['scene'],settings[r['setting']]);stub=SimpleNamespace(c=Config(**sc['config']),scene=sc)
  r['effective_configuration']=asdict(stub.c);write(OUT/'runs'/f"{r['run_id']}.json",r)
  for k,v in settings[r['setting']]['overrides'].items():assert json.loads(json.dumps(r['effective_configuration'][k]))==v
  pd=json.loads((OUT/'plans'/f"{r['run_id']}.json").read_text());assert pd['search']['evaluations']==48 and pd['search']['seed']==r['seed']
  assert pd['metrics']['resource_violations']==0 and pd['metrics']['hold_violations']==0
  rr=[]
  for i,key in enumerate(r['target_results']):
   f=OUT/'chains'/f'{key}.json';assert f.exists()
   if key not in chains:chains[key]=json.loads(f.read_text())
   c=chains[key];rr.append(c);assert c['repetitions']==200 and c['truth_sha256']==r['truth_hashes'][i]
   calendar=[{k:v[k] for k in ['start','end','satellite']} for v in pd['plan'] if v['target']==i]
   assert key==digest(dict(physical=physical_target(stub,i,r['truth_scene']),measurements=calendar))
   assert c['full_interval_pass_count']==int(np.sum((np.array(c['max_position_error_per_realization_m'])<=500)&(np.array(c['max_velocity_error_per_realization_m_s'])<=25)))
  assert np.isfinite(r['position_rmse_m']) and np.isfinite(r['velocity_rmse_m_s'])
  assert np.isclose(r['position_rmse_m'],np.sqrt(np.mean([c['time_position_rmse_m']**2 for c in rr])),rtol=1e-13)
  assert abs(r['compliance_pct']-100*sum(c['full_interval_pass_count'] for c in rr)/(r['target_count']*200))<1e-10
 for sid in PROTOCOL['scenes']:
  cs=[scene_for(sid,settings[f'population-{n}']) for n in [24,48,72,96]]
  for sc in cs:
   for k in ['targets','initial_sigma_r','initial_sigma_v']:assert sc[k]==cs[-1][k][:len(sc['targets'])]
   assert sc['observers']==cs[-1]['observers']
 changed={}
 for setting in PROTOCOL['settings']:
  if setting['group']!='heuristic':continue
  count=0
  for sid in PROTOCOL['scenes']:
   for seed in PROTOCOL['search_seeds']:
    a=json.loads((OUT/'plans'/f"{setting['id']}-{sid}-s{seed}-IP-AD.json").read_text())['plan']
    b=json.loads((OUT/'plans'/f'hold-200-{sid}-s{seed}-IP-AD.json').read_text())['plan']
    count+=a!=b
  changed[setting['id']]=count
 assert any(v>0 for k,v in changed.items() if k.startswith('weights'))
 write(ROOT/'checks/heuristic-calendar-effects.json',dict(paired_calendars_per_variant=40,changed_plan_counts=changed,interpretation='Same physical scene and search seed, one configured heuristic changed. Zero changes are retained and do not invalidate the independent candidate-interface checks.'))
 jrows=[json.loads(f.read_text()) for f in (OUT/'j2/runs').glob('*.json')];assert len(jrows)==160,len(jrows)
 from j2_revision import mask_calendar,run_accuracy,inspect
 jengine=digest(dict(j2=hashlib.sha256((ROOT/'j2_revision.py').read_bytes()).hexdigest(),original=inspect.getsource(run_accuracy.chain)))
 jexpected={(sid,seed,m,mode) for sid in PROTOCOL['scenes'] for seed in PROTOCOL['search_seeds'] for m in ['IP-AD','IP-F20-H'] for mode in ['J2-two-body','J2-J2']}
 assert {(r['scene'],r['seed'],r['method'],r['mode']) for r in jrows}==jexpected
 jt={}
 for r in jrows:
  assert r['engine']==jengine and r['protocol']==PROTOCOL['fingerprint']
  assert r['target_count']==r['timely_targets']+r['late_targets']+r['unobserved_targets']
  jp=json.loads((OUT/'j2/plans'/f"hold-200-{r['scene']}-s{r['seed']}-{r['method']}.json").read_text())
  assert r['execution_mask']==jp['mask_sha256']==digest(jp['plan'])
  assert r['observations']+r['canceled_visits']==jp['original_visits']
  jscene=scene_for(r['scene'],settings['hold-200'])
  rr=[json.loads((OUT/'j2/chains'/f'{key}.json').read_text()) for key in r['target_results']]
  for i,(key,c) in enumerate(zip(r['target_results'],rr)):
   calendar=[{k:v[k] for k in ['start','end','satellite']} for v in jp['plan'] if v['target']==i]
   assert c['measurement_schedule']==calendar
   assert key==c['signature']==digest(dict(scene=jscene,target=i,calendar=calendar,mode=r['mode'],engine=jengine,protocol=PROTOCOL['fingerprint'],reps=200))
  assert all(c['repetitions']==200 for c in rr)
  assert np.isfinite(r['position_rmse_m']) and np.isfinite(r['velocity_rmse_m_s'])
  assert np.isclose(r['position_rmse_m'],np.sqrt(np.mean([c['time_position_rmse_m']**2 for c in rr])),rtol=1e-13)
  assert abs(r['compliance_pct']-100*sum(c['full_interval_pass_count'] for c in rr)/(r['target_count']*200))<1e-10
  for c in rr:
   assert c['full_interval_pass_count']==int(np.sum((np.array(c['max_position_error_per_realization_m'])<=500)&(np.array(c['max_velocity_error_per_realization_m_s'])<=25)))
  for i,h in enumerate(r['truth_hashes']):
   key=(r['scene'],i)
   if key in jt:assert jt[key]==h
   jt[key]=h
 for sid in PROTOCOL['scenes']:
  for seed in PROTOCOL['search_seeds']:
   for method in ['IP-AD','IP-F20-H']:
    rr=[r for r in jrows if (r['scene'],r['seed'],r['method'])==(sid,seed,method)];assert len(rr)==2
    assert rr[0]['execution_mask']==rr[1]['execution_mask'];assert rr[0]['truth_hashes']==rr[1]['truth_hashes']
 assert sum(r.get('unperturbed_execution_mask_exact',False) for r in rows)==80
 write(ROOT/'checks/full-matrix.json',dict(status='PASS',scheduling_runs=len(rows),j2_replays=len(jrows),unique_two_body_chains=len(chains),common_truth_pairs=len(truths),nested_populations=True,j2_common_masks=True,all_targets_retained=True,unperturbed_execution_masks_exact=80,full_default_plan_reproductions=sum(r.get('full_default_plan_reproduction',False) for r in rows)))
 print('FULL MATRIX PASS',len(rows),len(jrows),flush=True)
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--baseline',action='store_true');args=ap.parse_args()
 if args.baseline:
  with ProcessPoolExecutor(max_workers=4) as pool:
   for x in pool.map(baseline_check,[(s,k) for s in PROTOCOL['scenes'] for k in PROTOCOL['search_seeds']]):print('BASELINE CHECK',x,flush=True)
 else:verify()
