"""Predeclared paired controls, fixed-duration sweep and search-budget study.

Run from a checkout with PYTHONPATH=src. Each target retains 200 complete
nonlinear Monte Carlo trajectories. Cached evaluation is keyed by the physical
target schedule, scene, estimator code and repetition count, never by outcome.
"""
from __future__ import annotations
import argparse, hashlib, json, os, platform, time
from pathlib import Path
from dataclasses import asdict, replace
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from .paths import ROOT, DATA
from .scenes import load_scene
from .model import Physics
from .scheduler import optimize, replay
from .run_accuracy import chain

OUT=ROOT/'results/controls'
SEEDS=list(range(3101,3106))
SCENES=list(range(2101,2109))
METHODS={
    'F20':dict(mode='fixed',deadline_priority=False,fixed_dwell=20.),
    'IP-F20':dict(mode='fixed',deadline_priority=True,fixed_dwell=20.),
    'F20-H':dict(mode='fixed_hold',deadline_priority=False,fixed_dwell=20.),
    'IP-F20-H':dict(mode='fixed_hold',deadline_priority=True,fixed_dwell=20.),
    'AD':dict(mode='hold',deadline_priority=False,fixed_dwell=20.),
    'IP-AD':dict(mode='hold',deadline_priority=True,fixed_dwell=20.),
    **{f'IP-F{d}-H':dict(mode='fixed_hold',deadline_priority=True,fixed_dwell=float(d)) for d in (10,30,40)},
}
CORE=list(METHODS)[:6]
CONDITIONS={'P1':{},'P4M1':dict(epsilon_r=2000.,epsilon_v=100.,sigma_angle=50e-6,sigma_range=30.)}

def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2,default=lambda x:x.tolist() if hasattr(x,'tolist') else str(x)))
    tmp.replace(path)

def jobs():
    result=[]
    for cond in CONDITIONS:
        for scene in SCENES:
            for seed in SEEDS:
                methods=CORE+(['IP-F10-H','IP-F30-H','IP-F40-H'] if cond=='P1' and seed==3101 else [])
                result.append(dict(condition=cond,scene=scene,seed=seed,budget=48,methods=methods))
    for scene in (2101,2106):
        for seed in SEEDS:
            for budget in (24,96):
                result.append(dict(condition='P1',scene=scene,seed=seed,budget=budget,methods=CORE))
    return result

def evaluate_target(p,scene_id,target,plan,physical_hash,estimator_hash):
    target_plan=[v for v in plan if v['target']==target]
    # Only measurement times/satellite affect a fixed target's EKF replay.
    schedule=[{k:v[k] for k in ('start','end','satellite')} for v in target_plan]
    signature=digest(dict(scene=physical_hash,target=target,schedule=schedule,repetitions=200,estimator=estimator_hash))
    path=OUT/'chains'/f'{signature}.json'
    if path.exists():return json.loads(path.read_text())
    r=chain(p,target,target_plan,scene_id,200)
    r.pop('curves_1s')
    r.update(signature=signature,physical_scene_sha256=physical_hash,schedule=schedule)
    write(path,r)
    return r

def run_job(job):
    sc=load_scene(DATA/'scenes'/f"{job['scene']}.json")
    sc['config'].update(CONDITIONS[job['condition']])
    t=time.perf_counter();p=Physics(sc);preprocess=time.perf_counter()-t
    physical={**sc,'config':{k:v for k,v in sc['config'].items() if k not in ('epsilon_r','epsilon_v','fixed_dwell')}}
    physical_hash=digest(physical)
    estimator_hash=digest({n:hashlib.sha256((Path(__file__).parent/n).read_bytes()).hexdigest() for n in ('nonlinear_math.py','run_accuracy.py','model.py')})
    rows=[]
    for method in job['methods']:
        run_id=f"{job['condition']}-{job['scene']}-s{job['seed']}-b{job['budget']}-{method}"
        result_file=OUT/'runs'/f'{run_id}.json'
        if result_file.exists():
            rows.append(json.loads(result_file.read_text()));continue
        cfg=METHODS[method];p.c=replace(p.c,fixed_dwell=cfg['fixed_dwell'])
        start=time.perf_counter();cpu=time.process_time()
        ans=optimize(p,seed=job['seed'],budget=job['budget'],mode=cfg['mode'],deadline_priority=cfg['deadline_priority'])
        plan_wall=time.perf_counter()-start;plan_cpu=time.process_time()-cpu
        audit=replay(p,ans['plan'])
        assert audit['objective']==ans['metrics']['objective']
        assert audit['resource_violations']==0
        if cfg['mode']!='fixed':assert audit['hold_violations']==0
        if cfg['mode'] in ('fixed','fixed_hold'):
            assert all(v.end-v.start==p.tick(cfg['fixed_dwell']) for v in ans['plan'])
        plan=[asdict(v) for v in ans['plan']]
        write(OUT/'plans'/f'{run_id}.json',dict(run_id=run_id,job=job,method=method,configuration=cfg,
               plan=plan,metrics=audit,search=ans['search'],genes=ans['genes'],planning_wall_s=plan_wall,
               planning_cpu_s=plan_cpu,preprocess_wall_s=preprocess))
        start=time.perf_counter()
        rr=[evaluate_target(p,job['scene'],i,plan,physical_hash,estimator_hash) for i in range(p.N)]
        row=dict(run_id=run_id,condition=job['condition'],scene=job['scene'],seed=job['seed'],budget=job['budget'],
                 method=method,target_count=p.N,repetitions_per_target=200,
                 position_rmse_m=float(np.sqrt(np.mean([r['time_position_rmse_m']**2 for r in rr]))),
                 velocity_rmse_m_s=float(np.sqrt(np.mean([r['time_velocity_rmse_m_s']**2 for r in rr]))),
                 compliance_pct=100*sum(r['full_interval_pass_count'] for r in rr)/(p.N*200),
                 laser_s=audit['laser_s'],occupied_s=audit['occupied_s'],hold_violations=audit['hold_violations'],
                 observations=len(plan),initialized_targets=audit['initialized_targets'],persistent_targets=audit['persistent_targets'],
                 observed_targets=len(set(v['target'] for v in plan)),planning_wall_s=plan_wall,planning_cpu_s=plan_cpu,
                 preprocess_wall_s=preprocess,evaluation_or_cache_wall_s=time.perf_counter()-start,
                 target_results=[r['signature'] for r in rr],truth_hashes=[r['truth_sha256'] for r in rr])
        write(result_file,row);rows.append(row)
        print(run_id,'Er',round(row['position_rmse_m'],1),'Ev',round(row['velocity_rmse_m_s'],2),
              'plan_s',round(plan_wall,2),'eval_s',round(row['evaluation_or_cache_wall_s'],2),flush=True)
    return rows

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--smoke',action='store_true');args=parser.parse_args()
    js=jobs()
    if args.smoke:js=[dict(condition='P1',scene=2101,seed=3101,budget=48,methods=CORE)]
    source_hashes={n:hashlib.sha256((Path(__file__).parent/n).read_bytes()).hexdigest() for n in ('model.py','scheduler.py','run_accuracy.py','nonlinear_math.py','run_controls.py')}
    protocol=dict(seeds=SEEDS,scenes=SCENES,conditions=CONDITIONS,methods=METHODS,
        primary_budget=48,scan_durations_s=[10,20,30,40],scan_seed=3101,
        budget_sensitivity=dict(scenes=[2101,2106],condition='P1',budgets=[24,48,96],seeds=SEEDS),
        repetitions_per_target=200,workers=args.workers,blas_threads=1,python=platform.python_version(),numpy=np.__version__,
        cpu_model=next((line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),'unknown'),
        logical_cpus=os.cpu_count(),source_sha256=source_hashes,
        timing='perf_counter around optimize only; cold dwell cache per run; Physics preprocessing and independent replay/EKF excluded; concurrent worker wall times and process CPU times both retained',
        fixed_hold_policy='fixed: resource/visibility gates only; fixed_hold: additionally enforce full Eq.(7) at exact prescribed duration; no shortening or extension',
        scan_policy='10,20,30,40 s specified before obtaining new results; descriptive sensitivity, no test-set best-duration selection',
        scope='same eight synthetic orbital configurations; seeds and parameter settings are not independent orbital scenes')
    fingerprint=digest(protocol)
    protocol['fingerprint']=fingerprint
    path=OUT/'protocol.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old['fingerprint']!=fingerprint:raise ValueError('Protocol changed; use fresh output')
    write(path,protocol)
    print('Jobs',len(js),'plans',sum(len(j['methods']) for j in js),'workers',args.workers,flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_job,j) for j in js]
        for done,f in enumerate(as_completed(futures),1):
            f.result();print('Completed jobs',done,'/',len(js),flush=True)
    rows=[json.loads(p.read_text()) for p in sorted((OUT/'runs').glob('*.json'))]
    write(OUT/'all-runs.json',rows)

if __name__=='__main__':main()
