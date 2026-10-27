"""All-target nonlinear accuracy evaluation with common random truth/noise.

Physical position/velocity errors are the final endpoints. The scheduling
covariance score is never used to classify Monte Carlo trajectory accuracy.
"""
import json, sys, time, platform, hashlib, os
from pathlib import Path
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import numpy as np
from .evaluation import THRESHOLDS, BASELINE, key, rescore_record
from .model import Physics, MU
from .scheduler import optimize, replay

from .paths import ROOT, DATA, REFERENCE
from .scenes import load_scene
OUT=Path(os.environ.get('AEROSPACE_ACCURACY_OUTPUT',REFERENCE))
SOURCE=DATA
from .nonlinear_math import rk4, transition, measurement

def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2,default=lambda a:a.tolist() if hasattr(a,'tolist') else str(a)))

def chain(p,i,plan,scene_id,reps=200):
    # Separate streams prevent a changed observation schedule from changing truth.
    seed=9210000+scene_id*100+i
    streams=[np.random.default_rng(np.random.SeedSequence([seed,j])) for j in range(3)]
    ri,rp,rm=streams;dt=p.c.dt;end=p.tick(p.c.horizon);init=p.tick(p.c.init_deadline)
    truth=p.x[i,0]+ri.normal(size=(reps,6))@np.linalg.cholesky(p.p0[i]).T
    estimate=np.broadcast_to(p.x[i,0],(reps,6)).copy()
    P=np.broadcast_to(p.p0[i],(reps,6,6)).copy();ql=np.linalg.cholesky(p.Q)
    schedule={k:v['satellite'] for v in plan for k in range(v['start']+1,v['end']+1)}
    passed=np.ones(reps,dtype=bool);sum_r2=sum_v2=0.;peak_r=peak_v=0.;count=0
    curves=[];nees=[];truth_digest=hashlib.sha256()
    max_r2=np.zeros(reps);max_v2=np.zeros(reps)
    for k in range(1,end+1):
        truth=rk4(truth,dt)+rp.normal(size=(reps,6))@ql.T
        noise=rm.normal(size=(reps,3)) # Draw also on non-observation steps.
        F=transition(estimate,dt);estimate=rk4(estimate,dt);P=F@P@F.swapaxes(-1,-2)+p.Q
        if k in schedule:
            s=schedule[k];C,R=p.H_R(i,s,k);rho=np.linalg.norm(p.x[i,k,:3]-p.obs[s,k,:3])
            basis=C[:,:3].copy();basis[:2]*=rho
            z,_=measurement(truth,p.obs[s,k],basis);z+=noise*np.sqrt(np.diag(R))
            zhat,H=measurement(estimate,p.obs[s,k],basis)
            S=H@P@H.swapaxes(-1,-2)+R;K=np.linalg.solve(S,H@P).swapaxes(-1,-2)
            residual=z-zhat;residual[:,:2]=(residual[:,:2]+np.pi)%(2*np.pi)-np.pi
            estimate+=np.einsum('nij,nj->ni',K,residual)
            B=np.eye(6)-K@H;P=B@P@B.swapaxes(-1,-2)+K@R@K.swapaxes(-1,-2)
        error=truth-estimate;rr=np.sum(error[:,:3]**2,axis=1);vv=np.sum(error[:,3:]**2,axis=1)
        rmse_r=float(np.sqrt(rr.mean()));rmse_v=float(np.sqrt(vv.mean()))
        if k>=init:
            max_r2=np.maximum(max_r2,rr);max_v2=np.maximum(max_v2,vv)
            passed &= (rr<=BASELINE[0]**2)&(vv<=BASELINE[1]**2)
            sum_r2+=float(rr.mean());sum_v2+=float(vv.mean());count+=1
            peak_r=max(peak_r,rmse_r);peak_v=max(peak_v,rmse_v)
        if k%5==0:curves.append([k*dt,rmse_r,rmse_v])
        if k in (init,end):
            truth_digest.update(truth.tobytes())
            ne=np.einsum('ni,ni->n',error,np.linalg.solve(P,error[...,None])[...,0])
            nees.append({'time':k*dt,'mean_nees':float(ne.mean())})
    return {'target':i,'repetitions':reps,'seed':seed,'truth_sha256':truth_digest.hexdigest(),
            'time_position_rmse_m':(sum_r2/count)**.5,'time_velocity_rmse_m_s':(sum_v2/count)**.5,
            'peak_position_rmse_m':peak_r,'peak_velocity_rmse_m_s':peak_v,
            'final_position_rmse_m':rmse_r,'final_velocity_rmse_m_s':rmse_v,
            'full_interval_pass_count':int(passed.sum()),'full_interval_pass_fraction':float(passed.mean()),
            'max_position_error_per_realization_m':np.sqrt(max_r2),
            'max_velocity_error_per_realization_m_s':np.sqrt(max_v2),
            'threshold_pass_pct':{f'{er:g}/{ev:g}':100*float(np.mean((max_r2<=er**2)&(max_v2<=ev**2)))
                for er,ev in THRESHOLDS},
            'curves_1s':curves,'nees_checkpoints':nees,'observation_count':len(plan)}

def guard_existing_output(output, fingerprint):
    previous=Path(output)/'protocol.json'
    if previous.exists():
        if json.loads(previous.read_text()).get('input_fingerprint')!=fingerprint:
            raise ValueError('Output contains a different experiment; choose a fresh output directory')
    elif any((Path(output)/name).exists() for name in ['plans','chains']):
        raise ValueError('Existing plans/chains lack a protocol; choose a fresh output directory')


def run_scene(scene_id, output, protocol):
    OUT=Path(output)
    summaries=[]
    model_scene=load_scene(SOURCE/'scenes'/f'{scene_id}.json')
    model_scene['config'].update(protocol.get('config_overrides',{}))
    scene_hash=hashlib.sha256(json.dumps(model_scene,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    p=Physics(model_scene)
    truth_hashes={}
    for method in protocol['methods']:
        planfile=OUT/'plans'/f'{scene_id}-{method}.json'
        if planfile.exists():
            data=json.loads(planfile.read_text())
            if data.get('scene_sha256')!=scene_hash:
                raise ValueError('Cached plan does not match current scene; use a fresh output directory')
        else:
            ans=optimize(p,seed=3101,budget=48,mode='fixed' if method=='Fixed20' else 'hold',
                         deadline_priority=(method=='Deadline'))
            audit=replay(p,ans['plan'])
            assert audit['objective']==ans['metrics']['objective']
            assert audit['resource_violations']==0
            if method!='Fixed20': assert audit['hold_violations']==0
            data={'scene':scene_id,'scene_sha256':scene_hash,'hold_s':p.c.hold,
                  'plan':[asdict(v) for v in ans['plan']],'metrics':audit,
                  'search':ans['search'],'genes':ans['genes']}
            write(planfile,data)
            print(protocol.get('condition_id',''),scene_id,method,'new plan',len(data['plan']),'observations',flush=True)
        plan_hash=hashlib.sha256(json.dumps(data['plan'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
        rr=[];t=time.perf_counter()
        for i in range(p.N):
            outfile=OUT/'chains'/f'{scene_id}-{method}-{i}.json'
            result=json.loads(outfile.read_text()) if outfile.exists() else None
            if result is not None and (result.get('scene_sha256')!=scene_hash or result.get('plan_sha256')!=plan_hash):
                raise ValueError('Cached chain input mismatch; use a fresh output directory')
            if result is None or 'threshold_pass_pct' not in result:
                previous=result
                result=chain(p,i,[v for v in data['plan'] if v['target']==i],scene_id,protocol['repetitions_per_target'])
                if previous:
                    assert result['truth_sha256']==previous['truth_sha256']
                    assert abs(result['time_position_rmse_m']-previous['time_position_rmse_m'])<1e-9
                result.update(scene_sha256=scene_hash,plan_sha256=plan_hash,hold_s=p.c.hold)
                write(outfile,result)
            result=rescore_record(result);write(outfile,result)
            if i in truth_hashes:assert result['truth_sha256']==truth_hashes[i],(scene_id,method,i)
            truth_hashes[i]=result['truth_sha256'];rr.append(result)
        curves=np.array([r['curves_1s'] for r in rr]);mean_curve=np.sqrt((curves[:,:,1:]**2).mean(axis=0))
        summary={'scene':scene_id,'method':method,'target_count':p.N,'repetitions':200,'target_trajectories':p.N*200,
                 'position_rmse_m':float(np.sqrt(np.mean([r['time_position_rmse_m']**2 for r in rr]))),
                 'velocity_rmse_m_s':float(np.sqrt(np.mean([r['time_velocity_rmse_m_s']**2 for r in rr]))),
                 'final_position_rmse_m':float(np.sqrt(np.mean([r['final_position_rmse_m']**2 for r in rr]))),
                 'final_velocity_rmse_m_s':float(np.sqrt(np.mean([r['final_velocity_rmse_m_s']**2 for r in rr]))),
                 'actual_full_interval_pass_pct':100*sum(r['full_interval_pass_count'] for r in rr)/(p.N*200),
                 'threshold_pass_pct':{key:float(np.mean([r['threshold_pass_pct'][key] for r in rr])) for key in rr[0]['threshold_pass_pct']},
                 'laser_s':data['metrics']['laser_s'],'observations':len(data['plan']),
                 'mean_curve_1s':np.column_stack([curves[0,:,0],mean_curve]),'seconds':time.perf_counter()-t}
        summaries.append(summary);write(OUT/'scene-summaries'/f'{scene_id}.json',summaries)
        print(protocol.get('condition_id',''),scene_id,method,'RMSE',round(summary['position_rmse_m'],2),round(summary['velocity_rmse_m_s'],3),
              'pass%',round(summary['actual_full_interval_pass_pct'],2),'sec',round(summary['seconds'],1),flush=True)
    del p
    return summaries


def prepare_protocol(output, condition, workers=1):
    OUT=Path(output)
    protocol={'scenes':list(range(2101,2109)),'methods':['Fixed20','Hold','Deadline'],
              'solver':'H-SADE','solver_seed':3101,'complete_plan_evaluations':48,
              'repetitions_per_target':200,'truth_seed_rule':'SeedSequence([9210000 + scene*100 + target, stream]); streams 0=initial,1=process,2=measurement',
              'truth':'Cartesian two-body RK4, dt=0.2 s, additive process Q, nonlinear tangent angles and true range',
              'estimator':'EKF, relinearized each time; Joseph update; wrapped angular residuals',
              'evaluation_interval_s':[120,900],'all_targets_in_denominator':True,
              'primary_endpoints':['position RMSE (m)','velocity RMSE (m/s)'],
              'secondary_endpoint':f'fraction of target-realization pairs whose actual position error <={BASELINE[0]:g} m and actual velocity error <={BASELINE[1]:g} m/s at every maintenance grid',
              'paired_randomness':'identical truth and standardized noise for each scene/target across all methods',
              'assumptions':'perfect observer ephemerides and association; nominally planned observations always available',
              'evaluation_baseline':list(BASELINE),
              'evaluation_threshold_pairs':[list(pair) for pair in THRESHOLDS],
              'threshold_sweep':'Each condition has its own newly optimized plans; all four actual-error acceptance thresholds are evaluated on each trajectory batch.',
              'condition_id':condition['id'],'config_overrides':condition['overrides'],
              'planning_thresholds':[condition['planning_position_m'],condition['planning_velocity_m_s']],
              'measurement_standard_deviations':[condition['sigma_angle_urad']/1e6,condition['sigma_range_m']],
              'hold_s':200.0,'plan_policy':'All three methods newly optimized; no old baseline reuse',
              'workers':workers,'python':platform.python_version(),'numpy':np.__version__}
    inputs={str(sid):load_scene(SOURCE/'scenes'/f'{sid}.json') for sid in protocol['scenes']}
    for sc in inputs.values():sc['config'].update(condition['overrides'])
    assert all(sc['config']['hold']==protocol['hold_s'] for sc in inputs.values())
    core={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
          for name in ['model.py','scheduler.py','nonlinear_math.py','run_accuracy.py']}
    fingerprint=hashlib.sha256(json.dumps({'scenes':inputs,'core':core,'seed':3101,'budget':48,
                             'repetitions':200,'numpy':np.__version__},sort_keys=True).encode()).hexdigest()
    protocol.update(input_fingerprint=fingerprint,core_sha256=core)
    guard_existing_output(OUT,fingerprint)

    write(OUT/'protocol.json',protocol)
    return protocol


def main(workers=1, condition=None):
    if condition is None:
        spec=json.loads((ROOT/'config/experiments.json').read_text())
        condition=next(c for c in spec['conditions'] if c['id']==spec['reference_condition'])
    protocol=prepare_protocol(OUT,condition,workers)
    summaries=[]
    job=partial(run_scene,output=OUT,protocol=protocol)
    if workers==1:
        batches=map(job,protocol['scenes'])
        for batch in batches:
            summaries.extend(batch);write(OUT/'scene-results.json',summaries)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for batch in pool.map(job,protocol['scenes']):
                summaries.extend(batch);write(OUT/'scene-results.json',summaries)
    summarize(OUT,protocol,summaries)


def summarize(output,protocol,summaries):
    OUT=Path(output)
    write(OUT/'scene-results.json',summaries)
    total=[]
    for method in protocol['methods']:
        rr=[r for r in summaries if r['method']==method]
        fields=['position_rmse_m','velocity_rmse_m_s','final_position_rmse_m','final_velocity_rmse_m_s','actual_full_interval_pass_pct','laser_s']
        total.append({'method':method,'scenes':len(rr),'aggregation':'equal-weight mean of 8 scene metrics',
                      **{f:float(np.mean([r[f] for r in rr])) for f in fields},
                      'threshold_pass_pct':{key:float(np.mean([r['threshold_pass_pct'][key] for r in rr])) for key in rr[0]['threshold_pass_pct']}})
    write(OUT/'summary.json',total)
    print(json.dumps(total,ensure_ascii=False),flush=True)

if __name__=='__main__':main()
