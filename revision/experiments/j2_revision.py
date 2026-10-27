"""J2 deterministic propagation and fixed-calendar paired tracking sensitivity."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,json,hashlib,math,time,copy,argparse
from pathlib import Path
from types import SimpleNamespace
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from run_revision import ROOT,OUT,PROTOCOL,digest,write,scene_for,ARCHIVE
from space_cluster_scheduling.model import MU,RE,Config,Physics,circular
from space_cluster_scheduling.nonlinear_math import rk4 as rk4_two,transition as transition_two
from space_cluster_scheduling import run_accuracy
from space_cluster_scheduling.scheduler import pointing,same_field,slew
J2=PROTOCOL['J2']['constant']

def acceleration(r,j2=J2):
 rho=np.linalg.norm(r,axis=-1);out=-MU*r/rho[...,None]**3
 if j2:
  b=np.full_like(r,-1.);b[...,2]=-3.;b+=5*(r[...,2]**2/rho**2)[...,None]
  out+=1.5*j2*MU*RE**2*r*b/rho[...,None]**5
 return out

def gradient(r,j2=J2):
 rho=np.linalg.norm(r,axis=-1);rr=r[..., :,None]*r[...,None,:]
 g=MU*(3*rr/rho[...,None,None]**5-np.eye(3)/rho[...,None,None]**3)
 if j2:
  b=np.full_like(r,-1.);b[...,2]=-3.;b+=5*(r[...,2]**2/rho**2)[...,None]
  db=-10*r[...,2,None]**2*r/rho[...,None]**4;db[...,2]+=10*r[...,2]/rho**2
  g+=1.5*j2*MU*RE**2*(np.eye(3)*b[..., :,None]/rho[...,None,None]**5+r[..., :,None]*db[...,None,:]/rho[...,None,None]**5-5*(r*b)[..., :,None]*r[...,None,:]/rho[...,None,None]**7)
 return g

def rhs(x,j2=J2):return np.concatenate((x[...,3:],acceleration(x[...,:3],j2)),axis=-1)
def rk4(x,dt,j2=J2):
 a=rhs(x,j2);b=rhs(x+dt*a/2,j2);c=rhs(x+dt*b/2,j2);d=rhs(x+dt*c,j2)
 return x+dt*(a+2*b+2*c+d)/6

def transition(x,dt,j2=J2):
 g=gradient(x[...,:3],j2);f=np.broadcast_to(np.eye(6),x.shape[:-1]+(6,6)).copy()
 f[...,:3,3:]=dt*np.eye(3);f[...,3:,:3]=dt*g
 f[...,:3,:3]+=.5*dt**2*g;f[...,3:,3:]+=.5*dt**2*g
 return f

def propagate(initial,dt,steps,j2=J2):
 arr=np.empty((len(initial),steps+1,6));arr[:,0]=initial
 for k in range(steps):arr[:,k+1]=rk4(arr[:,k],dt,j2)
 return arr

def visible(x,obs,c):
 result=np.empty((len(x),len(obs),x.shape[1]),bool)
 for i in range(len(x)):
  d=x[i,None,:,:3]-obs[:,:,:3];rho=np.linalg.norm(d,axis=-1);o=obs[:,:,:3]
  u=d/rho[...,None];nadir=-o/np.linalg.norm(o,axis=-1)[...,None]
  angle=np.degrees(np.arccos(np.clip(np.sum(u*nadir,axis=-1),-1,1)))
  frac=np.clip(-np.sum(o*d,axis=-1)/rho**2,0,1)
  clear=np.linalg.norm(o+frac[...,None]*d,axis=-1)>=RE
  dv=x[i,None,:,3:]-obs[:,:,3:];angular=np.linalg.norm(np.cross(d,dv),axis=-1)/rho**2*180/np.pi
  result[i]=(rho<=c.max_range)&(angle<=c.off_nadir)&clear&(angular<=c.omega)
 return result

class Geometry:
 tick=Physics.tick;ceil_tick=Physics.ceil_tick;H_R=Physics.H_R
 def __init__(self,scene):
  self.scene=scene;self.c=Config(**scene['config']);self.N=len(scene['targets']);self.S=len(scene['observers']);dt=self.c.dt
  times=np.arange(round(self.c.horizon/dt)+1)*dt
  self.two_x=np.array([circular(*e,times) for e in scene['targets']]);self.two_obs=np.array([circular(*e,times) for e in scene['observers']])
  allx=propagate(np.concatenate([self.two_x[:,0],self.two_obs[:,0]]),dt,len(times)-1)
  self.x=allx[:self.N];self.obs=allx[self.N:];self.visible=visible(self.x,self.obs,self.c)
  self.p0=np.array([np.diag([r*r]*3+[v*v]*3) for r,v in zip(scene['initial_sigma_r'],scene['initial_sigma_v'])])
  self.Q=self.c.q*np.block([[dt**3/3*np.eye(3),dt**2/2*np.eye(3)],[dt**2/2*np.eye(3),dt*np.eye(3)]])

def mask_calendar(p,plan):
 """Keep original timestamps; assess each visit against the last retained visit."""
 kept=[];rejected=[]
 for s in range(p.S):
  prev=None;used=0
  for v in sorted([SimpleNamespace(**v) for v in plan if v['satellite']==s],key=lambda v:v.start):
   reasons=[]
   if not p.visible[v.target,s,v.capture:v.end+1].all():reasons.append('visibility')
   if v.end-v.start>p.tick(p.c.continuous_limit):reasons.append('continuous_limit')
   if used+v.end-v.start>p.tick(p.c.laser_budget):reasons.append('laser_budget')
   shared=bool(v.shared and prev and prev.target!=v.target and v.capture-prev.end<=p.tick(10) and same_field(p,v.target,prev.target,s,prev.end,v.end))
   continued=bool(prev and prev.target==v.target and v.capture==prev.end)
   if v.start-v.capture<(0 if shared or continued else p.tick(p.c.capture)):reasons.append('acquisition')
   if prev:
    angle=math.degrees(math.acos(float(np.clip(pointing(p,prev.target,s,prev.end)@pointing(p,v.target,s,v.capture),-1,1))))
    need=0 if prev.target==v.target else p.ceil_tick(slew(angle,p.c)+p.c.settle)
    if v.capture-prev.end<need:reasons.append('transition_gap')
    if v.capture-v.setup<need:reasons.append('reserved_setup_interval')
    if v.setup<prev.end:reasons.append('overlap')
   if reasons:rejected.append(dict(visit=vars(v),reasons=reasons))
   else:kept.append(vars(v));prev=v;used+=v.end-v.start
 return sorted(kept,key=lambda v:(v['end'],v['satellite'])),rejected

# Preserve every random draw and measurement update from the submitted evaluator;
# only explicitly parameterize truth and estimator dynamics in a separate function.
import inspect
_source=inspect.getsource(run_accuracy.chain).replace('def chain(', 'def j2_chain(')
_source=_source.replace('truth=rk4(truth,dt)', 'truth=truth_step(truth,dt)').replace('estimate=rk4(estimate,dt)','estimate=estimate_step(estimate,dt)')
_ns=dict(run_accuracy.__dict__);_ns.update(truth_step=rk4,estimate_step=rk4_two,transition=transition_two)
exec(compile(_source,__file__+':adapted_chain','exec'),_ns)

def checks():
 r=np.array([[7e6,1e6,2e6],[-2e6,5e6,-4e6],[7e6,0,0],[0,0,7e6]])
 g=gradient(r);fd=np.empty_like(g)
 for j in range(3):
  d=np.zeros(3);d[j]=1;fd[:,:,j]=(acceleration(r+d)-acceleration(r-d))/2
 rel=float(np.linalg.norm(g-fd)/np.linalg.norm(g));assert rel<1e-8
 x=np.column_stack([r,np.tile([0,6000,3000],(len(r),1))]);assert np.array_equal(rk4(x,.2,0),rk4_two(x,.2));assert np.array_equal(transition(x,.2,0),transition_two(x,.2))
 # Independent potential derivative: a = grad(mu/r * [1-J2 (Re/r)^2 P2(z/r)]).
 def potential(y):
  rad=np.linalg.norm(y,axis=-1);return MU/rad*(1-J2*(RE/rad)**2*(3*(y[...,2]/rad)**2-1)/2)
 af=np.empty_like(r)
 for j in range(3):
  d=np.zeros(3);d[j]=10;af[:,j]=(potential(r+d)-potential(r-d))/20
 ar=float(np.linalg.norm(acceleration(r)-af)/np.linalg.norm(af));assert ar<1e-8
 a=propagate(x,.2,4500)[:,-1];b=propagate(x,.1,9000)[:,-1];er=float(np.max(np.linalg.norm(a[:,:3]-b[:,:3],axis=1)));ev=float(np.max(np.linalg.norm(a[:,3:]-b[:,3:],axis=1)));assert er<1e-3 and ev<1e-5
 write(ROOT/'checks/j2-numerics.json',dict(status='PASS',jacobian_relative_error=rel,acceleration_potential_relative_error=ar,zero_J2_exact=True,step_halving_position_m=er,step_halving_velocity_m_s=ev))
 print('J2 NUMERICAL CHECKS PASS',rel,ar,er,ev,flush=True)

def deterministic(p,sid):
 records={}
 rawdir=OUT/'j2/propagation';rawdir.mkdir(parents=True,exist_ok=True)
 np.savez_compressed(rawdir/f'{sid}.npz',two_body_targets=p.two_x,j2_targets=p.x,two_body_observers=p.two_obs,j2_observers=p.obs,j2_visible=p.visible,dt=p.c.dt)
 for kind,a,b in [('targets',p.x,p.two_x),('observers',p.obs,p.two_obs)]:
  dr=np.linalg.norm(a[:,:,:3]-b[:,:,:3],axis=-1);dv=np.linalg.norm(a[:,:,3:]-b[:,:,3:],axis=-1)
  records[kind]={str(t):dict(position_max_m=float(dr[:,p.tick(t)].max()),position_mean_m=float(dr[:,p.tick(t)].mean()),velocity_max_m_s=float(dv[:,p.tick(t)].max()),velocity_mean_m_s=float(dv[:,p.tick(t)].mean())) for t in [200,900]}
  records[kind]['full_horizon']=dict(position_max_m=float(dr.max()),velocity_max_m_s=float(dv.max()))
 base=visible(p.two_x,p.two_obs,p.c);flip=base!=p.visible
 records['visibility']=dict(changed_samples=int(flip.sum()),all_samples=int(flip.size),baseline_visible=int(base.sum()),j2_visible=int(p.visible.sum()),gained=int((~base&p.visible).sum()),lost=int((base&~p.visible).sum()),checkpoints={str(t):int(flip[:,:,p.tick(t)].sum()) for t in [200,900]})
 write(OUT/'j2/geometry'/f'{sid}.json',records)
 return records

def run_scene(sid):
 sc=scene_for(sid,PROTOCOL['settings'][1]);p=Geometry(sc);deterministic(p,sid)
 engine=digest(dict(j2=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),original=inspect.getsource(run_accuracy.chain)))
 for seed in PROTOCOL['search_seeds']:
  for method in ['IP-AD','IP-F20-H']:
   name=f'hold-200-{sid}-s{seed}-{method}';original=json.loads((ARCHIVE/'plans'/f'P1-{sid}-s{seed}-b48-{method}.json').read_text())['plan']
   plan,dropped=mask_calendar(p,original);mask=digest(plan)
   write(OUT/'j2/plans'/f'{name}.json',dict(plan=plan,canceled=dropped,mask_sha256=mask,original_visits=len(original)))
   truths={}
   for mode in ['J2-two-body','J2-J2']:
    out=OUT/'j2/runs'/f'{name}-{mode}.json'
    if out.exists():
     row=json.loads(out.read_text());assert row['engine']==engine and row['protocol']==PROTOCOL['fingerprint'];continue
    _ns.update(estimate_step=rk4 if mode=='J2-J2' else rk4_two,transition=transition if mode=='J2-J2' else transition_two)
    rows=[];start=time.perf_counter()
    for i in range(p.N):
     visits=[v for v in plan if v['target']==i]
     key=digest(dict(scene=sc,target=i,calendar=[{k:v[k] for k in ['start','end','satellite']} for v in visits],mode=mode,engine=engine,protocol=PROTOCOL['fingerprint'],reps=200))
     file=OUT/'j2/chains'/f'{key}.json'
     if file.exists():r=json.loads(file.read_text())
     else:
      r=_ns['j2_chain'](p,i,visits,sid,200);r.pop('curves_1s',None);r.update(signature=key,engine=engine,origin_calendar_mask=mask,measurement_schedule=[{k:v[k] for k in ['start','end','satellite']} for v in visits]);write(file,r)
     if i in truths:assert truths[i]==r['truth_sha256']
     truths[i]=r['truth_sha256'];rows.append(r)
    first={i:min(v['end'] for v in plan if v['target']==i) for i in set(v['target'] for v in plan)}
    row=dict(run_id=name+'-'+mode,mode=mode,scene=sid,seed=seed,method=method,protocol=PROTOCOL['fingerprint'],engine=engine,execution_mask=mask,position_rmse_m=float(np.sqrt(np.mean([r['time_position_rmse_m']**2 for r in rows]))),velocity_rmse_m_s=float(np.sqrt(np.mean([r['time_velocity_rmse_m_s']**2 for r in rows]))),compliance_pct=100*sum(r['full_interval_pass_count'] for r in rows)/(p.N*200),target_count=p.N,observations=len(plan),canceled_visits=len(dropped),laser_s=sum(v['end']-v['start'] for v in plan)*p.c.dt,occupied_s=sum(v['end']-v['setup'] for v in plan)*p.c.dt,timely_targets=sum(t<=p.tick(120) for t in first.values()),late_targets=sum(t>p.tick(120) for t in first.values()),unobserved_targets=p.N-len(first),evaluation_s=time.perf_counter()-start,target_results=[r['signature'] for r in rows],truth_hashes=[r['truth_sha256'] for r in rows])
    write(out,row);print(row['run_id'],'Er',round(row['position_rmse_m'],2),'pass',round(row['compliance_pct'],2),'cancel',len(dropped),flush=True)
 return sid
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');ap.add_argument('--workers',type=int,default=4);ap.add_argument('--scene',type=int);args=ap.parse_args()
 if args.check:checks()
 else:
  assert (ROOT/'checks/j2-numerics.json').exists()
  with ProcessPoolExecutor(max_workers=args.workers) as pool:
   for x in pool.map(run_scene,[args.scene] if args.scene else PROTOCOL['scenes']):print('SCENE COMPLETE',x,flush=True)
