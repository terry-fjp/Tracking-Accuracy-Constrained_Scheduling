"""Two-body nominal orbits and linearized covariance, SI units throughout.

No window-midpoint covariance labels. Every visit uses its actual prior and time.
Discrete feasibility at 0.2 s; no claimed continuous-time or flight guarantee.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math
import numpy as np

MU=3.986004418e14
RE=6378137.0

@dataclass(frozen=True)
class Config:
    dt: float=0.2
    horizon: float=900.0
    init_deadline: float=120.0
    hold: float=200.0
    epsilon_r: float=500.0
    epsilon_v: float=25.0
    q: float=1e-4
    sigma_angle: float=1e-4
    sigma_range: float=50.0
    max_range: float=2e6
    off_nadir: float=60.0
    fov: float=1.0
    omega: float=5.0
    acceleration: float=2.5
    settle: float=2.0
    capture: float=5.0
    continuous_limit: float=60.0
    laser_budget: float=120.0
    warning: float=0.8
    fixed_dwell: float=20.0
    duration_multipliers: tuple=(1.0,1.5)


def circular(radius,inc,raan,phase,times):
    inc,raan,phase=np.radians([inc,raan,phase])
    n=np.sqrt(MU/radius**3)
    u=phase+n*times
    a=np.array([np.cos(raan),np.sin(raan),0.])
    b=np.array([-np.sin(raan)*np.cos(inc),np.cos(raan)*np.cos(inc),np.sin(inc)])
    r=radius*(np.cos(u)[:,None]*a+np.sin(u)[:,None]*b)
    v=radius*n*(-np.sin(u)[:,None]*a+np.cos(u)[:,None]*b)
    return np.concatenate([r,v],axis=1)


def make_scene(seed=100,n_targets=12,**changes):
    """Independent synthetic circular-orbit families, not TLE replicas.

    All 72 Walker observers are retained. Targets populate two randomly selected
    planes near a selected observer; every target has independent orbital elements.
    No candidate/target is discarded from the reported denominator.
    """
    rng=np.random.default_rng(seed)
    cfg=Config(**changes)
    observers=[]
    phase0=float(rng.uniform(0,360));raan0=float(rng.uniform(0,360))
    for p in range(8):
        for s in range(9):
            observers.append([RE+800e3,60.,raan0+p*45,phase0+s*40+p*5])
    anchors=rng.choice(72,2,replace=False)
    targets=[]
    for i in range(n_targets):
        parent=observers[int(anchors[i%2])]
        # Independently sampled co-orbit families deliberately exercise persistent
        # observation windows; these are conceptual objects, not real Starlinks.
        alt=float(rng.uniform(610e3,690e3))
        targets.append([RE+alt,parent[1]+rng.uniform(-.025,.025),
                        parent[2]+rng.uniform(-.025,.025),parent[3]+rng.uniform(-1.1,-.55)])
    return {'seed':seed,'config':asdict(cfg),'observers':observers,'targets':targets,
            'initial_sigma_r':rng.uniform(800,1200,n_targets).tolist(),
            'initial_sigma_v':rng.uniform(70,110,n_targets).tolist(),
            'epoch_utc':'2026-09-08T00:00:00Z',
            'frame':'Earth-centered inertial; circular two-body elements',
            'target_source':'synthetic independent elements; no real-object identity'}


class Physics:
    def __init__(self,scene):
        self.scene=scene;self.c=Config(**scene['config'])
        self.times=np.arange(round((self.c.horizon+self.c.hold+65)/self.c.dt)+1)*self.c.dt
        self.N=len(scene['targets']);self.S=len(scene['observers'])
        self.x=np.array([circular(*e,self.times) for e in scene['targets']])
        self.obs=np.array([circular(*e,self.times) for e in scene['observers']])
        self.p0=np.array([np.diag([sr*sr]*3+[sv*sv]*3) for sr,sv in
                          zip(scene['initial_sigma_r'],scene['initial_sigma_v'])])
        nt=len(self.times);dt=self.c.dt
        self.F=np.zeros((self.N,nt-1,6,6));self.F[:,:,:3,:3]=np.eye(3)
        self.F[:,:,3:,3:]=np.eye(3);self.F[:,:,:3,3:]=dt*np.eye(3)
        r=self.x[:,:-1,:3];norm=np.linalg.norm(r,axis=-1)
        g=MU*(3*r[:,:,:,None]*r[:,:,None,:]/norm[:,:,None,None]**5-np.eye(3)/norm[:,:,None,None]**3)
        self.F[:,:,3:,:3]=dt*g
        self.F[:,:,:3,:3]+=.5*dt*dt*g
        self.F[:,:,3:,3:]+=.5*dt*dt*g
        self.Q=self.c.q*np.block([[dt**3/3*np.eye(3),dt**2/2*np.eye(3)],
                                  [dt**2/2*np.eye(3),dt*np.eye(3)]])
        # Fundamental matrices and transformed process-noise integrals, enabling
        # arbitrary coast propagation without resetting the posterior.
        self.phi=np.zeros((self.N,nt,6,6));self.phi[:,0]=np.eye(6)
        for k in range(nt-1):self.phi[:,k+1]=self.F[:,k]@self.phi[:,k]
        self.inv=np.linalg.inv(self.phi)
        increments=self.inv[:,1:]@self.Q@self.inv[:,1:].swapaxes(-1,-2)
        self.W=np.zeros_like(self.phi);self.W[:,1:]=np.cumsum(increments,axis=1)
        self.Gr=self.phi[:,:,:3,:].swapaxes(-1,-2)@self.phi[:,:,:3,:]
        self.Gv=self.phi[:,:,3:,:].swapaxes(-1,-2)@self.phi[:,:,3:,:]
        self.wr=np.sum(self.W*self.Gr,axis=(-1,-2));self.wv=np.sum(self.W*self.Gv,axis=(-1,-2))
        # Geometry evaluated on the actual 0.2 s grid, including Earth occultation.
        self.visible=np.zeros((self.N,self.S,nt),dtype=bool)
        for i in range(self.N):
            d=self.x[i,None,:,:3]-self.obs[:,:,:3];rho=np.linalg.norm(d,axis=-1)
            u=d/rho[:,:,None];o=self.obs[:,:,:3]
            nadir=-o/np.linalg.norm(o,axis=-1)[:,:,None]
            angle=np.degrees(np.arccos(np.clip(np.sum(u*nadir,axis=-1),-1,1)))
            frac=np.clip(-np.sum(o*d,axis=-1)/(rho*rho),0,1)
            clear=np.linalg.norm(o+frac[:,:,None]*d,axis=-1)>=RE
            # Tracking angular rate independent of finite differences.
            dv=self.x[i,None,:,3:]-self.obs[:,:,3:]
            angular=np.linalg.norm(np.cross(d,dv),axis=-1)/(rho*rho)*180/np.pi
            self.visible[i]=(rho<=self.c.max_range)&(angle<=self.c.off_nadir)&clear&(angular<=self.c.omega)
        self.sat_options=[np.flatnonzero(self.visible[i,:,:self.tick(self.c.horizon)+1].any(axis=1)).tolist() for i in range(self.N)]
        self.cache={};self.hits=0;self.misses=0

    def tick(self,t):return int(round(t/self.c.dt))
    def ceil_tick(self,t):return int(math.ceil(t/self.c.dt-1e-8))
    def risk(self,p):return max(math.sqrt(max(0,np.trace(p[:3,:3])))/self.c.epsilon_r,
                                math.sqrt(max(0,np.trace(p[3:,3:])))/self.c.epsilon_v)
    def coast(self,i,p,a,b):
        if a==b:return p.copy()
        assert b>=a
        C=self.inv[i,a]@p@self.inv[i,a].T+self.W[i,b]-self.W[i,a]
        out=self.phi[i,b]@C@self.phi[i,b].T
        return (out+out.T)*.5
    def coast_risks(self,i,p,a,b):
        C=self.inv[i,a]@p@self.inv[i,a].T-self.W[i,a]
        vr=np.einsum('kij,ij->k',self.Gr[i,a:b+1],C)+self.wr[i,a:b+1]
        vv=np.einsum('kij,ij->k',self.Gv[i,a:b+1],C)+self.wv[i,a:b+1]
        return np.maximum(np.sqrt(np.maximum(0,vr))/self.c.epsilon_r,np.sqrt(np.maximum(0,vv))/self.c.epsilon_v)
    def hold_risk(self,i,p,k,endpoint=False):
        end=k+self.tick(self.c.hold)
        if endpoint:return self.risk(self.coast(i,p,k,end))
        return float(self.coast_risks(i,p,k,end).max())
    def H_R(self,i,s,k,noise_scale=1.):
        d=self.x[i,k,:3]-self.obs[s,k,:3];rho=np.linalg.norm(d);u=d/rho
        axis=np.array([1.,0,0]) if abs(u[0])<.8 else np.array([0.,1,0])
        e1=np.cross(u,axis);e1/=np.linalg.norm(e1);e2=np.cross(u,e1)
        H=np.zeros((3,6));H[:,:3]=np.array([e1/rho,e2/rho,u])
        R=np.diag([self.c.sigma_angle**2,self.c.sigma_angle**2,self.c.sigma_range**2])*noise_scale**2
        return H,R
    def update(self,i,s,p,k):
        # Cartesian equivalent of two tangent-plane angles plus range. Orthogonal
        # whitening followed by Joseph scalar updates preserves covariance PSD.
        H,R=self.H_R(i,s,k)
        for h,var in zip(H,np.diag(R)):
            ph=p@h;den=float(h@ph+var);K=ph/den
            A=np.eye(6)-np.outer(K,h)
            p=A@p@A.T+var*np.outer(K,K)
        return (p+p.T)*.5
    def step(self,i,p,k):return self.F[i,k-1]@p@self.F[i,k-1].T+self.Q
    def observe(self,i,s,p,a,b,drop_every=0):
        for k in range(a+1,b+1):
            p=self.step(i,p,k)
            if not drop_every or k%drop_every:p=self.update(i,s,p,k)
        return p
    def dwell(self,i,s,p,prior_tick,start,mode='hold',limit=None,drop_every=0):
        limit=self.tick(self.c.continuous_limit) if limit is None else limit
        key=(i,s,prior_tick,start,mode,limit,drop_every,p.tobytes())
        if key in self.cache:self.hits+=1;return self.cache[key]
        self.misses+=1
        current=self.coast(i,p,prior_tick,start)
        stop=None;result=[];first=None
        if mode in ('fixed','fixed_hold'):stop=self.tick(self.c.fixed_dwell)
        for length in range(1,limit+1):
            k=start+length
            if not self.visible[i,s,k]:break
            current=self.step(i,current,k)
            if not drop_every or k%drop_every:current=self.update(i,s,current,k)
            if mode in ('fixed','fixed_hold'):
                passed=length==stop
                if passed and mode=='fixed_hold':passed=self.hold_risk(i,current,k)<=1
            elif mode=='instant':passed=self.risk(current)<=1
            else:
                passed=self.hold_risk(i,current,k,True)<=1
                if passed and mode!='endpoint':passed=self.hold_risk(i,current,k)<=1
            if first is None and passed:
                first=length
                result.append((length,current.copy(),self.hold_risk(i,current,k)))
                if mode in ('fixed','fixed_hold','instant'):break
                stop=min(limit,int(math.ceil(first*1.5)))
                if stop==first:break
            elif first is not None and length==stop:
                hold=self.hold_risk(i,current,k)
                if hold<=1:result.append((length,current.copy(),hold))
                break
            if mode in ('fixed','fixed_hold') and length>=stop:break
        self.cache[key]=result
        return result
    def next_release(self,i,p,k):
        end=self.tick(self.c.horizon)
        if k>=end:return end+1
        risks=self.coast_risks(i,p,k,end)
        hit=np.flatnonzero(risks>=self.c.warning)
        # Always advance by at least one measurement interval after a visit.
        return k+max(1,int(hit[0])) if len(hit) else end+1
