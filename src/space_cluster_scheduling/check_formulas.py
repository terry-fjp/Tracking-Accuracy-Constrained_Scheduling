"""Independent numerical checks of equations; no claim of flight validation."""
import json, math
from pathlib import Path
import numpy as np
from .model import Physics, MU

from .paths import ROOT, DATA, REFERENCE
from .scenes import load_scene
OUT=ROOT/'results'

def expm_series(a):
    # Scaling and Taylor series, independent of the propagation implementation.
    scale=max(0,int(np.ceil(np.log2(max(np.linalg.norm(a,np.inf),1)))))
    b=a/(2**scale);term=np.eye(len(a));out=term.copy()
    for n in range(1,60):
        term=term@b/n;out+=term
        if np.linalg.norm(term,np.inf)<1e-18:break
    for _ in range(scale):out=out@out
    return out

def main():
    p=Physics(load_scene(DATA/'scenes/2101.json'))
    i=2;s=p.sat_options[i][0];k=100;dt=p.c.dt
    r=p.x[i,k,:3];g=MU*(3*np.outer(r,r)/np.linalg.norm(r)**5-np.eye(3)/np.linalg.norm(r)**3)
    accel=lambda x:-MU*x/np.linalg.norm(x)**3
    eps=1.;numeric=np.column_stack([(accel(r+eps*np.eye(3)[j])-accel(r-eps*np.eye(3)[j]))/(2*eps) for j in range(3)])
    gravity_error=np.linalg.norm(numeric-g)/np.linalg.norm(g)
    C,R=p.H_R(i,s,k);rho=np.linalg.norm(r-p.obs[s,k,:3]);basis=C[:,:3].copy();basis[:2]*=rho
    def h(x):
        d=x-p.obs[s,k,:3];z=basis@d
        return np.array([np.arctan2(z[0],z[2]),np.arctan2(z[1],z[2]),np.linalg.norm(d)])
    hh=np.column_stack([(h(r+eps*np.eye(3)[j])-h(r-eps*np.eye(3)[j]))/(2*eps) for j in range(3)])
    measurement_errors=[float(np.linalg.norm(hh[j]-C[j,:3])/np.linalg.norm(C[j,:3])) for j in range(3)]
    prior=p.step(i,p.p0[i],k);K=np.linalg.solve(C@prior@C.T+R,C@prior).T
    B=np.eye(6)-K@C;batch=B@prior@B.T+K@R@K.T;sequential=p.update(i,s,prior,k)
    joseph_error=np.linalg.norm(batch-sequential)/np.linalg.norm(batch)
    A=np.block([[np.zeros((3,3)),np.eye(3)],[g,np.zeros((3,3))]])
    LQL=np.zeros((6,6));LQL[3:,3:]=p.c.q*np.eye(3)
    aug=np.block([[A,LQL],[np.zeros((6,6)),-A.T]])
    E=expm_series(aug*dt);Fexact=E[:6,:6];Qexact=E[:6,6:]@Fexact.T
    F2=np.eye(6)+A*dt+.5*A@A*dt**2
    ferr=np.linalg.norm(F2-Fexact)/np.linalg.norm(Fexact)
    qerr=np.linalg.norm(p.Q-Qexact)/np.linalg.norm(Qexact)
    c=prior.copy()
    span=p.tick(p.c.hold)
    for tick in range(k+1,k+span+1):c=p.step(i,c,tick)
    fast=p.coast(i,prior,k,k+span)
    coast_error=np.linalg.norm(c-fast)/np.linalg.norm(c)
    rng=np.random.default_rng(921);min_eig=math.inf
    for _ in range(100):
        a=rng.normal(size=(6,6));prior=a@a.T+np.eye(6)*1e-6
        post=p.update(i,s,prior,k);min_eig=min(min_eig,float(np.linalg.eigvalsh(post).min()))
    checks={'gravity_jacobian_relative_error':float(gravity_error),'angle_range_jacobian_relative_errors':measurement_errors,
        'batch_vs_sequential_joseph_relative_error':float(joseph_error),
        'F_second_order_vs_matrix_exponential_relative_error':float(ferr),
        'Q_constant_velocity_vs_frozen_gravity_exact_relative_error':float(qerr),
        'hold_s':p.c.hold,'hold_fast_vs_stepwise_coast_relative_error':float(coast_error),
        'minimum_eigenvalue_100_random_Joseph_updates':min_eig}
    assert gravity_error<1e-7 and max(measurement_errors)<1e-6
    assert joseph_error<1e-10 and ferr<1e-7 and qerr<1e-6 and coast_error<1e-7 and min_eig>0
    checks['status']='PASS';checks['interpretation']='local numerical consistency at specified states and dt; Q is an approximation for gravity dynamics, not its exact discretization'
    OUT.mkdir(exist_ok=True);(OUT/'formula-checks.json').write_text(json.dumps(checks,indent=2))
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
