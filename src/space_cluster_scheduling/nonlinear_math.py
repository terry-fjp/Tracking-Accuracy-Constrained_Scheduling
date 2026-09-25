"""Nonlinear two-body truth and angle-range model, preserved from 20260909.
Functions copied without mathematical changes; see provenance in README.
"""
import numpy as np
from .model import MU

def rhs(x):
    r=x[...,:3];a=-MU*r/np.linalg.norm(r,axis=-1)[...,None]**3
    return np.concatenate((x[...,3:],a),axis=-1)

def rk4(x,dt):
    a=rhs(x);b=rhs(x+dt*a/2);c=rhs(x+dt*b/2);d=rhs(x+dt*c)
    return x+dt*(a+2*b+2*c+d)/6

def transition(x,dt):
    r=x[...,:3];rho=np.linalg.norm(r,axis=-1)
    g=MU*(3*r[..., :,None]*r[...,None,:]/rho[...,None,None]**5-np.eye(3)/rho[...,None,None]**3)
    f=np.broadcast_to(np.eye(6),x.shape[:-1]+(6,6)).copy()
    f[...,:3,3:]=dt*np.eye(3);f[...,3:,:3]=dt*g
    f[...,:3,:3]+=.5*dt**2*g;f[...,3:,3:]+=.5*dt**2*g
    return f

def measurement(x,obs,basis):
    d=x[...,:3]-obs[:3];z=d@basis.T;rho=np.linalg.norm(d,axis=-1)
    y=np.stack((np.arctan2(z[...,0],z[...,2]),np.arctan2(z[...,1],z[...,2]),rho),axis=-1)
    h=np.zeros(x.shape[:-1]+(3,6))
    for j in range(2):h[...,j,:3]=(z[...,2,None]*basis[j]-z[...,j,None]*basis[2])/(z[...,2]**2+z[...,j]**2)[...,None]
    h[...,2,:3]=d/rho[...,None]
    return y,h
