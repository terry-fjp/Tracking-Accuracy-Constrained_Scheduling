"""Three-key search with causal append decoding and independent full replay."""
from __future__ import annotations
from dataclasses import dataclass,asdict
import math,time
import numpy as np
from .model import Physics

@dataclass
class Visit:
    target:int
    satellite:int
    setup:int
    capture:int
    start:int
    end:int
    hold_risk:float
    shared:bool=False


def slew(deg,c):
    return 2*math.sqrt(deg/c.acceleration) if deg<=c.omega**2/c.acceleration else deg/c.omega+c.omega/c.acceleration

def pointing(p,i,s,k):
    d=p.x[i,k,:3]-p.obs[s,k,:3]
    return d/np.linalg.norm(d)

def same_field(p,i,j,s,a,b):
    # Use one common candidate boresight (normalized vector sum). For two members,
    # its tangent coordinates enclose both in the actual rectangular 1 x 1 field.
    if i==j:return True
    if not (p.visible[i,s,a:b+1].all() and p.visible[j,s,a:b+1].all()):return False
    u=p.x[i,a:b+1,:3]-p.obs[s,a:b+1,:3];u/=np.linalg.norm(u,axis=1)[:,None]
    v=p.x[j,a:b+1,:3]-p.obs[s,a:b+1,:3];v/=np.linalg.norm(v,axis=1)[:,None]
    bore=u+v;bore/=np.linalg.norm(bore,axis=1)[:,None]
    # Conservative circular inner field is sufficient for the rectangular field.
    coshalf=math.cos(math.radians(p.c.fov/2))
    return bool(np.all(np.sum(u*bore,axis=1)>=coshalf) and np.all(np.sum(v*bore,axis=1)>=coshalf))

def setup_start(p,i,s,release,prev,grouping):
    a=max(release,prev.end if prev else 0)
    visible=p.visible[i,s]
    end=p.tick(p.c.horizon)
    while a<end:
        good=np.flatnonzero(visible[a:end+1])
        if not len(good):return None
        a+=int(good[0])
        shared=False
        if prev:
            old=pointing(p,prev.target,s,prev.end)
            for _ in range(20):
                if a>=end:return None
                new=pointing(p,i,s,a)
                angle=math.degrees(math.acos(float(np.clip(old@new,-1,1))))
                motion=0 if i==prev.target else slew(angle,p.c)+p.c.settle
                required=prev.end+p.ceil_tick(motion)
                if required<=a:break
                a=required
            else:return None
            if a>=end:return None
            shared=grouping!='none' and prev.target!=i and a-prev.end<=p.tick(10) and same_field(p,i,prev.target,s,prev.end,a)
            if grouping=='static':shared=shared and same_field(p,i,prev.target,s,0,0)
        capture=0 if shared or (prev and prev.target==i and a==prev.end) else p.tick(p.c.capture)
        start=a+capture
        if start<end and visible[a:start+1].all():return a,start,shared
        a=start+1
    return None


def decode(p:Physics,genes,mode='hold',opportunity='dynamic',grouping='dynamic',prefix=(),now=0,
           fault=None,drop_every=0,base_plan=(),deadline_priority=False):
    n=p.N;end=p.tick(p.c.horizon)
    accepted=genes[:n]>=.5;order=genes[n:2*n];satkey=genes[2*n:]
    plans=list(prefix);last=[None]*p.S;laser=np.zeros(p.S,dtype=int)
    P=p.p0.copy();pt=np.zeros(n,dtype=int);release=np.full(n,now,dtype=int);visits=np.zeros(n,dtype=int)
    # Reconstruct inherited posteriors from actually executed/frozen observations.
    for v in sorted(prefix,key=lambda a:a.end):
        P[v.target]=p.coast(v.target,P[v.target],int(pt[v.target]),v.start)
        P[v.target]=p.observe(v.target,v.satellite,P[v.target],v.start,v.end,drop_every)
        pt[v.target]=v.end;visits[v.target]+=1;accepted[v.target]=True
        laser[v.satellite]+=v.end-v.start
        if last[v.satellite] is None or v.end>last[v.satellite].end:last[v.satellite]=v
    for i in range(n):
        if visits[i]:release[i]=max(now,p.next_release(i,P[i],int(pt[i])))
    active={i for i in range(n) if accepted[i] and release[i]<end}
    static={};counts={'candidate_requests':0,'candidate_actions':0,'shared_visits':0,'no_candidate':0}
    # Across one decode the target's actual posterior changes only on its visit;
    # other satellites need refreshing only when their own append calendar changes.
    action_cache={}
    while active:
        options={}
        for i in sorted(active):
            actions=[];sats=p.sat_options[i]
            preferred=sats[min(len(sats)-1,int(satkey[i]*len(sats)))] if sats else -1
            for s in sats:
                prev=last[s]
                if fault is not None and s==fault[0]:
                    # Candidate may finish before fault. Entire affected satellite
                    # unavailable from the specified grid tick onward.
                    cap_end=min(end,fault[1])
                else:cap_end=end
                key=(i,s,int(pt[i]),int(release[i]),prev.end if prev else -1,
                     prev.target if prev else -1,int(laser[s]))
                if key in action_cache:
                    actions.extend(action_cache[key]);continue
                found=[];counts['candidate_requests']+=1
                setup=setup_start(p,i,s,int(release[i]),prev,grouping)
                if setup is not None:
                    a,start,shared=setup
                    limit=min(p.tick(p.c.continuous_limit),p.tick(p.c.laser_budget)-int(laser[s]),cap_end-start)
                    if limit>0 and start>=pt[i]:
                        candidates=p.dwell(i,s,P[i],int(pt[i]),start,mode,limit,drop_every)
                        for length,post,hold in candidates:
                            finish=start+length
                            # Shared acquisition must remain supported throughout
                            # the complete pair's current execution interval.
                            if shared and not same_field(p,i,prev.target,s,prev.end,finish):continue
                            if finish>cap_end:continue
                            if not p.visible[i,s,a:finish+1].all():continue
                            gain=max(0.,p.risk(p.coast(i,P[i],int(pt[i]),finish))-p.risk(post))
                            density=gain/max(1,finish-(prev.end if prev else a))
                            pref=0 if s==preferred else .1
                            score=density-pref-.0001*(finish-start)
                            v=Visit(i,s,prev.end if prev and a==prev.end else max(0,a-p.tick(p.c.settle)),a,start,finish,hold,shared)
                            # Occupancy begins when previous attitude leaves its
                            # old direction; waiting before that is not charged.
                            if prev:
                                angle=math.degrees(math.acos(float(np.clip(pointing(p,prev.target,s,prev.end)@pointing(p,i,s,a),-1,1))))
                                v.setup=max(prev.end,a-p.ceil_tick(slew(angle,p.c)+p.c.settle))
                            else:v.setup=a
                            found.append((score,v,post))
                action_cache[key]=found;actions.extend(found)
            if actions:
                # Deduplicate physical (satellite, actual measurement interval).
                actions=list({(v.satellite,v.start,v.end):(sc,v,post) for sc,v,post in actions}.values())
                # Prefer on-time initial observations before comparing information
                # density. Retain late observations to improve actual estimation
                # accuracy for every target, including those missing the deadline.
                if deadline_priority and visits[i]==0:
                    ontime_actions=[z for z in actions if z[1].end<=p.tick(p.c.init_deadline)]
                    if ontime_actions:actions=ontime_actions
                actions.sort(key=lambda z:(-z[0],z[1].end,z[1].satellite))
                options[i]=actions;counts['candidate_actions']+=len(actions)
            else:counts['no_candidate']+=1
        if not options:break
        for i,acts in options.items():
            static.setdefault(i,(len(acts)==1,acts[0][0]-acts[1][0] if len(acts)>1 else 0))
        def priority(i):
            acts=options[i]
            single,regret=(len(acts)==1,acts[0][0]-acts[1][0] if len(acts)>1 else 0)
            if opportunity=='static':single,regret=static[i]
            if opportunity=='none':single,regret=False,0.
            # Temporal blocks avoid selecting a distant future demand at the
            # expense of an imminent initialization; explicit discrete heuristic.
            earliest=acts[0][1].start
            block=earliest//p.tick(30)
            urgent=p.risk(p.coast(i,P[i],int(pt[i]),earliest))
            ontime=int(deadline_priority and visits[i]==0 and acts[0][1].end<=p.tick(p.c.init_deadline))
            return (ontime,-block,float(single),p.c.ranking_weights[0]*order[i]+p.c.ranking_weights[1]*min(urgent,3)/3+p.c.ranking_weights[2]*min(regret,1),-i)
        i=max(options,key=priority);_,v,post=options[i][0]
        plans.append(v);P[i]=post;pt[i]=v.end;visits[i]+=1;last[v.satellite]=v
        laser[v.satellite]+=v.end-v.start;counts['shared_visits']+=int(v.shared)
        release[i]=p.next_release(i,post,v.end)
        active.remove(i)
        if release[i]<end and v.end>v.start:active.add(i)
        # Finite horizon plus strictly positive dwell bounds the loop.
        if len(plans)>n*100:raise RuntimeError('nonterminating revisit construction')
    metrics=replay(p,plans,drop_every=drop_every,detail=False,base_plan=base_plan,now=now)
    if metrics['resource_violations']:raise AssertionError('resource audit failed')
    if mode in ('hold','endpoint') and metrics['hold_violations']:raise AssertionError('hold audit failed')
    return {'plan':plans,'metrics':metrics,'objective':tuple(metrics['objective']),'decoder':counts}


def replay(p,plan,drop_every=0,detail=False,base_plan=(),now=0):
    """Fresh evaluation from original p0, without trusting decoder hold labels."""
    end=p.tick(p.c.horizon);init=p.tick(p.c.init_deadline)
    rates=np.zeros((p.N,end+1));initialized=np.zeros(p.N,dtype=bool)
    hold_violations=0;max_label_error=0.;resource_violations=0
    actual_holds=[]
    bytarget=[[v for v in plan if v.target==i] for i in range(p.N)]
    for i,vs in enumerate(bytarget):
        current=p.p0[i].copy();k0=0
        for v in sorted(vs,key=lambda z:z.start):
            if v.start<k0:raise ValueError('concurrent measurements for one target')
            coast=p.coast_risks(i,current,k0,v.start);rates[i,k0:v.start+1]=coast
            current=p.coast(i,current,k0,v.start)
            for k in range(v.start+1,v.end+1):
                current=p.step(i,current,k)
                if not drop_every or k%drop_every:current=p.update(i,v.satellite,current,k)
                rates[i,k]=p.risk(current)
            hold=p.hold_risk(i,current,v.end)
            actual_holds.append(hold);hold_violations+=int(hold>1+1e-8)
            max_label_error=max(max_label_error,abs(hold-v.hold_risk))
            if v.end<=init and p.risk(current)<=1:initialized[i]=True
            k0=v.end
        rates[i,k0:]=p.coast_risks(i,current,k0,end)
    laser=0;occupied=0
    for s in range(p.S):
        cal=sorted([v for v in plan if v.satellite==s],key=lambda v:v.start)
        use=sum(v.end-v.start for v in cal);laser+=use
        resource_violations+=int(use>p.tick(p.c.laser_budget))
        for j,v in enumerate(cal):
            resource_violations+=int(not p.visible[v.target,s,v.capture:v.end+1].all())
            resource_violations+=int(v.end-v.start>p.tick(p.c.continuous_limit))
            resource_violations+=int(v.start-v.capture < (0 if v.shared else p.tick(p.c.capture)))
            if j:
                prev=cal[j-1]
                angle=math.degrees(math.acos(float(np.clip(pointing(p,prev.target,s,prev.end)@pointing(p,v.target,s,v.capture),-1,1))))
                required=0 if prev.target==v.target else p.ceil_tick(slew(angle,p.c)+p.c.settle)
                resource_violations+=int(v.capture<prev.end+required)
            if v.shared and j:resource_violations+=int(not same_field(p,v.target,cal[j-1].target,s,cal[j-1].end,v.end))
            occupied+=v.end-v.setup
    valid=initialized & (rates[:,init:].max(axis=1)<=1+1e-8)
    bad=rates[:,init:]>1+1e-8
    # Discrete right-endpoint interval accounting, not count*dt including both ends.
    exposure=bad[:,1:].sum(axis=1)*p.c.dt
    longest=[];recoveries=[]
    for row,point_states in zip(bad[:,1:],bad):
        intervals=np.diff(np.r_[False,row,False].astype(int));starts=np.flatnonzero(intervals==1);ends=np.flatnonzero(intervals==-1)
        longest.append(float(max(ends-starts,default=0)*p.c.dt))
        recoveries.append(int(np.count_nonzero(np.diff(point_states.astype(int))==-1)))
    current_by={(v.target,j):v for i in range(p.N) for j,v in enumerate(sorted(bytarget[i],key=lambda z:z.start)) if v.target==i and v.start>=now}
    old_by={}
    for i in range(p.N):
        for j,v in enumerate(sorted([v for v in base_plan if v.target==i],key=lambda z:z.start)):
            if v.start>=now:old_by[(i,j)]=v
    cancels=len(set(old_by)-set(current_by));reassign=0;shift=0.
    for key in set(old_by)&set(current_by):
        reassign+=int(old_by[key].satellite!=current_by[key].satellite)
        shift+=abs(old_by[key].start-current_by[key].start)*p.c.dt
    maxrisk=float(rates[:,init:].max());avgrisk=float(rates[:,init:].mean())
    J4=cancels+reassign+shift/max(p.c.horizon,1)
    objective=[-int(valid.sum()),maxrisk,avgrisk,occupied*p.c.dt,laser*p.c.dt,J4]
    out={'persistent_targets':int(valid.sum()),'target_count':p.N,'persistent_rate':float(valid.mean()),
         'initialized_targets':int(initialized.sum()),'max_risk':maxrisk,'mean_risk':avgrisk,
         'laser_s':laser*p.c.dt,'occupied_s':occupied*p.c.dt,'visits':len(plan),
         'mean_exceedance_s':float(exposure.mean()),'max_exceedance_s':float(max(longest,default=0)),
         'recovery_count':sum(recoveries),'hold_violations':hold_violations,
         'max_label_error':max_label_error,'resource_violations':resource_violations,
         'cancellations':cancels,'reassignments':reassign,'shift_s':shift,'objective':objective,
         'shared_visits':sum(v.shared for v in plan)}
    if detail:out.update(risk_trajectories=rates.tolist(),valid_targets=valid.tolist(),
                         actual_hold_risks=actual_holds,exceedance_s=exposure.tolist())
    return out


def optimize(p,method='H-SADE',seed=1,budget=48,seconds=None,**decode_args):
    """Same three genes, population, guided initialization, decoder and evaluator.

    H-SADE differs solely through self-adaptive differential strategies. GA and
    LNS-ref are reference implementations; no claim to reproduce full ALNS.
    """
    start=time.perf_counter();rng=np.random.default_rng(seed);n=p.N;d=3*n
    size=min(12,budget);X=rng.random((size,d));X[:,:n]=.5+.5*X[:,:n]
    X[0,:n]=1;X[0,n:2*n]=np.linspace(1,0,n);X[0,2*n:]=.5
    if 'warm' in decode_args:
        warm=decode_args.pop('warm')
        if warm is not None:X[0]=warm
    p.cache={};p.hits=p.misses=0 # same cold-cache policy per run
    evaluations=[];history=[];count=0
    def evaluate(x):
        nonlocal count
        ans=decode(p,x,**decode_args);count+=1
        ans['genes']=x.copy();return ans
    def expired():return seconds is not None and time.perf_counter()-start>=seconds and count>0
    for j in range(size):
        if expired():break
        evaluations.append(evaluate(X[j]))
    X=X[:len(evaluations)];size=len(X)
    best=min(evaluations,key=lambda z:z['objective'])
    history.append([count,time.perf_counter()-start,best['metrics']['persistent_targets']])
    Fmem=.55;CRmem=.75;succ=np.ones(3);trials=np.ones(3)
    while count<budget and not expired():
        j=(count-size)%size
        if method=='Random':
            candidate=rng.random(d);candidate[:n]=.5+.5*candidate[:n];strategy=None
        elif method=='GA':
            def tournament():
                a,b=rng.integers(size,size=2);return a if evaluations[a]['objective']<evaluations[b]['objective'] else b
            a,b=tournament(),tournament();candidate=np.where(rng.random(d)<.5,X[a],X[b])
            mutation=rng.random(d)<1.5/d;candidate[mutation]=rng.random(mutation.sum())
            strategy=None
        elif method=='LNS-ref':
            candidate=best['genes'].copy();fraction=(.05,.1,.2)[count%3]
            ix=rng.choice(d,max(1,int(fraction*d)),replace=False);candidate[ix]=rng.random(len(ix));strategy=None
        else:
            strategy=0 if method=='DE' else int(rng.choice(3,p=(succ/trials)/(succ/trials).sum()))
            F=.55 if method=='DE' else float(np.clip(rng.normal(Fmem,.12),.15,.95))
            CR=.75 if method=='DE' else float(np.clip(rng.normal(CRmem,.1),.1,.98))
            pool=[k for k in range(size) if k!=j]
            if len(pool)<4:pool=list(range(size))
            ix=rng.choice(pool,4,replace=len(pool)<4)
            if strategy==0:mutant=X[ix[0]]+F*(X[ix[1]]-X[ix[2]])
            elif strategy==1:mutant=X[j]+F*(best['genes']-X[j])+F*(X[ix[0]]-X[ix[1]])
            else:mutant=best['genes']+F*(X[ix[0]]-X[ix[1]]+X[ix[2]]-X[ix[3]])
            mask=rng.random(d)<CR;mask[rng.integers(d)]=True
            candidate=np.clip(np.where(mask,mutant,X[j]),0,1)
        result=evaluate(candidate)
        if strategy is not None:trials[strategy]+=1
        if result['objective']<evaluations[j]['objective']:
            X[j]=candidate;evaluations[j]=result
            if strategy is not None:
                succ[strategy]+=1;Fmem=.9*Fmem+.1*F;CRmem=.9*CRmem+.1*CR
        if result['objective']<best['objective']:best=result
        history.append([count,time.perf_counter()-start,best['metrics']['persistent_targets']])
    best['search']={'method':method,'seed':seed,'evaluations':count,'elapsed_s':time.perf_counter()-start,
                    'cache_hits':p.hits,'cache_misses':p.misses,'history':history,
                    'time_budget_s':seconds,'budget':budget,
                    'initialization_evaluations':size,'search_evaluations':count-size,
                    'initial_best_targets':history[0][2]}
    return best
