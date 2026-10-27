import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from run_revision import *
from j2_revision import Geometry,visible,mask_calendar
from space_cluster_scheduling.model import Config
from types import SimpleNamespace
import inspect
# Isolated dwell fixture: first feasible at 3 ticks, no geometry cutoff.
class Fixture:
 dwell=Physics.dwell
 def __init__(self,m,limit=10):
  self.c=Config(duration_multipliers=(1,m));self.cache={};self.hits=0;self.misses=0;self.visible=np.ones((1,1,20),bool)
 def coast(self,*args):return np.zeros((1,1))
 def step(self,i,p,k):return np.array([[k]])
 def update(self,i,s,p,k):return p
 def hold_risk(self,i,p,k,endpoint=False):return .9 if k>=3 else 2
for m,expected in [(1.25,[3,4]),(1.5,[3,5]),(1.75,[3,6])]:
 p=Fixture(m);a=p.dwell(0,0,np.zeros((1,1)),0,0,limit=10);assert [x[0] for x in a]==expected
p=Fixture(1.75);assert [x[0] for x in p.dwell(0,0,np.zeros((1,1)),0,0,limit=4)]==[3,4]
p=Fixture(1.75);p.visible[:,:,5:]=False;assert [x[0] for x in p.dwell(0,0,np.zeros((1,1)),0,0,limit=10)]==[3]
# A direct fixture for the lexicographic priority verifies all three config weights
# affect the same expression used by the submitted decoder, without touching constraints.
from space_cluster_scheduling import scheduler
source=inspect.getsource(scheduler.decode)
assert all(f'p.c.ranking_weights[{j}]' in source for j in range(3))
features=np.array([.2,.8,.4]);weights=PROTOCOL['settings'][8]['overrides']['ranking_weights'];assert abs(np.dot(features,weights)-np.dot(features,[.5,.3,.2]))>1e-10
# Check the execution-mask policy against original, unperturbed geometry.
checks=[]
for sid in [2101,2106]:
 p=Geometry(scene_for(sid,PROTOCOL['settings'][1]));p.x=p.two_x;p.obs=p.two_obs;p.visible=visible(p.x,p.obs,p.c)
 for seed in [3101,3105]:
  for method in ['IP-AD','IP-F20-H']:
   plan=json.loads((ARCHIVE/'plans'/f'P1-{sid}-s{seed}-b48-{method}.json').read_text())['plan'];keep,drop=mask_calendar(p,plan)
   assert not drop,(sid,seed,method,drop)
   checks.append([sid,seed,method,len(keep)])
write(ROOT/'checks/interface-semantics.json',dict(status='PASS',multiplier_candidate_ticks={'1.25':[3,4],'1.5':[3,5],'1.75':[3,6]},cap_and_visibility=True,weight_interface=True,unperturbed_mask_retains_all=checks))
print('INTERFACE AND MASK CHECKS PASS')
