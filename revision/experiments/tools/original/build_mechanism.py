"""Recompute the dwell diagnostic at the current default forecast interval."""
import json
import numpy as np
from space_cluster_scheduling.model import Physics, make_scene
from space_cluster_scheduling.paths import DATA
OUT=DATA/'diagnostics'
def write_json(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2)+'\n')

def mechanisms():
    rows=[];curves=[]
    for seed in [4101,4102,4103]:
        for scale in [.5,1.,2.]:
            scene=make_scene(seed,2)
            scene['initial_sigma_r']=[v*scale for v in scene['initial_sigma_r']]
            scene['initial_sigma_v']=[v*scale for v in scene['initial_sigma_v']]
            p=Physics(scene);i=0;s=p.sat_options[i][0]
            visible=p.visible[i,s];start=max(p.tick(5),int(np.flatnonzero(visible)[0])+p.tick(5))
            current=p.coast(i,p.p0[i],0,start);trace=[]
            first_full=first_end=first_instant=None
            for L in range(1,p.tick(45)+1):
                k=start+L
                if not visible[k]:break
                current=p.update(i,s,p.step(i,current,k),k)
                hold=p.hold_risk(i,current,k);endpoint=p.hold_risk(i,current,k,True);instant=p.risk(current)
                if first_full is None and hold<=1:first_full=L*p.c.dt
                if first_end is None and endpoint<=1:first_end=L*p.c.dt
                if first_instant is None and instant<=1:first_instant=L*p.c.dt
                trace.append([L*p.c.dt,instant,endpoint,hold])
            rows.append({'scenario':seed,'prior_sigma_multiplier':scale,'satellite':s,'start_s':start*p.c.dt,
                         'range_km':float(np.linalg.norm(p.x[i,start,:3]-p.obs[s,start,:3]))/1000,
                         'instant_dwell_s':first_instant,'endpoint_dwell_s':first_end,'hold_dwell_s':first_full,
                         'max_envelope_endpoint_difference':max(abs(x[2]-x[3]) for x in trace)})
            curves.append({'scenario':seed,'scale':scale,'hold_s':p.c.hold,'curve':trace})
    write_json(OUT/'mechanism.json',rows);write_json(OUT/'mechanism-curves.json',curves)

if __name__=='__main__':
    mechanisms()
