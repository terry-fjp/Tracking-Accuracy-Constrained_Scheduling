from run_revision import *
import importlib.util
spec=importlib.util.spec_from_file_location('mechanism',ROOT/'tools/original/build_mechanism.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
module.OUT=OUT/'figure3-check';module.mechanisms()
old=json.loads((ROOT/'data/diagnostics/mechanism-curves.json').read_text());new=json.loads((module.OUT/'mechanism-curves.json').read_text())
assert old==new
rows=json.loads((module.OUT/'mechanism.json').read_text());inputs=[]
for seed in [4101,4102,4103]:
 sc=make_scene(seed,2);row=next(r for r in rows if r['scenario']==seed and r['prior_sigma_multiplier']==1)
 inputs.append(dict(seed=seed,target=0,observer=row['satellite'],start_s=row['start_s'],range_km=row['range_km'],initial_sigma_r_m=sc['initial_sigma_r'][0],initial_sigma_v_m_s=sc['initial_sigma_v'][0]))
assert [r['observer'] for r in inputs]==[8,23,26] and all(r['start_s']==5 for r in inputs)
write(ROOT/'checks/figure3-reproduction.json',dict(status='PASS',curves_exact=True,inputs=inputs,scope='Recomputed original diagnostic with default current code; not a new main-experiment scene. Original manuscript figure remains unchanged.'))
print('FIGURE 3 NUMERICAL REPRODUCTION PASS')
