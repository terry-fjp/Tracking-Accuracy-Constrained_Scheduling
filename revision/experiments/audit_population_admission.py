"""Read saved nested-workload plans to distinguish admission from actual service."""
from pathlib import Path
import json,csv,statistics
R=Path(__file__).resolve().parent
rows=[]
for n in [24,48,72,96]:
 for method in ['IP-AD','IP-F20-H']:
  for sid in range(2101,2109):
   for seed in range(3101,3106):
    name=f'population-{n}-{sid}-s{seed}-{method}'
    plan=json.loads((R/'results/plans'/f'{name}.json').read_text())
    run=json.loads((R/'results/runs'/f'{name}.json').read_text())
    admitted={i for i,x in enumerate(plan['genes'][:n]) if x>=.5}
    served={v['target'] for v in plan['plan']}
    assert served<=admitted and n-len(served)==run['unobserved_targets']
    rows.append(dict(population=n,method=method,scene=sid,seed=seed,admitted=len(admitted),excluded_by_admission=n-len(admitted),admitted_without_observation=len(admitted-served),observed=len(served),unobserved=n-len(served)))
summary=[]
for n in [24,48,72,96]:
 for method in ['IP-AD','IP-F20-H']:
  rs=[x for x in rows if x['population']==n and x['method']==method]
  entry=dict(population=n,method=method)
  for metric in ['admitted','excluded_by_admission','admitted_without_observation','observed','unobserved']:
   vals=[statistics.mean(x[metric] for x in rs if x['seed']==seed) for seed in range(3101,3106)]
   entry[metric]=dict(mean=statistics.mean(vals),sample_sd=statistics.stdev(vals))
  summary.append(entry)
with (R/'results/population-admission-audit.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(R/'results/population-admission-summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps([x for x in summary if x['population']==96],indent=2))
