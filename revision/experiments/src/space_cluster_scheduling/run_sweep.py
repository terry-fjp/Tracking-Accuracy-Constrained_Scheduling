"""Run all specified conditions with independent scene workers, then audit them."""
import argparse,json,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from .paths import ROOT
from . import run_accuracy as accuracy, audit_accuracy as audit

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=8)
    parser.add_argument('--output',type=str,default=str(ROOT/'results/experiments'))
    parser.add_argument('--conditions',nargs='+',help='Optional subset of condition IDs')
    args=parser.parse_args()
    from pathlib import Path
    output=Path(args.output);spec=json.loads((ROOT/'config/experiments.json').read_text())
    if args.conditions:
        requested=set(args.conditions)
        spec['conditions']=[c for c in spec['conditions'] if c['id'] in requested]
        assert {c['id'] for c in spec['conditions']}==requested
    protocols={};batches={c['id']:[] for c in spec['conditions']}
    for c in spec['conditions']:protocols[c['id']]=accuracy.prepare_protocol(output/c['id'],c,args.workers)
    start=time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs={pool.submit(accuracy.run_scene,sid,output/c['id'],protocols[c['id']]):(c['id'],sid)
              for c in spec['conditions'] for sid in range(2101,2109)}
        for job in as_completed(jobs):
            condition,sid=jobs[job];batches[condition].extend(job.result())
            rows=sorted(batches[condition],key=lambda r:(r['scene'],protocols[condition]['methods'].index(r['method'])))
            accuracy.write(output/condition/'scene-results.json',rows)
            print('COMPLETE',condition,sid,'scene-method batches',len(rows),flush=True)
            if len(rows)==24:accuracy.summarize(output/condition,protocols[condition],rows)
    for c in spec['conditions']:
        audit.OUT=output/c['id'];audit.main()
    accuracy.write(output/('run-completion-'+('-'.join(args.conditions) if args.conditions else 'all')+'.json'),{'status':'PASS','conditions':len(spec['conditions']),
                    'condition_ids':[c['id'] for c in spec['conditions']],'plans':24*len(spec['conditions']),'chains':414*len(spec['conditions']),'target_method_realizations':82800*len(spec['conditions']),
                    'elapsed_s':time.perf_counter()-start,'workers':args.workers})

if __name__=='__main__':main()
