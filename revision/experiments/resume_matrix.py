"""Resume frozen matrix, skipping complete jobs before costly geometry setup."""
from run_revision import *
if __name__ == '__main__':
 jobs=[]
 for s in PROTOCOL['scenes']:
  for k in PROTOCOL['search_seeds']:
   for c in PROTOCOL['settings']:
    files=[OUT/'runs'/f"{c['id']}-{s}-s{k}-{m}.json" for m in c['methods']]
    for f in files:
     if f.exists():assert json.loads(f.read_text())['protocol']==PROTOCOL['fingerprint']
    if not all(f.exists() for f in files):jobs.append((s,k,c))
 write(ROOT/'provenance/resume-matrix.json',dict(protocol=PROTOCOL['fingerprint'],pending_jobs=len(jobs),runner_sha256=hashlib.sha256((ROOT/'run_revision.py').read_bytes()).hexdigest(),workers=4))
 print('PENDING JOBS',len(jobs),flush=True)
 with ProcessPoolExecutor(max_workers=4) as pool:
  fs=[pool.submit(run,j) for j in jobs]
  for i,f in enumerate(as_completed(fs),1):f.result();print('RESUMED JOBS',i,'/',len(jobs),flush=True)
