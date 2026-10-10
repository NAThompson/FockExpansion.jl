"""Regenerate and check every residual from source, without cached results."""
from pathlib import Path
import subprocess,sys,hashlib,json,time
HERE=Path(__file__).resolve().parent
sectors=['Z0','Z1','Z2','Z3','E','EZ','a21','a21Z']
started=time.time()
def run(script,*args):
 print('Running',script,*args,flush=True)
 result=subprocess.run([sys.executable,str(HERE/script),*args],cwd=HERE,check=True,capture_output=True,text=True)
 (HERE/(script+'.log')).write_text(result.stdout+result.stderr)
 return result.stdout
run('panel_operator.py')
run('exact_residual.py','--emit-only')
out=run('flint_check.py',*sectors)
for tag in sectors:assert tag+' EXACT ZERO' in out
assert 'NONZERO' not in out
files=sorted(HERE.glob('*.py'))+sorted((HERE/'inputs').rglob('*.py'))+sorted((HERE/'inputs').rglob('*.json'))
report={'status':'passed','sectors':sectors,'seconds':time.time()-started,'sha256':{str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
(HERE/'certificate.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS: all eight sectors vanish exactly, freshly regenerated.',flush=True)
