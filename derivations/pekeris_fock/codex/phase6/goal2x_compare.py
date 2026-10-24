"""Read-only Claude comparison, run only after goal2x_step1 was saved."""
import sys,os,json,hashlib
sys.dont_write_bytecode=True
os.environ['GOAL2X_PURE']='1' # C accelerator is intentionally specific to the independent prime.
import goal2x_probe as G
from pathlib import Path
assert (G.OUT/'goal2x_step1_hashes.json').exists()
paths=['claude6/guess_probe.py','claude6/guess_line.py','claude6/guess_line2.py','claude6/ratrec.py','claude6/guess_line.json','claude6/line_a.json']
hashes={}
for name in paths:
 p=G.ROOT/name;blob=p.read_bytes();hashes[name]=hashlib.sha256(blob).hexdigest()
 (G.OUT/('goal2x_claude_'+p.name)).write_bytes(blob)
d=json.loads((G.OUT/'goal2x_claude_line_a.json').read_text());targets=[tuple(t) for t in d['targets']]
G.PRIME=998244353;results=[]
z=pow(3,(G.PRIME-1)//8,G.PRIME);rt=(z+pow(z,-1,G.PRIME))%G.PRIME
sr=1 if int(G.sp.sqrt_mod(2,G.PRIME))==rt else -1
si=1 if int(G.sp.sqrt_mod(-1,G.PRIME))==pow(3,(G.PRIME-1)//4,G.PRIME) else -1
for k in [-2,-1]:
 u=d['xs'][k];point=[3+u,4+2*u]
 got=G.relations(*point,pi=101,si=si,sr=sr,targets=targets)
 expected=d['vals'][k]
 assert got==expected,(u,[(a,b) for a in range(len(got)) for b in range(45) if got[a][b]!=expected[a][b]])
 results.append(dict(u=u,point=point,entries=len(targets)*45,match=True))
 print('PASS Claude corrected target samples',point,len(targets)*45,'coefficients',flush=True)
(G.OUT/'goal2x_comparison.json').write_text(json.dumps(dict(input_hashes=hashes,results=results,note='No reconstructed bivariate psi30 coefficients were available; cross-checked saved modular solves instead.'),indent=2)+'\n')
