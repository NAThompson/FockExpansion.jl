"""Independent adaptive-mpmath check of the Julia x16-state fit cache."""
import sys,ast,json,time,resource
from pathlib import Path
sys.dont_write_bytecode=True
import mpmath as m
out=Path(__file__).resolve().parent;root=out.parents[1];sys.path.insert(0,str(root/'codex/phase8'))
import rational_panel as R
m.mp.dps=62
node=next(n for n in ast.parse((root/'codex/phase8/verify.py').read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='actual')
ns={'m':m,'R':R};exec(compile(ast.Module(body=[node],type_ignores=[]),'phase8_actual','exec'),ns)
actual=ns['actual'];ids=[int(l) for l in (out/'x19_fit40_big_ids.txt').read_text().splitlines()]
data=(root/'claude9/spec40.txt').read_text().splitlines();cache=(out/'x19_fit40_big_atoms.txt').read_text().splitlines()
names=json.loads((root/'codex/phase8/x16_summary.json').read_text())['1']['reachable']
worst=m.mpf(0);count=0;start=time.time()
for row in [0,12]:
 a,t,*_=map(m.mpf,data[ids[row]-1].split());x=m.cos(a/2);y=m.sin(a/2);c=m.sqrt(1-m.sin(a)*m.cos(t));Q=x*x+y*y;e=2*Q-c*c
 maps=[(x*x/(y*y),(c*c-Q)/(y*y)),(y*y/(x*x),(c*c-Q)/(x*x)),(c*c/e,2*(x*x-y*y)/e),(c*c/e,2*(y*y-x*x)/e)]
 vals=list(map(m.mpf,cache[row].split()))
 for p,(H,k) in enumerate(maps):
  ref=dict(zip(R.names,actual(H,k)))
  for j,n in enumerate(names):
   err=abs(vals[23+p*23+j]-ref[n])/(1+abs(ref[n]));worst=max(worst,err);count+=1
assert worst<m.mpf('1e-51'),worst
result=dict(comparisons=count,max_relative_scale_error=str(worst),elapsed_seconds=time.time()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
(out/'x19_panelcheck.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
