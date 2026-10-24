"""Deterministic independent reference for the Julia port; no source writes."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'phase3'))
import mpmath as m
from panel_connection import rhs,moments_from_state
m.mp.dps=40
out=[]
for b,y in [(m.mpf('.3'),m.mpf('-.6')),(m.mpf('-.27'),m.mpf('-.7'))]:
 X=[[[m.mpf((f*39+i*3+k)%17-8)/19 for k in range(3)] for i in range(13)] for f in range(2)]
 for axis in (0,1):
  vals=[]
  for f in range(2):
   vals.extend(v for row in rhs(b,y,axis,'A' if f==0 else 'D',X[f]) for v in row)
  out.append(vals)
 out.append(list(moments_from_state(b,y,*X)))
p=Path(__file__).with_name('goal2_port_reference.txt')
p.write_text('\n'.join(' '.join(str(float(m.re(v)))+' '+str(float(m.im(v))) for v in row) for row in out)+'\n')
