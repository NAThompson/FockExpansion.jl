import mpmath as m,json
from pathlib import Path
from panel_connection import state,rhs,moments_from_state,direct_moments
m.mp.dps=42
results=[]
for bs,ys in [('0.3','-0.6'),('-0.27','-0.7')]:
 b,y=map(m.mpf,(bs,ys)); states={f:state(b,y,f) for f in ('A','D')}
 ref=direct_moments(b,y);actual=moments_from_state(b,y,states['A'],states['D'])
 err=max(abs(a-c) for a,c in zip(actual,ref));assert err<m.mpf('1e-36')
 print('MOMENTS',bs,ys,'max error',m.nstr(err,8),flush=True)
 for family in ('A','D'):
  for axis in (0,1):
   # Independent centered difference with two step sizes, not a Pfaffian integrator.
   errors=[]
   for step in (m.mpf('1e-9'),m.mpf('5e-10')):
    p=state(b+step if axis==0 else b,y+step if axis==1 else y,family)
    n=state(b-step if axis==0 else b,y-step if axis==1 else y,family)
    R=rhs(b,y,axis,family,states[family])
    err=max(abs((p[i][j]-n[i][j])/(2*step)-R[i][j])/(1+abs(R[i][j])) for i in range(13) for j in range(3))
    errors.append(err)
   assert errors[-1]<m.mpf('1e-17'),errors
   assert errors[-1]<errors[0]/3,errors
   print('CONNECTION',bs,ys,family,axis,'errors',*[m.nstr(e,8) for e in errors],flush=True)
   results.append(dict(b=bs,y=ys,family=family,axis=axis,errors=list(map(str,errors))))
Path(__file__).with_name('connection_checks.json').write_text(json.dumps(results,indent=2)+'\n')
print('PASS: 312 complex component/jet equations per step size; four moment reconstructions',flush=True)
