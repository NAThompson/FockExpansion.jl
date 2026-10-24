"""Independent panel evaluation + attributed earlier classical formula vs Julia source.
The classical formula is a copied read-only earlier-workspace formula, NOT a new
independent derivation. New Pfaffian state evaluation is in panel_connection.py.
"""
from pathlib import Path
import mpmath as m,json
from panel_connection import state,moments_from_state
from classical_reference import CL_compact
m.mp.dps=48
rows=[]
for line in Path(__file__).with_name('repo_values.tsv').read_text().splitlines():
 fields=line.split('\t');a,t,Z,E,c=map(m.mpf,fields[:5]); refs=list(map(m.mpf,fields[5:]))
 af=min(a,m.pi-a);r1,r2=m.cos(af/2),m.sin(af/2)
 v=m.sin(af)*m.cos(t);x=m.sqrt(1-v);s=r1+r2
 cl=CL_compact(af,m.acos(v),Z)
 st=E*(Z*s*(2+r1*r2)/18-(6-x*x)*x/72)-c*(Z*s*v/2-(6-5*x*x)*x/12)
 ell=0
 if Z:
  for k,(b,y) in enumerate([(m.cos(af),-x/m.sqrt(2)),(-m.cos(af),-x/m.sqrt(2)),(-v,-r1),(-v,-r2)]):
   ns=moments_from_state(b,y,state(b,y,'A'),state(b,y,'D'))
   if k<2:
    c0=-m.sqrt(2)*Z*(9*Z*b*b-18*Z-8*b*b+8)/576
    c2=-m.sqrt(2)*Z*(1-b)*(90*Z*y*y-45*Z-18*b-32*y*y+16)/288
   else:c0=Z*Z*(1-b*b)/18;c2=-2*Z*Z*(1-b)*(2*y*y-1)/9
   ell+=c0*ns[0]+c2*ns[1]
 vals=[cl,st,ell,cl+st+ell]
 errors=[abs(v-r) for v,r in zip(vals,refs[:4])]
 floaterr=abs(vals[-1]-refs[-1])
 print('REPO',fields[:5],'part errors',*[m.nstr(e,8) for e in errors],'Float64',m.nstr(floaterr,8),flush=True)
 assert max(errors)<m.mpf('1e-40'),errors
 assert floaterr<m.mpf('2e-12'),floaterr
 rows.append(dict(parameters=fields[:5],errors=list(map(str,errors)),float_error=str(floaterr)))
assert len(rows)==3
Path(__file__).with_name('repo_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
print('PASS: all four parts at three parameter sets, both distance orderings and Z=0')
