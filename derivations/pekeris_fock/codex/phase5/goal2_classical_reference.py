import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'phase3'))
import mpmath as m
from classical_reference import L,T,lc,CL_compact,E0p
m.mp.dps=40
def vals(a,p):
 out=[]
 for sign in (1,-1):
  w=(a+sign*p)/4
  out.extend(L(w+mm*m.pi/4) for mm in (-1,0,1));out.append(T(w))
 out.extend([L(p/2),L((a+m.pi/2)/2),L((m.pi/2-a)/2)])
 out.extend([lc((a+p)/4),lc((a-p)/4),a,p,a*a,p*p,m.mpf(1)])
 return out
rows=[]
for a,p in [(m.mpf('-.4'),m.mpf('-1.3')),(m.mpf('-.7'),m.mpf('-.9'))]:
 rows.append(vals(a,p))
 for ax in (0,1):rows.append([m.diff(lambda z: vals(z,p)[i] if ax==0 else vals(a,z)[i],a if ax==0 else p) for i in range(18)])
 al=a+m.pi/2;ps=p+m.pi;sg=m.cos(al/2)+m.sin(al/2);xi=m.sqrt(2)*m.sin(ps/2)
 rows.append([CL_compact(al,ps,2)-E0p(sg,xi,2,m.catalan,m.log(2))])
Path(__file__).with_name('goal2_classical_reference.txt').write_text('\n'.join(' '.join(str(float(v)) for v in r) for r in rows)+'\n')
