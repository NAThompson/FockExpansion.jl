"""Actual-function check of the quotient and old-atom pullback."""
import sys
from pathlib import Path
here=Path(__file__).resolve().parent
sys.path.insert(0,str(here.parent/'phase3'))
sys.path.insert(0,str(here.parent/'phase2'))
import mpmath as m
from panel_connection import state
from verify_atoms import atoms
from goal2_classical_reference import vals
m.mp.dps=35
lines=[]
for a,p in [(m.mpf('-.4'),m.mpf('-1.3')),(m.mpf('-.6'),m.mpf('-1.1'))]:
 r=m.cos(a/2+m.pi/4);t=m.sin(a/2+m.pi/4);c=m.sqrt(2)*m.cos(p/2)
 xx=vals(a,p);zz=list(xx)
 panels=[(-m.sin(a),-m.cos(p/2)),(m.sin(a),-m.cos(p/2)),(m.cos(p),-r),(m.cos(p),-t)]
 for b,y in panels:
  A,D=[state(b,y,f) for f in ('A','D')]
  xx.extend(v for X in (A,D) for row in X for v in row)
  for j in range(11):zz.extend([A[j][0],A[j][1],A[j][2],D[j][1],D[j][2]])
 def old(aa,pp):
  rr=m.cos(aa/2+m.pi/4);tt=m.sin(aa/2+m.pi/4);cc=m.sqrt(2)*m.cos(pp/2)
  return atoms((cc-rr+tt)/4,(cc+rr-tt)/4,(rr+tt-cc)/2)
 av=old(a,p);xx.extend([av[1],av[4]]);zz.extend([av[1],av[4]])
 lines.extend([xx,zz])
 for axis in (0,1):lines.append([m.diff(lambda q:old(q,p)[j] if axis==0 else old(a,q)[j],a if axis==0 else p) for j in (1,4)])
path=here/'goal2_combined_reference.txt'
path.write_text('\n'.join(' '.join(str(float(m.re(v)))+' '+str(float(m.im(v))) for v in row) for row in lines)+'\n')
print('Saved two actual 332-state vectors, two 240-state quotient vectors and eight old-atom derivatives')
