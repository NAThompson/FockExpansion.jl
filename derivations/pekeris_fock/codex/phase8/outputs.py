"""Exact moment and homogeneous order-three output rows."""
import json
from pathlib import Path
import sympy as S
from rational_panel import names
O=Path(__file__).resolve().parent
H,k,x,y,c,Z=S.symbols('H k x y c Z');pi=S.pi
v={n:S.Symbol('X'+str(i)) for i,n in enumerate(names)}
def I2(fam,n):
 suf=('A' if n==0 else fam)+str(n)
 a=v['E_'+suf]-k*v['J_'+suf]-(2*H-k*k)*v['Q2_'+suf]+k*H*v['Q0_'+suf]
 if n:
  prev=('A' if n==1 else fam)+str(n-1)
  a-=v['K_'+prev]-v['U1_'+prev] if fam=='A' else H*v['K_'+prev]-H*H*v['V1_'+prev]
 return S.expand(a)
F0=pi/2*v['J_A0']-2*v['J_A1']+4/pi*(v['J_A2']-v['J_B2']/H)
Im=-2*v['Q2_B2']-k*v['Q0_B2']+v['FP_B1']-v['V1_B1']-v['E_B2']
F2=4/pi*Im-pi/2*I2('A',0)+2*I2('A',1)-4/pi*I2('A',2)
Q=x*x+y*y;e=2*Q-c*c
out=Z*Z*c*e*F0/18-2*Z*Z*e*(x*x-y*y)*F2/(9*c)
out=out.subs({H:c*c/e,k:2*(x*x-y*y)/e},simultaneous=True)
def row(a):return [[i,str(S.cancel(S.expand(a).coeff(v[n])))] for i,n in enumerate(names) if S.expand(a).coeff(v[n])!=0]
(O/'outputs.json').write_text(json.dumps(dict(F0=row(F0),F2=row(F2),N0='sqrt(x^2+y^2)/c * F0',N2='sqrt(x^2+y^2)/c * F2',panel3_R3psi30=row(out),panel4='exchange x and y',Z='nuclear charge (helium Z=2)'),indent=2)+'\n')
print('Exported exact F0,F2 and rational homogeneous-degree-three panel output:',len(row(out)),'terms')
