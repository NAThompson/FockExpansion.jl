"""Exact quartic involution and actual fixed-endpoint physical moments."""
import sympy as S
import mpmath as m
from pathlib import Path
import json
m.mp.dps=50
h,s,Y=S.symbols('h s Y',nonzero=True);b=(h*h-1)/(h*h+1)
T=lambda x:(1-b)*x**4+(4*Y*Y-2)*x*x+1+b
assert S.cancel(T(h/s)-h*h/s**4*T(s))==0
print('PASS exact symbolic quartic involution',flush=True)
worst=m.mpf(0);cases=[]
for H in map(m.mpf,['.5','.8','1.2','2']):
 B=(H*H-1)/(H*H+1);yy=m.mpf('-.6')
 def g(q,k):
  A=m.pi/2-m.atan(q);D=m.atan(q/H);W=m.sqrt(1+B+(4*yy*yy-2)*q*q+(1-B)*q**4)
  return 2/m.pi*((A*A-D*D) if k==0 else ((H*D/q)**2-(q*A)**2))/W
 for k in [0,2]:
  for q in map(m.mpf,['.2','.7','1.3']):
   assert abs(g(H/q,k)*H/q**2+g(q,k))<m.mpf('1e-45')
  segment=m.quad(lambda q:g(q,k),[1,m.sqrt(H),H]);worst=max(worst,abs(segment))
  fixed=m.quad(lambda q:g(q,k),[0,1]);moving=m.quad(lambda q:g(q,k),[0,min(1,H)])
  assert abs(fixed-moving)<m.mpf('1e-44')
  assert abs(segment)<m.mpf('1e-44')
  cases.append(dict(h=str(H),moment=k,fixed=str(fixed),moving=str(moving),segment=str(segment)))
print('PASS 24 pointwise involution checks and 8 actual moment/endpoint comparisons; max segment',m.nstr(worst,6),flush=True)
# This also verifies that replacing the old primitive by fixed-endpoint J
# removes the first-derivative jump; physical moments did not change.
yy=m.mpf('-.6')
def J(B):return m.quad(lambda q:1/m.sqrt(1+B+(4*yy*yy-2)*q*q+(1-B)*q**4),[0,1])
e=m.mpf('1e-12');jump=(J(e)-2*J(0)+J(-e))/e
assert abs(jump)<m.mpf('1e-10')
print('PASS fixed-endpoint primitive has zero derivative jump (finite-difference error',m.nstr(jump,4)+')',flush=True)
Path(__file__).with_name('fixed_endpoint_check.json').write_text(json.dumps(dict(cases=cases,max_segment=str(worst),fixed_primitive_jump=str(jump)),indent=2)+'\n')
