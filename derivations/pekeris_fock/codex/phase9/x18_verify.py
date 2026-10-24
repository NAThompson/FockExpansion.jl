"""Small serial identity checks; no repository code execution or outside writes."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import sympy as S
import mpmath as m
import json, resource, time, hashlib
from pathlib import Path
start=time.time(); g=S.symbols('g', positive=True)
G=((g-S.pi)*S.cos(3*g)+S.sin(3*g)/6)/(4*S.pi**2*S.sin(g))
G2=((g*g-2*S.pi*g)/(48*S.pi**2)+S.Rational(1,72)+1/(864*S.pi**2))*S.sin(3*g)/S.sin(g)
A=lambda f:S.diff(f,g,2)+2*S.cos(g)/S.sin(g)*S.diff(f,g)+8*f
zero=lambda f:S.simplify(S.expand_trig(f))==0
checks={}
checks['kernel_off_diagonal']=zero(A(G)+3*(4*S.cos(g)**2-1)/(2*S.pi**2))
checks['iterated_kernel']=zero(A(G2)-G)
checks['delta_flux']=S.limit(4*S.pi*g*g*S.diff(G,g),g,0)==1
s,H,k,a,d=S.symbols('s H k a d'); T=s**4+k*s*s+H
checks['panel_discriminant']=S.factor(S.discriminant(T,s)-16*H*(k*k-4*H)**2)==0
I=k*k+12*H; J=72*k*H-2*k**3
checks['panel_j']=S.factor(6912*I**3/(4*I**3-J**2)-16*I**3/(H*(k*k-4*H)**2))==0
F=(s*s-4)*((s-2*a)**2-4*d*d)
checks['kernel_discriminant']=S.factor(S.discriminant(F,s)-65536*d*d*((a-d)**2-1)**2*((a+d)**2-1)**2)==0
m.mp.dps=45
Gm=S.lambdify(g,G,'mpmath'); G2m=S.lambdify(g,G2,'mpmath')
moments=[]
for fun in (Gm,G2m):
 for n in (0,1,2,3,4):
  v=4*m.pi*m.quad(lambda t:fun(t)*m.sin((n+1)*t)*m.sin(t),[0,m.pi/2,m.pi])
  expected=0 if n==2 else 1/(8-m.mpf(n)*(n+2))
  if fun==G2m:expected=expected**2
  moments.append(float(abs(v-(n+1)*expected)))
checks['ten_spectral_moments']=max(moments)<1e-40
# Exact generic smoothness at two finite-field points, not a rank measurement.
ff=[]
p=1000003
for av,dv in [(S.Rational(1,5),S.Rational(1,7)),(S.Rational(2,9),S.Rational(1,6))]:
 poly=S.Poly(F.subs({a:av,d:dv}),s)
 coeff=[int(v.p%p)*pow(int(v.q),-1,p)%p for v in poly.all_coeffs()]
 fp=S.Poly.from_list(coeff,s,modulus=p)
 ff.append(dict(a=str(av),d=str(dv),prime=p,gcd_degree=S.gcd(fp,fp.diff()).degree()))
checks['two_ff_smooth_fibres']=all(v['gcd_degree']==0 for v in ff)
root=Path(__file__).resolve().parents[2]; repo=Path('/Users/nathompson7/FockExpansion.jl')
inputs=[repo/'src/green_k4.jl',repo/'src/fourth_order.jl',repo/'src/feynman_k4.jl',repo/'derivations/k4/README.md',root/'codex/phase3/structure.md']
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
out=dict(checks=checks,max_moment_error=max(moments),finite_field_smoothness=ff,rank_measured=False,elapsed_seconds=time.time()-start,peak_rss_native=rss,input_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs})
Path(__file__).with_name('x18_verify.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));assert all(checks.values())
