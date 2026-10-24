#!/usr/bin/env python3
"""Independent direct actual-function check of Claude's saved stencils.
No Claude verifier, connection or operator generator is imported.
"""
import sys
sys.dont_write_bytecode=True
import ast,json,random,math
from pathlib import Path
import mpmath as mp
import sympy as S
import order2 as O
from verify_atoms import atoms
from verify_stencils import compile_coefficient
from fractions import Fraction as Q
ROOT=Path(__file__).resolve().parents[2]
mp.mp.dps=45
raw=json.loads((ROOT/'claude2/order2_recurrence.json').read_text())
# In the open folded physical chart U=Om+pi*lb, L=la,
# V=(beP+pi*P/2)/r12, Slog=-lR.
def ca(*z):
    a,b,d=2*z[1]+z[2],2*z[0]+z[2],2*(z[0]+z[1])
    q=a*a+b*b; p=2*q-d*d; P=mp.sqrt(p)
    lb=mp.log((a-b+d)**2/(d*P+(a+b)*(a-b)))
    one,L,al,V,U,Slog=atoms(*z)
    return dict(Om=U-mp.pi*lb,al=al,beP=d*V-mp.pi*P/2,lR=-Slog,la=L,lb=lb,one=one,P=P)
# eps=1/4 makes the physical distances exactly (2v+w,2u+w,2(u+v)).
eps=Q(1,4)
pars={O.eps:S.Rational(1,4),O.E:-S.Rational(1,16),O.Z:2,O.lam:1}
base={d:S.lambdify(O.xyz,f.subs(pars),'mpmath') for d,f in O.old.derive()['P'].items()}
D={k:S.lambdify(O.xyz,S.sympify(f,locals={**dict(zip(('u','v','w'),O.xyz)),'eps':O.eps}).subs(pars),'mpmath') for k,f in raw['D'].items()}
points=[tuple(map(mp.mpf,p)) for p in [('0.4','0.7','1.5'),('0.2','0.9','0.6')]]
sources=[((0,0,0),Q(2,3)),((1,0,1),Q(-3,7)),((0,2,0),Q(1,5))]
F=sum(S.Rational(c.numerator,c.denominator)*S.prod(S.laguerre(n,z) for n,z in zip(j,O.xyz)) for j,c in sources)
ff=S.lambdify(O.xyz,F,'mpmath')
# Pre-evaluate all expanded stencil coefficients exactly; then reconstruct.
compiled={}
for kind in ('sigma2P','C1','C2'):
    compiled[kind]={}
    for key,stencil in raw[kind].items():
        outputs={}
        for off,expr in stencil.items():
            off=ast.literal_eval(off)
            code=None
            for source,c in sources:
                target=tuple(source[i]-off[i] for i in range(3))
                if min(target)<0:continue
                if code is None:code=compile_coefficient(expr)
                val=eval(code,{'__builtins__':{}},dict(Q=Q,eps=eps,epsilon=eps,l=target[0],m=target[1],n=target[2]))*c
                outputs[target]=outputs.get(target,Q(0))+val
        compiled[kind][ast.literal_eval(key)]={j:c for j,c in outputs.items() if c}
print('Parsed',sum(len(v) for v in compiled.values()),'blocks;',sum(len(v) for k in ('sigma2P','C1','C2') for v in raw[k].values()),'offsets',flush=True)
maxerr=mp.mpf(0)
for point in points:
    g=ca(*point); log=mp.log(sum(point)); L={}
    for kind,blocks in compiled.items():
        L[kind]={}
        for (target,source),pol in blocks.items():
            value=sum(mp.mpf(c.numerator)/c.denominator*math.prod(mp.laguerre(jj,0,z) for jj,z in zip(j,point)) for j,c in pol.items())
            L[kind][target,source]=value
    for source in raw['atoms']:
        for power in (0,1,2):
            def total(*z):return ca(*z)[source]*ff(*z)*mp.log(sum(z))**power
            direct=sum(fun(*point)*mp.diff(total,point,d) for d,fun in base.items())*sum(point)**2
            predicted=mp.mpf(0)
            for kind,pref in [('sigma2P',log**power),('C1',power*log**(power-1)),('C2',power*(power-1)*log**(power-2))]:
                predicted+=pref*sum(g[target]*v/D[target](*point) for (target,src),v in L[kind].items() if src==source)
            err=abs(direct-predicted)/max(1,abs(direct));maxerr=max(maxerr,err)
            assert err<mp.mpf('1e-32'),(point,source,power,err)
        print('PASS actual atom',source,'three log powers at',tuple(map(str,point)),flush=True)
print('PASS 48 actual Hamiltonian/stencil actions; max relative error',mp.nstr(maxerr,12))
# Counts and clearing factors, unabridged audit.
for kind in ('sigma2P','C1','C2'):
    for key,st in raw[kind].items():
        shells=[sum(ast.literal_eval(o)) for o in st]
        print('COUNT',kind,key,len(st),min(shells),max(shells))
print('ROW FACTORS',raw['D'])
