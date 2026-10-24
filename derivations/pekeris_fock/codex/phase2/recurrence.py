#!/usr/bin/env python3
"""Clear rows and produce exact finite Laguerre stencils.

Polynomial-ring arithmetic avoids thousands of slow expression simplifications.
"""
import sys
sys.dont_write_bytecode=True
import ast,itertools,json
from functools import lru_cache
from pathlib import Path
import sympy as S
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
import order2 as O
OUT=Path(__file__).resolve().parent
u,v,w=O.xyz
loc={str(z):z for z in (*O.xyz,O.eps,O.E,O.Z)}
loc['lam']=O.lam
def parse(z):return S.sympify(z.replace('lambda','lam'),locals=loc)
def load():
    raw=json.loads((OUT/'atom_operators.json').read_text())
    return {kind:{key:{ast.literal_eval(d):parse(f) for d,f in op.items()}
                for key,op in blocks.items()} for kind,blocks in raw.items()}

Q2=2*u*u+2*u*w+2*v*v+2*v*w+w*w
P2=(u-v)**2+2*(u+v)*w+w*w
clear={'one':u*v*w*(u+v)*Q2, 'L':u*v*w, 'x':u*v*w,
       'V':u*v*w*(u+v)*P2,'U':u*v*w,'Slog':u*v*w}
R,*gens=ring('l,m,n,eps,E,Z,lam,PI',QQ)
ell,em,en,ee,energy,charge,coupling,piconst=gens
k=S.Symbol('k')

@lru_cache(None)
def finite1(p,d,axis):
    assert p>=d
    # Source n -> target n+h, using x^d D^d L_n
    # = n_under_d sum_j (-1)^j binomial(d,j) L_(n-j).
    source=S.Symbol('source')
    acc={}
    for j in range(d+1):
        val=S.ff(source,d)*(-1)**j*S.binomial(d,j)
        for h,weight in O.old.multiplication(p-d).items():
            shift=h-j
            f=val*weight.subs(k,source-j)
            delta=-shift
            equation=(O.old.indices[axis])
            f=S.expand(f.subs(source,equation+delta))
            acc[delta]=S.expand(acc.get(delta,0)+f)
    return {delta:R.from_expr(f.xreplace({z:S.Symbol(str(z)) for z in O.old.indices}))
            for delta,f in acc.items() if f!=0}

def rec(op):
    result={}
    for d,f in op.items():
        for powers,coef in S.Poly(f,*O.xyz).terms():
            assert all(p>=der for p,der in zip(powers,d)),(powers,d)
            coef=R.from_expr(coef.xreplace({S.pi:S.Symbol('PI'),O.lam:S.Symbol('lam')}))
            maps=[finite1(powers[i],d[i],i) for i in range(3)]
            for entries in itertools.product(*(z.items() for z in maps)):
                offset=tuple(t[0] for t in entries)
                value=coef*entries[0][1]*entries[1][1]*entries[2][1]
                result[offset]=result.get(offset,R.zero)+value
    return {d:f for d,f in result.items() if f}

def main():
    raw=load(); cleared={};stats={};allrecs={}
    # Exact restart: existing blocks are reusable only if their cleared
    # differential operators match coefficient by coefficient.
    previous_ops=json.loads((OUT/'cleared_operators.json').read_text()) if (OUT/'cleared_operators.json').exists() else {}
    previous_recs=json.loads((OUT/'codex_recurrence.json').read_text()) if (OUT/'codex_recurrence.json').exists() else {}
    previous_stats=json.loads((OUT/'counts.json').read_text()) if (OUT/'counts.json').exists() else {}
    for kind in ('S2P','C1','C2'):
        cleared[kind]={};stats[kind]={};allrecs[kind]={}
        for key,op in raw[kind].items():
            target=key.strip('()').split(',')[0]
            op={d:S.Poly(S.cancel(clear[target]*f),*O.xyz).as_expr() for d,f in op.items()}
            cleared[kind][key]=op
            saved=previous_ops.get(kind,{}).get(key)
            if saved is not None and saved=={str(d):str(f) for d,f in op.items()}:
                allrecs[kind][key]=previous_recs[kind][key]
                stats[kind][key]=previous_stats[kind][key]
                print('Reusing identical cleared operator',kind,key,flush=True)
                continue
            print('Deriving',kind,key,flush=True)
            stencil=rec(op)
            fixed={d:f.compose(energy,-ee**2) for d,f in stencil.items()}
            fixed={d:f for d,f in fixed.items() if f}
            stat={'general_count':len(stencil),'fixed_count':len(fixed),
                  'general_shell_range':[min(map(sum,stencil)),max(map(sum,stencil))],
                  'fixed_shell_range':[min(map(sum,fixed)),max(map(sum,fixed))],
                  'coordinate_bounds':[[min(d[i] for d in stencil),max(d[i] for d in stencil)] for i in range(3)],
                  'differential_monomials':sum(len(S.Poly(f,*O.xyz).terms()) for f in op.values())}
            stats[kind][key]=stat
            print(stat,flush=True)
            # Exact expanded polynomials; compatible offset convention with previous JSON.
            allrecs[kind][key]={str(d):str(f).replace('lam','lambda').replace('PI','pi')
                               for d,f in sorted(stencil.items())}
    (OUT/'cleared_operators.json').write_text(json.dumps({kind:{key:{str(d):str(f) for d,f in op.items()}
        for key,op in blocks.items()} for kind,blocks in cleared.items()},indent=2)+'\n')
    (OUT/'codex_recurrence.json').write_text(json.dumps(allrecs,indent=2)+'\n')
    (OUT/'counts.json').write_text(json.dumps(stats,indent=2)+'\n')
    (OUT/'clearing.json').write_text(json.dumps({key:str(S.factor(f)) for key,f in clear.items()},indent=2)+'\n')

if __name__=='__main__':main()
