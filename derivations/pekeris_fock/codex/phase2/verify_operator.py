#!/usr/bin/env python3
"""Direct numerical polylog differentiation against independent block action."""
import sys
sys.dont_write_bytecode=True
import ast,json,random
from pathlib import Path
import sympy as S
import mpmath as mp
import order2 as O
import recurrence as R
from verify_atoms import atoms
OUT=Path(__file__).resolve().parent
mp.mp.dps=40

def main():
    ops=R.load()
    pars={O.eps:S.Rational(3,2),O.E:S.Rational(-7,3),O.Z:2,O.lam:S.Rational(4,5)}
    sig=sum(O.xyz)
    maxerr=mp.mpf(0)
    for seed in (3,11,17):
        if seed==17:pars[O.E]=-pars[O.eps]**2
        pbase={d:S.lambdify(O.xyz,c.subs(pars),'mpmath') for d,c in O.old.derive()['P'].items()}
        rng=random.Random(seed)
        F=[]
        for k in range(3):
            row=[]
            for atom in O.names:
                row.append(sum(S.Rational(rng.randint(-3,3),rng.randint(1,5))*
                    S.prod(S.laguerre(n,z) for n,z in zip(source,O.xyz))
                    for source in [(0,0,0),(1,0,1),(0,2,0)]))
            F.append(row)
        ff=[[S.lambdify(O.xyz,f,'mpmath') for f in row] for row in F]
        def total(*z):
            g=atoms(*z); log=mp.log(sum(z))
            return sum(log**k*sum(g[i]*ff[k][i](*z) for i in range(O.rank)) for k in range(3))
        # Separate every log sector; block application uses symbolic polynomial derivatives.
        projected=[]
        for k in range(3):
            row=[]
            for i in range(O.rank):
                value=0
                for j in range(O.rank):
                    key=f'({O.names[i]},{O.names[j]})'
                    for kind,shift,prefactor in [('S2P',0,1),('C1',1,k+1),('C2',2,(k+1)*(k+2))]:
                        if k+shift>=3:continue
                        for d,c in ops[kind].get(key,{}).items():
                            value+=prefactor*c*S.diff(F[k+shift][j],*sum(([z,n] for z,n in zip(O.xyz,d)),[]))
                row.append(S.lambdify(O.xyz,value.subs(pars),'mpmath'))
            projected.append(row)
        point=(mp.mpf(2)/5,mp.mpf(7)/10,mp.mpf(3)/2)
        direct=sum(fun(*point)*mp.diff(total,point,d) for d,fun in pbase.items())*sum(point)**2
        g=atoms(*point); log=mp.log(sum(point))
        expected=sum(log**k*sum(g[i]*projected[k][i](*point) for i in range(O.rank)) for k in range(3))
        err=abs(direct-expected)/max(mp.mpf(1),abs(direct))
        maxerr=max(maxerr,err)
        print('seed',seed,'full polylog/log/Laguerre Hamiltonian relative error',mp.nstr(err,8),flush=True)
        assert err<mp.mpf('1e-32')
    print('PASS: 3 full',O.rank,'atom, three-log-sector random actions (one E=-eps^2); max',mp.nstr(maxerr,8))

if __name__=='__main__':main()
