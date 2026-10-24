#!/usr/bin/env python3
"""Independent order-two differential module; no Claude inputs.

All output is confined to this directory. Imports the *previous independent*
P derivation read-only; prevents Python bytecode writes outside phase2.
"""
import sys
sys.dont_write_bytecode = True
import importlib.util
import itertools
import json
from pathlib import Path
import sympy as S

OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('old', OUT.parent/'codex_derive.py')
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
u,v,w = old.x
xyz = (u,v,w)
eps,E,Z,lam = old.eps,old.E,old.Z,old.lam
names = ('one','L','x','V','U','Slog')
rank=len(names)
h = S.Symbol('h')
pi = S.pi
can = S.cancel

def geometry(coords=xyz):
    if coords == xyz:
        a,b,d = 2*v+w,2*u+w,2*(u+v)
    else:
        a,b,d = coords
    q=a*a+b*b
    p=2*q-d*d
    sn=2*a*b/q
    cs=(a*a-b*b)/q
    V=(q-d*d)/q
    t=d*d/q
    f=p/d**2  # h^2
    den=sn*sn-V*V
    A=cs+(a+b)*d/q
    return a,b,d,q,p,sn,cs,V,t,f,den,A

def connection(coords=xyz):
    a,b,d,q,p,sn,cs,V,t,f,den,A=geometry(coords)
    mats=[]
    for z in coords:
        D=lambda F:S.diff(F,z)
        kh=D(p)/(2*p)-D(d)/d
        kg=D(p)/(2*p)+D(d)/d-D(q)/q
        J0=2*D(A)/A-D(sn)/(1+sn)-(t*t*f*kg-cs*D(cs))/den
        J1=-t*(D(cs)-cs*kg)/den
        T_B=2*t*(cs*kg-D(cs))/den
        T_x=2*(sn*D(V)-V*D(sn))/den
        C=S.zeros(6)
        C[1,1]=kh
        C[0,2]=D(a+b+d)/(a+b+d)-D(q)/(2*q)
        C[0,3]=2*(a*D(b)-b*D(a))/q
        C[4,4]=kh
        C[0,4]=q*D(V)/(d*d)
        C[0,5]=pi*J0
        C[1,5]=pi*J1
        C[3,5]=T_x
        C[4,5]=T_B
        C=C.applyfunc(can)
        # V=h*(beta+pi/2): the explicit h generator cancels identically.
        assert can(C[1,5]-pi*C[4,5]/2)==0
        small=S.zeros(rank)
        small[0,1]=C[0,2]; small[0,2]=C[0,3]
        small[3,3]=C[4,4]; small[0,3]=C[0,4]
        small[0,4]=C[0,5]; small[2,4]=C[3,5]; small[3,4]=C[4,5]
        # log(R) -> log(sigma) adds Slog=log((a+b)/R).
        small[0,5]=can(D(a+b)/(a+b)-D(q)/(2*q))
        mats.append(small)
    return mats

def atom_operators(C):
    """P(g_j F) = sum_i g_i P_ij F; column-oriented connection."""
    base=old.derive()['P']
    M=S.zeros(3)
    drift=[]
    for i in range(3):
        e=tuple(int(k==i) for k in range(3))
        drift.append(base[e])
        for j in range(3):
            d=tuple(int(k==i)+int(k==j) for k in range(3))
            M[i,j]=base[d] if i==j else base[d]/2
    second={(i,j):(C[j].diff(xyz[i])+C[i]*C[j]).applyfunc(can)
            for i in range(3) for j in range(3)}
    blocks={}
    for target in range(rank):
        for source in range(rank):
            op={d:f for d,f in base.items()} if target==source else {}
            for j in range(3):
                d=tuple(int(k==j) for k in range(3))
                val=2*sum(M[i,j]*C[i][target,source] for i in range(3))
                op[d]=can(op.get(d,0)+val)
            val=sum(drift[i]*C[i][target,source] for i in range(3))
            val+=sum(M[i,j]*second[i,j][target,source]
                     for i in range(3) for j in range(3))
            op[(0,0,0)]=can(op.get((0,0,0),0)+val)
            op={d:can(f) for d,f in op.items() if f!=0}
            if op:blocks[(target,source)]=op
    sig=u+v+w
    # C1 extended: sigma sum drift - sum M + 2 sigma sum M C.
    base1=old.derive()['C1']; base2=old.derive()['C2']
    c1={};c2={}
    for target in range(rank):
        for source in range(rank):
            op=dict(base1) if target==source else {}
            val=2*sig*sum(M[i,j]*C[j][target,source]
                          for i in range(3) for j in range(3))
            op[(0,0,0)]=can(op.get((0,0,0),0)+val)
            op={d:f for d,f in op.items() if f!=0}
            if op:c1[target,source]=op
            if target==source:c2[target,source]=dict(base2)
    return {'P':blocks,'S2P':{key:{d:can(sig**2*f) for d,f in op.items()}
                             for key,op in blocks.items()},'C1':c1,'C2':c2}

def serial(ops):
    return {name:{f'({names[i]},{names[j]})':{str(d):str(S.factor(f))
                for d,f in op.items()} for (i,j),op in blocks.items()}
            for name,blocks in ops.items()}

def main():
    print('Constructing rational connection',flush=True)
    C=connection()
    for i in range(3):
        for j in range(i):
            curvature=C[j].diff(xyz[i])-C[i].diff(xyz[j])+C[i]*C[j]-C[j]*C[i]
            assert all(can(z)==0 for z in curvature)
    print('Flat connection: 108 entries checked exactly',flush=True)
    (OUT/'connection.json').write_text(json.dumps({str(z):{
        f'({names[i]},{names[j]})':str(S.factor(c[i,j]))
        for i in range(rank) for j in range(rank) if c[i,j]!=0}
        for z,c in zip(xyz,C)},indent=2)+'\n')
    print('Constructing atom operators',flush=True)
    ops=atom_operators(C)
    (OUT/'atom_operators.json').write_text(json.dumps(serial(ops),indent=2)+'\n')
    print({k:len(b) for k,b in ops.items()},flush=True)
    for kind,blocks in ops.items():
        for key,op in blocks.items():
            assert all(not f.has(h) for f in op.values())
    print('All block coefficients rational over Q(pi,eps,E,Z,lambda,u,v,w)',flush=True)

if __name__=='__main__':main()
