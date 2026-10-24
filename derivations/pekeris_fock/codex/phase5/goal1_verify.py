#!/usr/bin/env python3
"""Small serial checks of the physical Pekeris pencil and backward selector.
Only Codex recurrence JSON is input; no Claude assembly or large arrays.
Run: PYTHONDONTWRITEBYTECODE=1 venv/bin/python codex/phase5/goal1_verify.py
"""
import ast
import hashlib
import json
import resource
import sys
from pathlib import Path
import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
mp.mp.dps = 45
EPS = mp.mpf('1.704')
E0 = mp.mpf('-2.903724377034119598311159245194404446696925310')
l, m, n, eps, E, Z, lam = sp.symbols('l m n eps E Z lambda')
source = ROOT / 'codex/codex_recurrence.json'
raw = json.loads(source.read_text())['P']
ops = []
for off, expr in raw.items():
    p = sp.sympify(expr.replace('lambda','lam'), locals={**{str(x): x for x in (l,m,n,eps,E,Z)}, 'lam':lam})
    p = p.subs({Z: 2, lam: 1})
    ops.append((ast.literal_eval(off), sp.lambdify((l,m,n,eps),p.subs(E,0),'mpmath'),
                sp.lambdify((l,m,n,eps),sp.diff(p,E),'mpmath')))

def indices(N):
    return [(i,j,k-i-j) for k in range(N+1) for i in range(k+1)
            for j in range(i,k-i+1)]

def pencil(N):
    idx = indices(N); pos = {x:i for i,x in enumerate(idx)}
    H = mp.matrix(len(idx)); G = mp.matrix(len(idx))
    for row,(i,j,k) in enumerate(idx):
        mult = 1 if i==j else 2
        for off,f0,f1 in ops:
            a,b,c = i+off[0],j+off[1],k+off[2]
            if min(a,b,c)<0: continue
            t=(min(a,b),max(a,b),c)
            if t in pos:
                col=pos[t]
                H[row,col] -= mult*f0(i,j,k,EPS)
                G[row,col] += mult*f1(i,j,k,EPS)
    assert mp.norm(H-H.T,mp.inf)<mp.mpf('1e-35')
    assert mp.norm(G-G.T,mp.inf)<mp.mpf('1e-35')
    mp.cholesky(G); mp.cholesky(H+5*G)
    # Three shells form a block; rectangular off-diagonal couplings allowed.
    groups=[]
    for start in range(0,N+1,3):
        groups.append([a for a,t in enumerate(idx) if start<=sum(t)<start+3])
    for i,t in enumerate(idx):
        for j,s in enumerate(idx):
            if abs(sum(t)-sum(s))>3:
                assert H[i,j]==0 and G[i,j]==0
    return idx,H,G,groups

def sub(A,rows,cols):
    return mp.matrix([[A[i,j] for j in cols] for i in rows])

def schur_solve(A,b,groups):
    # Keep only tiny finite matrices. No off-diagonal inverse.
    S=[sub(A,g,g) for g in groups]
    d=[sub(b,g,[0]) for g in groups]
    U=[sub(A,groups[i],groups[i+1]) for i in range(len(groups)-1)]
    V=[sub(A,groups[i+1],groups[i]) for i in range(len(groups)-1)]
    for i in range(len(groups)-2,-1,-1):
        S[i] -= U[i]*(S[i+1]**-1)*V[i]
        d[i] -= U[i]*mp.lu_solve(S[i+1],d[i+1])
    x=[mp.lu_solve(S[0],d[0])]
    for i in range(1,len(groups)):
        x.append(mp.lu_solve(S[i],d[i]-V[i-1]*x[-1]))
    out=mp.matrix(A.rows,1)
    for g,v in zip(groups,x):
        for i,j in enumerate(g): out[j]=v[i]
    return out

def fmt(x): return mp.nstr(x,14)

print('input_sha256',hashlib.sha256(source.read_bytes()).hexdigest())
print('precision',mp.mp.dps,'cutoffs 2,4,6; largest dimension 50; serial')
for N in (2,4,6):
    idx,H,G,groups=pencil(N)
    e=mp.matrix(len(idx),1); e[0]=1
    for z in (mp.mpf(-5),E0+mp.mpc(0,'.05')):
        A=H-z*G
        x=schur_solve(A,e,groups)
        residual=mp.norm(A*x-e,mp.inf)
        direct=mp.lu_solve(A,e)
        agreement=mp.norm(x-direct,mp.inf)
        assert residual<mp.mpf('1e-32') and agreement<mp.mpf('1e-32')
        u=x/x[0]
        tail=mp.lu_solve(A[1:,1:],-A[1:,0:1])
        assert mp.norm(u[1:,0:1]-tail,mp.inf)<mp.mpf('1e-32')
        print('N',N,'z',fmt(z),'Schur_residual',fmt(residual),
              'direct_difference',fmt(agreement),'normalized_source',fmt(1/x[0]))
    A=H-E0*G
    mp.cholesky(A[1:,1:])
    u=mp.matrix(len(idx),1);u[0]=1
    tail=mp.lu_solve(A[1:,1:],-A[1:,0:1])
    for i in range(1,len(idx)):u[i]=tail[i-1]
    defect=(u.T*A*u)[0]
    source0=(A*u)[0]
    assert defect>0 and abs(defect-source0)<mp.mpf('1e-32')
    assert mp.norm((A*u)[1:,0:1],mp.inf)<mp.mpf('1e-32')
    if N>2: assert defect<previous
    previous=defect
    print('N',N,'size',len(idx),'E0_tail_positive definite; defect',fmt(defect),
          'raw_P_dropped_row',fmt(-source0))
print('PASS: Hermiticity, Gram/shifted/tail positivity, bandwidth, backward Schur,')
print('normalized-tail equivalence, all retained tail equations, decreasing E0 defect.')
print('Finite checks support algebra only; the infinite-limit proof is in goal1_theorem.md.')

peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print('peak_RSS_MiB',round(peak/(1024**2 if sys.platform=='darwin' else 1024),2))
