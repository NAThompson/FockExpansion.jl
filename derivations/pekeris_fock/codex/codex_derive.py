#!/usr/bin/env python3
"""Independent distance-coordinate derivation. No input from sym/.

Equation-index convention: coefficient at (l,m,n) multiplies A(l+dl,m+dm,n+dn).
We use exact finite differences of the standard Laguerre identities, including
all derivative tails through a certified polynomial tail region.
"""
import sympy as s
import itertools, json
from pathlib import Path
from functools import lru_cache
u,v,w=s.symbols('u v w'); x=(u,v,w)
a,b,c=s.symbols('r1 r2 r12', positive=True)
eps,E,Z,lam=s.symbols('eps E Z lambda')
l,m,n=s.symbols('l m n', integer=True, nonnegative=True); indices=(l,m,n)
OUT=Path(__file__).resolve().parent

def derive():
    # Delta_1+Delta_2 in independent interparticle distances, symmetric matrix.
    G=s.Matrix([[1,0,(a*a+c*c-b*b)/(2*a*c)],
                [0,1,(b*b+c*c-a*a)/(2*b*c)],
                [(a*a+c*c-b*b)/(2*a*c),(b*b+c*c-a*a)/(2*b*c),2]])
    h=s.Matrix([2/a,2/b,4/c])
    J=s.Matrix([[-eps,eps,eps],[eps,-eps,eps],[2*eps,2*eps,-2*eps]])
    inv={a:(2*v+w)/(4*eps),b:(2*u+w)/(4*eps),c:(u+v)/(2*eps)}
    Q=(2*v+w)*(2*u+w)*(u+v)
    M=(Q*(J*G*J.T)/2).applyfunc(lambda z:s.cancel(z.subs(inv)).expand())
    drift=(Q*J*h/2).applyfunc(lambda z:s.cancel(z.subs(inv)).expand())
    pot=s.cancel((Q*(Z/a+Z/b-lam/c+E)).subs(inv)).expand()
    # e^(sigma/2) D_i e^(-sigma/2)=D_i-1/2.
    B=[s.expand(drift[i]-sum(M[i,j] for j in range(3))) for i in range(3)]
    C=s.expand(pot-sum(drift)/2+sum(M)/4)
    P={}
    def add(t,d,coef):
        t[d]=s.expand(t.get(d,0)+coef)
    for i in range(3):
        for j in range(3):
            d=[0]*3;d[i]+=1;d[j]+=1;add(P,tuple(d),M[i,j])
        d=[0]*3;d[i]=1;add(P,tuple(d),B[i])
    add(P,(0,0,0),C)
    sig=u+v+w; T=s.expand(sum(M))
    C1={(0,0,0):s.expand(sig*sum(B)-T)}
    for j in range(3):
        d=[0]*3;d[j]=1;C1[tuple(d)]=s.expand(2*sig*sum(M[i,j] for i in range(3)))
    ops={'P':P,'S2P':{d:s.expand(sig**2*f) for d,f in P.items()},'C1':C1,'C2':{(0,0,0):T}}
    assert s.expand(T-4*eps**2*Q)==0
    return ops

@lru_cache(None)
def multiplication(power):
    # L_k -> sum_h weight(k) L_(k+h).
    k=s.Symbol('k'); out={0:s.Integer(1)}
    for _ in range(power):
        nxt={}
        for h,f in out.items():
            for shift,g in ((1,-(k+h+1)),(0,2*(k+h)+1),(-1,-(k+h))):
                nxt[h+shift]=s.expand(nxt.get(h+shift,0)+f*g)
        out=nxt
    return out

@lru_cache(None)
def oned(power,deriv,index,bound):
    # x^power D^deriv L_source projected at equation index.
    k=s.Symbol('k'); result={}
    for d in range(-power,bound+1):
        val=0
        for h,f in multiplication(power).items():
            # derivative output k=index-h, source=index+d.
            gap=d+h
            if deriv==0:
                weight=1 if gap==0 else 0
            elif gap>=deriv:
                weight=(-1)**deriv*s.binomial(gap-1,deriv-1)
            else:
                weight=0
            val+=f.subs(k,index-h)*weight
        val=s.expand(val)
        if val!=0:result[d]=val
    return result

def recurrence(op,bound=8):
    acc={}
    for deriv,f in op.items():
        for powers,coef in s.Poly(f,*x).terms():
            maps=[oned(powers[i],deriv[i],indices[i],bound) for i in range(3)]
            for entries in itertools.product(*(list(t.items()) for t in maps)):
                d=tuple(z[0] for z in entries)
                val=coef*s.prod(z[1] for z in entries)
                acc[d]=acc.get(d,0)+val
    result={}
    for d,f in acc.items():
        f=s.expand(f)
        if f!=0: result[d]=s.factor(f)
    # For d_i > max coordinate degree, derivative tails are polynomial of degree
    # <=1 (operator order <=2). Two consecutive zero slices certify each tail.
    maxdegree=max(s.Poly(f,*x).degree(z) for f in op.values() for z in x)
    assert maxdegree+2<=bound
    assert not any(max(d)>maxdegree for d in result), 'nonzero derivative tail'
    return result

def apply(op,f):
    return s.Poly(s.expand(sum(coef*s.diff(f,u,d[0],v,d[1],w,d[2]) for d,coef in op.items())),*x)

def laguerre_coeffs(poly):
    # Exact monomial projection: integral exp(-x)x^p L_i = (-1)^i p! binomial(p,i).
    out={}
    for powers,coef in poly.terms():
        for target in itertools.product(*(range(p+1) for p in powers)):
            q=coef*s.prod((-1)**j*s.factorial(p)*s.binomial(p,j) for p,j in zip(powers,target))
            out[target]=out.get(target,0)+q
    return {k:s.expand(v) for k,v in out.items() if s.expand(v)!=0}

def verify(ops,recs):
    values={eps:s.Rational(3,2),E:s.Rational(-7,3),Z:2,lam:s.Rational(4,5)}
    count=0
    for name,op in ops.items():
        numeric={d:f.subs(values) for d,f in op.items()}
        for source in [(0,0,0),(1,2,0),(3,1,2)]:
            f=s.prod(s.laguerre(j,z) for j,z in zip(source,x))
            expected=laguerre_coeffs(apply(numeric,f))
            actual={}
            for d,coef in recs[name].items():
                target=tuple(source[i]-d[i] for i in range(3))
                if min(target)<0:continue
                val=s.expand(coef.subs(dict(zip(indices,target))).subs(values))
                if val!=0:actual[target]=val
            assert expected==actual,(name,source,[(k,expected.get(k,0)-actual.get(k,0)) for k in expected.keys()|actual.keys() if expected.get(k,0)!=actual.get(k,0)][:5])
            count+=1
    # Independent product-rule check on monomials, using t=log(sigma).
    sig=sum(x); f=1+u*v+w**2
    log=s.log(sig)
    for k in (1,2,3):
        left=s.expand(sum(coef*s.diff(log**k*f,u,d[0],v,d[1],w,d[2]) for d,coef in ops['P'].items()))*sig**2
        right=log**k*apply(ops['S2P'],f).as_expr()+k*log**(k-1)*apply(ops['C1'],f).as_expr()
        if k>=2:right+=k*(k-1)*log**(k-2)*apply(ops['C2'],f).as_expr()
        assert s.cancel(s.expand(left-right))==0
    return count

def main():
    ops=derive(); recs={}
    for name,op in ops.items():
        print('Deriving',name,flush=True)
        recs[name]=recurrence(op)
        print(name,len(recs[name]),flush=True)
    count=verify(ops,recs)
    save=lambda name,data:(OUT/name).write_text(json.dumps(data,indent=2)+'\n')
    save('codex_recurrence.json',{name:{str(d):str(f) for d,f in sorted(rec.items())} for name,rec in recs.items()})
    save('codex_operator.json',{name:{str(d):str(f) for d,f in sorted(op.items())} for name,op in ops.items()})
    stats={}
    for name,rec in recs.items():
        fixed={d:s.factor(f.subs(E,-eps**2)) for d,f in rec.items() if s.expand(f.subs(E,-eps**2))!=0}
        stats[name]={'general_count':len(rec),'general_shell_range':[min(map(sum,rec)),max(map(sum,rec))],'fixed_count':len(fixed),'fixed_shell_range':[min(map(sum,fixed)),max(map(sum,fixed))]}
    stats['verification']={'polynomial_action_cases':count,'log_product_rule_cases':3,'tail_bound':8}
    save('codex_counts.json',stats); print(json.dumps(stats,indent=2),flush=True)

if __name__=='__main__': main()
