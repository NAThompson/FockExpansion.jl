# Run with Sage; ore_algebra must be on PYTHONPATH. All exact work over QQ(u).
from sage.all import *
from ore_algebra import OreAlgebra
R=PolynomialRing(QQ,'u'); u=R.gen(); K=R.fraction_field()
OA=OreAlgebra(K,'Du'); D=OA.gen()
E=4*u**2*(u**2-1)*D**2+4*u*(3*u**2-1)*D+4*u**2+1
Et=4*(u**2-1)*D**2+8*u*D+1

def primitive(L):
    cs=list(L); den=lcm([c.denominator() for c in cs]); ps=[R(c*den) for c in cs]
    g=gcd(ps); ps=[p//g for p in ps]
    return OA(ps)

def companion(L):
    n=L.order(); B=matrix(K,n)
    for i in range(n-1): B[i,i+1]=1
    for j in range(n): B[n-1,j]=-L[j]/L[n]
    return B

def rational_horizontal(B):
    """Complete rational horizontal space, using a cyclic observable and rational_solutions.
    Failure to find an observable raises, never returns a false negative.
    """
    n=B.nrows(); C=None
    candidates=[vector(K,[1 if i==j else 0 for i in range(n)]) for j in range(n)]
    candidates += [vector(K,[(u+c)**j for j in range(n)]) for c in range(1,9)]
    for w in candidates:
        rows=[w]
        for j in range(n): rows.append(rows[-1].derivative(u)+rows[-1]*B)
        C=matrix(K,rows[:n])
        if C.det()!=0: break
    else: raise RuntimeError('No cyclic observable found; undecided')
    rel=rows[n]*C.inverse(); S=primitive(D**n-sum(rel[j]*D**j for j in range(n)))
    sols=S.rational_solutions(); result=[]
    for sol in sols:
        q=K(sol[0]); v=C.solve_right(vector(K,[q.derivative(u,j) for j in range(n)]))
        assert v.derivative(u)==B*v
        result.append(v)
    return result

def gauge(F,T):
    """Rational gauge H mapping solution columns of T to those of F."""
    if F.order()!=2: return None
    A=companion(T); B=companion(F)
    C=matrix(K,4)
    for k in range(4):
        H=matrix(K,2,[1 if i==k else 0 for i in range(4)])
        col=B*H-H*A
        C.set_column(k,vector(K,col.list()))
    basis=rational_horizontal(C)
    # For irreducible rank-two T every nonzero Hom is invertible (Schur).
    for v in basis:
        H=matrix(K,2,list(v))
        if H.det()!=0:
            assert H.derivative(u)==B*H-H*A
            return H
    return None

def ftest(L, factors=None):
    """Factor and test all Jordan-Hoelder constituents over QQ(u).
    A supplied factor list must have exact product L; use it to avoid costly
    seminumeric factorization when an independently proved composition series
    is known. Exceptions/constant extensions mean UNDECIDED, not absence.
    This recognizes rational gauge equivalence, including isogeny-induced
    rational de Rham maps; it does not cover arbitrary algebraic pullbacks.
    """
    L=OA(L)
    complete = factors is None
    if factors is None: factors=L.factor()
    factors=[OA(F) for F in factors]
    assert prod(factors)==L
    hits=[]
    for i,F in enumerate(factors):
        if F.order()!=2: continue
        for name,T in [('E',E),('sqrt_u_twist',Et)]:
            H=gauge(F,T)
            if H is not None: hits.append((i,name,H))
    return {'hits':hits,'factors':factors,
            'status':('PRESENT' if hits else ('ABSENT' if complete or all(F.order()<=2 for F in factors) else 'UNDECIDED: supplied factors need refinement'))}

if __name__=='__main__':
    print('Use load(...x29_ftest.sage), then ftest(L). Exceptions mean UNDECIDED.')
