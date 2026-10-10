"""PROOF (computer algebra, exact) of the Clausen/log closed forms of the conic integrals
   pi*Ktilde_{s0}(b,y) = Phi_{s0}(alpha,beta),  b = cos(alpha), y = sin(beta),
   Ktilde_inf = sqrt(1-b) K_inf,  Ktilde_1 = 2 sqrt(1-y^2) K_1,  Ktilde_{-1} = 2|y| K_{-1},
   K_{s0} = Reg int_0^1 Ahat(t) t dt/((t^2-s0) sqrt(P(t)))   (K_inf: no factor 1/(t^2-s0)),
on the diamond Dm = {|alpha-pi/2| + |2beta+pi/2| < pi/2} = {(b,y): b^2+4y^4-4y^2 < 0, y<0}.
Step 1: exact gradient of Ktilde (Hermite reduction on the conic + partial fractions: the integrands reduce to
        1/(t^2+1), 1/((1-b)t^2+1+b) and, for s0=1, 1/(t+-1)).  Step 2: exact gradient of Phi.  Step 3: compare
        exactly (trig -> exponentials, reduction modulo the cyclotomic polynomial).  Step 4: value at the vertex (b,y)=(0,0).
Run: PYTHONPATH=.../python-deps python3 conic_proof.py"""
import sympy as sp, mpmath as mp
t,b,y,s=sp.symbols('t b y s'); al,be=sp.symbols('alpha beta',real=True)
AL,BE,LOG2,PI=sp.symbols('AL BE LOG2 PI')     # free occurrences of alpha, beta, log2, pi (outside trig)
pi=sp.pi
# ---------------- Step 1: integral side -----------------
P=(1+b)+(2-4*y**2)*t**2+(1-b)*t**4; p=(1+b)+(2-4*y**2)*s+(1-b)*s**2
Nn=(1+t**2)*(1+b+(1-b)*t**2)
dA={'t':2*y*(1+b-(1-b)*t**4)/Nn,'y':2*t,'b':-y*t*(1-t**2)/(1+b+(1-b)*t**2)}  # dA/dx * sqrt(P), A = asin(2yt/sqrt N)
# check the three derivative formulas of A exactly
Aexpr=sp.asin(2*y*t/sp.sqrt(Nn))
for x,X in (('t',t),('y',y),('b',b)):
    e=sp.diff(Aexpr,X)-dA[x]/sp.sqrt(P)
    # A' formulas verified by squaring-free substitution: 1-4y^2t^2/N = P/N
    e=sp.simplify(e.subs(sp.sqrt(1-4*y**2*t**2/Nn),sp.sqrt(P/Nn)))
    vals=[e.subs({t:sp.Rational(1,3),b:sp.Rational(1,5),y:-sp.Rational(1,2)}).evalf(30)]
    assert abs(vals[0])<1e-25, (x,vals)
def nfac(s0):
    return sp.sqrt(1-b) if s0=='inf' else sp.sqrt(p.subs(s,s0))
def hermite(s0,x):
    r0,r1=sp.symbols('r0 r1'); R=r0+r1*s; n=nfac(s0)
    g=n/sp.sqrt(p) if s0=='inf' else n/((s-s0)*sp.sqrt(p))
    lhs=sp.diff(g,x); rhs=sp.diff(n*R/sp.sqrt(p),s)
    e=sp.together(sp.simplify((lhs-rhs)*p**sp.Rational(3,2)/n*(1 if s0=='inf' else (s-s0)**2)))
    sol=sp.solve(sp.Poly(sp.expand(sp.numer(e)),s).coeffs(),[r0,r1],dict=True)
    R=sp.factor(R.subs(sol[0]))
    chk=sp.simplify(sp.diff(g,x)-sp.diff(n*R/sp.sqrt(p),s)); assert chk==0, chk
    return n,R
def grad_K(s0,x):
    """exact d/dx Ktilde as a polynomial in AL (=alpha), BE (=beta), LOG2, with coefficients algebraic in b,y."""
    X=sp.Symbol(x); n,R=hermite(s0,X)
    den=1 if s0=='inf' else (t**2-s0)
    rho=sp.cancel((2*dA[x]*t/(den*P)-R.subs(s,t**2)*dA['t']/P))    # = (pi/n) * integrand
    ap=sp.apart(rho,t); tot=0
    for term in sp.Add.make_args(ap):
        num,dn=sp.fraction(sp.factor(term)); dn=sp.Poly(dn,t)
        c=sp.cancel(num/dn.LC()); m=sp.monic(dn.as_expr(),t)
        if sp.expand(m-(t**2+1))==0: tot+=c*PI/4                                   # int_0^1 dt/(1+t^2)
        elif sp.expand(m-(t**2+(1+b)/(1-b)))==0: tot+=c*AL/(2*sp.sin(al))*(1-b)       # int dt/(t^2+rho^2) = (alpha/2)/(rho) ; (1-b)rho = sin(alpha)
        elif sp.expand(m-(t+1))==0: tot+=c*LOG2                                     # int_0^1 dt/(t+1)
        elif sp.expand(m-(t-1))==0: tot+=0                                          # Reg int_0^1 dt/(t-1) = 0
        else: raise ValueError(term)
    # boundary term (1/2)[Ahat n R/sqrt P]_0^1 * pi :  Ahat(1) = 1+2beta/pi, Ahat(0)=1, sqrtP(1)=2sqrt(1-y^2), sqrtP(0)=sqrt(1+b)
    Rs=R.subs(s,1); R0=R.subs(s,0)
    bnd=sp.Rational(1,2)*n*((PI+2*BE)*Rs/(2*sp.sqrt(1-y**2))-PI*R0/sp.sqrt(1+b))
    return sp.expand(n*tot/PI*PI + bnd)            # this is pi * d_x Ktilde  (tot already carries the pi of 2/pi*... see rho)
# NOTE on normalisation: d_x Ktilde = (n/pi) int rho dt + (1/2)[Ahat n R/sqrtP]_0^1; we return pi*d_x Ktilde.
def grad_K_pi(s0,x):
    X=sp.Symbol(x); n,R=hermite(s0,X)
    den=1 if s0=='inf' else (t**2-s0)
    rho=sp.cancel((2*dA[x]*t/(den*P)-R.subs(s,t**2)*dA['t']/P))
    ap=sp.apart(rho,t); tot=0
    for term in sp.Add.make_args(ap):
        num,dn=sp.fraction(sp.factor(term))
        dt_=sp.Mul(*[f for f in sp.Mul.make_args(dn) if f.has(t)]); rest=sp.cancel(dn/dt_)
        lc=sp.Poly(dt_,t).LC(); m=sp.expand(dt_/lc); c=sp.cancel(num/(lc*rest))
        if sp.cancel(m-(t**2+1))==0: tot+=c*PI/4
        elif sp.cancel(m-(t**2+(1+b)/(1-b)))==0: tot+=c*AL*(1-b)/(2*sp.sin(al))
        elif sp.expand(m-(t+1))==0: tot+=c*LOG2
        elif sp.expand(m-(t-1))==0: pass
        else: raise ValueError((term,m))
    Rs=R.subs(s,1); R0=R.subs(s,0)
    bnd=sp.Rational(1,2)*n*((PI+2*BE)*Rs/(2*sp.sqrt(1-y**2))-PI*R0/sp.sqrt(1+b))
    return n*tot+bnd
# ---------------- Step 2: closed forms Phi -----------------
class Cl2(sp.Function):
    def fdiff(self,argindex=1): return -sp.log(2*sp.sin(self.args[0]/2))
g=[pi/4+al/2,pi/4-al/2,-pi/4+al/2,-pi/4-al/2]
ls=[sp.log(sp.sin((be-x)/2)) for x in g]; lc=[sp.log(sp.cos((be-x)/2)) for x in g]
Cm=[Cl2(be-x) for x in g]; Cp=[Cl2(be-x+pi) for x in g]
def Phi(s0):
    if s0=='inf':
        S=[-1,-1,1,1]; u=[1,-1,-1,1]
        v=sum(al/4*S[k]*(ls[k]-lc[k])+be/2*u[k]*(ls[k]-lc[k])+u[k]*(Cm[k]-Cp[k])/2 for k in range(4))
        v+=pi/8*(-ls[0]-3*lc[0]+ls[1]+3*lc[1]-ls[2]+5*lc[2]+ls[3]-5*lc[3]); return v
    if s0==1:
        S=[-1,1,-1,1]
        v=sum(al/4*S[k]*(ls[k]+lc[k])+be/2*(ls[k]+lc[k])+(Cm[k]+Cp[k])/2 for k in range(4))
        v+=-2*be*sp.log(sp.cos(be))-Cl2(2*be+pi)
        v+=pi/8*(-ls[0]+3*lc[0]-ls[1]+3*lc[1]+ls[2]+5*lc[2]+ls[3]+5*lc[3])-pi*sp.log(sp.cos(be))
        v+=sp.log(2)*be+pi*sp.log(2)/2; return v
    S=[1,-1,-1,1]; u=[-1,-1,1,1]
    # cos(al/2)+sin(al/2) = sqrt2 sin(al/2+pi/4),  cos(al/2)-sin(al/2) = sqrt2 sin(pi/4-al/2)  (identical; written so that logs canonicalise)
    cs=sp.log(sp.sqrt(2)*sp.sin(al/2+pi/4)); dd=sp.sqrt(2)*sp.sin(pi/4-al/2)
    v=sum(al/4*S[k]*(ls[k]+lc[k])+be/2*u[k]*(ls[k]+lc[k])+u[k]*(Cm[k]+Cp[k])/2 for k in range(4))
    v+=(al/2-pi/4)*sp.log(dd)-(al/2+pi/4)*cs-Cl2(al+pi/2)/2+Cl2(al-pi/2)/2
    v+=pi/8*(ls[0]-3*lc[0]+ls[1]-3*lc[1]+ls[2]+5*lc[2]+ls[3]+5*lc[3])+3*pi*sp.log(2)/4
    return v
# ---------------- Step 3: exact comparison -----------------
a8,b4,X,Y,zt=sp.symbols('a8 b4 X Y zeta')      # alpha = 8 a8, beta = 4 b4, X = e^{i a8}, Y = e^{i b4}, zeta = e^{i pi/16}
def free_to_symbols(e):
    """replace alpha, beta, log 2, pi occurring outside trig/log-of-trig by AL, BE, LOG2, PI."""
    e=sp.expand(sp.expand_log(e,force=True))
    out=0
    for term in sp.Add.make_args(e):
        fac=sp.Mul.make_args(term); new=1
        for f in fac:
            if f==al: new*=AL
            elif f==be: new*=BE
            elif f==pi: new*=PI
            elif f.is_Pow and f.base in (al,be,pi) and f.exp.is_Integer:
                new*={al:AL,be:BE,pi:PI}[f.base]**f.exp
            elif f==sp.log(2): new*=LOG2
            else: new*=f
        out+=new
    return out
def to_exp_poly(e):
    """trig expression (no free alpha/beta) -> rational function in X, Y, zeta."""
    e=e.subs({al:8*a8,be:4*b4})
    e=e.rewrite(sp.exp)
    e=sp.expand(e)
    def rep(ex):
        arg=sp.expand(ex.args[0]/sp.I)
        ca=arg.coeff(a8); cb=arg.coeff(b4); cp=sp.expand(arg-ca*a8-cb*b4)/pi
        assert sp.expand(arg-ca*a8-cb*b4-cp*pi)==0, ex
        k=cp*16; assert k.is_Integer and ca.is_Integer and cb.is_Integer, ex
        return X**int(ca)*Y**int(cb)*zt**int(k)
    e=e.replace(lambda ex: isinstance(ex,sp.exp), rep)
    e=e.subs(sp.sqrt(2),zt**4+zt**-4)
    return e
A4,B2,sa,ca,sb,cb,c8,r2s=sp.symbols('A4 B2 sa ca sb cb c8 r2s')
GB=sp.groebner([sa**2+ca**2-1,sb**2+cb**2-1,r2s**2-2,4*c8**2-2-r2s],r2s,c8,sa,ca,sb,cb,order='lex')
def is_zero_trig(e):
    """exact zero test: alpha = 4 A4, beta = 2 B2; sin/cos(A4) = sa, ca; sin/cos(B2) = sb, cb; cos(pi/8) = c8, sqrt2 = r2s.
    Clear denominators and reduce the numerator modulo the prime ideal
    (sa^2+ca^2-1, sb^2+cb^2-1, r2s^2-2, 4c8^2-2-r2s); remainder 0 proves the identity (denominators are nonzero on Dm)."""
    e=sp.expand_trig(sp.sympify(e).subs({al:4*A4,be:2*B2}))
    # |f| -> +-f with the sign of f on the diamond (f has constant sign there; test point alpha=pi/2, beta=-pi/4)
    def sgn_abs(ex):
        f=ex.args[0]; v=sp.N(f.subs({A4:sp.pi/8,B2:-sp.pi/8}),30)
        assert abs(v)>1e-10, ex
        return f if v>0 else -f
    e=e.replace(lambda ex: isinstance(ex,sp.Abs), sgn_abs)
    e=e.subs({sp.sqrt(2+sp.sqrt(2)):2*c8, sp.sqrt(2-sp.sqrt(2)):sp.sqrt(2)/(2*c8)})
    e=e.subs(sp.sqrt(2),r2s)
    def rad(ex):
        if ex.is_Pow and ex.exp.is_Rational and ex.exp.q==2:
            base=sp.expand(ex.base)
            if base==r2s+2: v=2*c8
            elif base==2-r2s: v=r2s/(2*c8)
            elif sp.expand(4*base)==r2s+2: v=c8                 # cos(pi/8)
            elif sp.expand(4*base)==2-r2s: v=r2s/(4*c8)         # sin(pi/8) = sqrt2/(4 cos(pi/8))
            else: return ex
            return v**int(2*ex.exp)
        return ex
    e=e.replace(lambda ex: ex.is_Pow, rad)
    e=e.subs({sp.sin(A4):sa,sp.cos(A4):ca,sp.sin(B2):sb,sp.cos(B2):cb})
    assert not e.has(sp.sin) and not e.has(sp.cos) and not e.has(A4) and not e.has(B2), e.atoms(sp.Function)
    num=sp.numer(sp.together(e))
    r=GB.reduce(sp.expand(num))[1]
    return sp.expand(r)==0
def canon_logs(e):
    """log(sin L), log(cos L) -> log(sin L') with L' canonical (log|.| semantics: sign and period pi irrelevant)."""
    def can(L):
        L=sp.expand(L)
        lead=[L.coeff(al),L.coeff(be)]
        if (lead[0]<0) or (lead[0]==0 and lead[1]<0): L=-L
        c=sp.expand(L-L.coeff(al)*al-L.coeff(be)*be)/pi
        c=c-sp.floor(c)
        return sp.expand(L.coeff(al)*al+L.coeff(be)*be+c*pi)
    def rep(ex):
        f=ex.args[0]
        if isinstance(f,sp.sin): return sp.log(sp.sin(can(f.args[0]),evaluate=False))
        if isinstance(f,sp.cos): return sp.log(sp.sin(can(sp.pi/2-f.args[0]),evaluate=False))
        return ex
    return e.replace(lambda ex: isinstance(ex,sp.log), rep)
def check(s0):
    F=Phi(s0)
    # log|.| semantics: expand_log turns log of a negative factor into log(-1) = i pi; these constants are dropped (I -> 0)
    dPa=sp.expand(sp.expand_log(sp.diff(F,al),force=True)).subs(sp.I,0); dPb=sp.expand(sp.expand_log(sp.diff(F,be),force=True)).subs(sp.I,0)
    sub={b:sp.cos(al),y:sp.sin(be)}
    sq={sp.sqrt(1-b):sp.sqrt(2)*sp.sin(al/2),sp.sqrt(1+b):sp.sqrt(2)*sp.cos(al/2),sp.sqrt(1-y**2):sp.cos(be)}
    res={}
    for x,dP,jac in (('b',dPa,-1/sp.sin(al)),('y',dPb,1/sp.cos(be))):
        G=grad_K_pi(s0,x)
        G=G.subs(sp.sqrt(p.subs(s,s0)) if s0!='inf' else sp.sqrt(1-b), nfac(s0))
        # algebraic prefactors in terms of alpha, beta on the diamond (alpha in (0,pi), beta in (-pi/2,0))
        G=G.subs({sp.sqrt(4-4*y**2):2*sp.sqrt(1-y**2),sp.sqrt(4*y**2):-2*y,sp.sqrt(y**2):-y})   # y<0 on the diamond
        G=G.subs(sq).subs(sub)
        G=G.subs({sp.sqrt(2-2*sp.cos(al)):2*sp.sin(al/2),sp.sqrt(1-sp.cos(al)):sp.sqrt(2)*sp.sin(al/2),
                  sp.sqrt(1+sp.cos(al)):sp.sqrt(2)*sp.cos(al/2),sp.sqrt(1-sp.sin(be)**2):sp.cos(be),
                  sp.sqrt(4-4*sp.sin(be)**2):2*sp.cos(be),sp.sqrt(4*sp.sin(be)**2):-2*sp.sin(be)})
        D=free_to_symbols(G)-free_to_symbols(sp.expand(jac*dP))
        D=sp.expand(canon_logs(D))
        # remaining logs of trig functions must cancel
        logs=[f for f in D.atoms(sp.log) if f!=sp.log(2)]
        D=sp.collect(D,[AL,BE,LOG2,PI]+logs,evaluate=False)
        ok=True
        for k,c in D.items():
            import time; t0=time.time(); z=is_zero_trig(c); print('      coeff',k,z,round(time.time()-t0,1),flush=True)
            ok&=z
            if not z: print('   nonzero coefficient of',k)
        res[x]=ok; print('s0=%s  d/d%s:  exact identity %s   (log atoms after cancellation: %d)'%(s0,x,ok,len(logs)))
    return res
if __name__=='__main__':
    import sys
    for s0 in (sys.argv[1:] and [eval(a) if a!='inf' else a for a in sys.argv[1:]] or ['inf',1,-1]): check(s0)
