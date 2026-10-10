"""Step 4 of the conic proof: the constant.  pi*Ktilde - Phi is constant on the diamond (conic_proof.py).  At the vertex
(b,y) = (0,0), i.e. (alpha,beta) -> (pi/2, 0^-), Ktilde is continuous (P > 0 on [0,1] near (0,0)) with
   Ktilde_inf(0,0) = int_0^1 t dt/(1+t^2) = log2/2,   Ktilde_1(0,0) = 2 Reg int_0^1 t dt/((t^2-1)(t^2+1)) = 0,   Ktilde_{-1}(0,0) = 0.
Since Phi = pi*Ktilde - const, lim Phi exists and may be computed along alpha = pi/2, beta -> 0^-."""
import sympy as sp, mpmath as mp, sys
sys.path.insert(0,'.')
from conic_proof import Phi, Cl2, al, be
t=sp.Symbol('t')
# values of Ktilde at the vertex (exact)
Kinf0=sp.integrate(t/(1+t**2),(t,0,1))
eps=sp.Symbol('eps',positive=True)
# antiderivative of t/((t^2-1)(t^2+1)) on [0,1): (1/4)log(1-t^2) - (1/4)log(1+t^2); Reg = constant term in log(eps), t = 1-eps
Fa=lambda tt: sp.log(1-tt**2)/4-sp.log(1+tt**2)/4
assert sp.simplify(sp.diff(Fa(t),t)-t/((t**2-1)*(t**2+1)))==0
up=sp.expand_log(sp.log(eps)/4+sp.log(2-eps)/4-sp.log(1+(1-eps)**2)/4,force=True)
K10=2*(sp.limit(up-sp.log(eps)/4,eps,0)-Fa(0))
print('Ktilde_inf(0,0) =',sp.simplify(Kinf0),'  Ktilde_1(0,0) =',sp.simplify(K10),'  Ktilde_-1(0,0) = 0 (factor 2|y|)')
target={'inf':sp.pi*Kinf0,1:sp.pi*K10,-1:0}
G=sp.Symbol('G')
def cl2_at_zero(arg):
    v=sp.nsimplify(arg.subs(be,0)/sp.pi)          # multiple of pi
    v=v-2*sp.floor(v/2)                            # mod 2pi
    tab={0:0,sp.Rational(1,2):G,1:0,sp.Rational(3,2):-G}
    return tab[v]
ok=True
for s0 in ['inf',1,-1]:
    # on alpha = pi/2 the term (alpha/2-pi/4) log|cos(alpha/2)-sin(alpha/2)| of Phi_{-1} is 0 (continuous extension, x log|x| -> 0)
    F=Phi(s0).subs(sp.log(sp.sqrt(2)*sp.sin(sp.pi/4-al/2)),0).subs(al,sp.pi/2)
    cl=[f for f in F.atoms(Cl2)]
    Fe=F.subs({c:cl2_at_zero(c.args[0]) for c in cl})        # Cl2 is continuous
    e_=sp.Symbol('e_',positive=True); Fe=Fe.subs(be,-e_)
    def absl(ex):       # log -> log|.| : flip sign of the argument where it is negative for small e_>0
        f=ex.args[0]; v=f.subs(e_,sp.Rational(1,10**6)).evalf(30)
        return sp.log(-f) if v<0 else ex
    Fe=Fe.replace(lambda ex: isinstance(ex,sp.log), absl)
    lim=sp.limit(Fe,e_,0,'+')
    lim=sp.nsimplify(sp.simplify(lim))
    d=sp.simplify(lim-target[s0]); ok&=(d==0)
    print('s0=%s: lim Phi = %s,  pi*Ktilde(vertex) = %s,  difference %s'%(s0,lim,target[s0],d))
    # numerical sanity: Phi at a point near the vertex inside the diamond
    mp.mp.dps=30
print('VERTEX CONSTANTS VANISH' if ok else 'MISMATCH')
