"""Small serial exact checks for the phase-11 mathematical audit.

No Claude code is imported. The displayed L5 is transcribed from COLLAB.
All writes use this file's directory. No numerical fit is a proof here.
"""
import sys
sys.dont_write_bytecode = True
import json, math, resource, time
from pathlib import Path
import sympy as s
import mpmath as mp

start = time.time()
out = Path(__file__).resolve().parent
checks, facts = {}, {}
u, v, k, a, b, q, r, z, ell, rho = s.symbols('u v k a b q r z ell rho')
zero = lambda f: s.cancel(s.expand(f)) == 0

# Scalar differential operators, with coefficient on the left.
def add(*ops):
    ans = {}
    for op in ops:
        for ij, c in op.items():
            ans[ij] = ans.get(ij, 0) + c
    return {ij: s.expand(c) for ij, c in ans.items() if s.expand(c) != 0}

def scale(c, op):
    return {ij: s.expand(c*f) for ij, f in op.items() if c*f != 0}

def mul(A, B):
    ans = {}
    for (i,j), f in A.items():
        for (m,n), g in B.items():
            for di in range(i+1):
                for dj in range(j+1):
                    ij = (m+i-di, n+j-dj)
                    c = f*s.binomial(i,di)*s.binomial(j,dj)*s.diff(g,u,di,v,dj)
                    ans[ij] = ans.get(ij, 0) + c
    return add(ans)

lap = {(2,0):1-u*u, (1,1):-2*u*v, (0,2):1-v*v,
       (1,0):-3*u, (0,1):-3*v}
L = lambda kk: add(lap, {(0,0):kk*(kk+2)})
Tu = lambda kk: {(1,0):1-u*u,(0,1):-u*v,(0,0):kk*u}
Tv = lambda kk: {(1,0):-u*v,(0,1):1-v*v,(0,0):kk*v}
for tag, T in [('u',Tu),('v',Tv)]:
    checks['lowering_'+tag] = not add(mul(L(k-1),T(k)), scale(-1,mul(T(k-2),L(k))))
R = mul(Tv(-3),mul(Tv(-2),Tv(-1)))
S = mul(Tv(-5),mul(Tv(-4),Tv(-3)))
checks['three_step_P_harmonic_intertwiner'] = not add(mul(L(-4),R),scale(-1,mul(S,L(-1))))
facts['Rv3'] = {str(ij):str(s.factor(c)) for ij,c in sorted(R.items())}

# For any fixed vector (a,b,0,0), Q=|a|^2-(a.x)^2.
Q = a*a+b*b-(a*u+b*v)**2
logdu = -s.diff(Q,u)/(2*Q)
logdv = -s.diff(Q,v)/(2*Q)
d2u = -s.diff(Q,u,2)/(2*Q)+3*s.diff(Q,u)**2/(4*Q**2)
d2v = -s.diff(Q,v,2)/(2*Q)+3*s.diff(Q,v)**2/(4*Q**2)
dudv = -s.diff(Q,u,v)/(2*Q)+3*s.diff(Q,u)*s.diff(Q,v)/(4*Q**2)
checks['K3_measure_Delta_minus_one'] = zero((1-u*u)*d2u-2*u*v*dudv+(1-v*v)*d2v-3*u*logdu-3*v*logdv-1)
checks['lowered_section_nonzero_leading_coefficient'] = s.prod(s.Rational(1,2)-j for j in range(3)) == s.Rational(3,8)
facts['Rv3_half_branch_leading_coefficient_at_u1'] = '-3*v**3/8 times (u-1)**(-5/2) times original leading coefficient'

# The supplied order-five scalar operator and its indicial polynomials.
cs = [8505*ell,144585*(ell**2-s.Rational(73,357)),
      272160*ell*(ell**2-s.Rational(107,189)),
      140616*(ell**4-s.Rational(2171,1953)*ell**2+s.Rational(344,1953)),
      24624*ell*(ell**2-s.Rational(164,171))*(ell**2-s.Rational(8,9)),
      1296*(ell**2-1)*(ell**2-s.Rational(8,9))**2]
fall = lambda rr,j: s.prod(rr-i for i in range(j))
def indicial(coeffs, var, at):
    terms=[]
    for j,c in enumerate(coeffs):
        if c == 0: continue
        poly=s.Poly(s.expand(c.subs(var,at+z)),z)
        val=min(mm[0] for mm,cc in poly.terms() if cc != 0)
        terms.append((val-j,poly.nth(val)*fall(rho,j)))
    mn=min(i for i,c in terms)
    return s.factor(sum(c for i,c in terms if i == mn))
expected_boundary=rho**2*(rho-1)**2*(rho-2)
expected_unit=rho*(rho-s.Rational(1,2))*(rho-1)*(rho-2)*(rho-3)
for at,expected,name in [(s.sqrt(8)/3,expected_boundary,'plus_boundary'),
                         (-s.sqrt(8)/3,expected_boundary,'minus_boundary'),
                         (s.Integer(1),expected_unit,'plus_one'),
                         (-s.Integer(1),expected_unit,'minus_one')]:
    ip=indicial(cs,ell,at)
    checks['indicial_'+name] = s.Poly(ip,rho).monic().as_expr() == s.expand(expected)
    facts['indicial_'+name] = str(ip)

# Replace ell by sqrt(z), f(ell)=g(ell^2); Lf/ell has rational coefficients.
# Chain-rule derivative coefficients obtained independently by recurrence.
ders=[{0:s.Integer(1)}]
for j in range(5):
    dd={}
    for n,c in ders[-1].items():
        dd[n]=dd.get(n,0)+s.diff(c,ell)
        dd[n+1]=dd.get(n+1,0)+2*ell*c
    ders.append(dd)
zc=[]
for n in range(6):
    e=s.expand(sum(cs[j]*ders[j].get(n,0) for j in range(6))/ell)
    zc.append(s.expand(e.subs(ell**2,z)))
checks['even_pullback_rational'] = all(not c.has(ell) for c in zc)
t=s.symbols('t'); z0=s.Rational(8,9)
pc=[s.Poly(c.subs(z,z0+t),t) for c in zc]
# Lowest Taylor equation t^n has a_(n+3) coefficient ip(n+3).
ip=indicial(zc,z,z0)
checks['analytic_recurrence_nonzero_all_orders'] = s.Poly(ip,rho).monic().as_expr() == s.expand(expected_boundary)
facts['analytic_recurrence_pivot'] = str(s.factor(ip.subs(rho,s.symbols('n',integer=True,nonnegative=True)+3)))
series=[]
for seed in range(3):
    aa=[s.Integer(int(i==seed)) for i in range(3)]
    for n in range(9):
        pivot=s.Integer(0); remainder=s.Integer(0)
        for j,p in enumerate(pc):
            for (d,),c in p.terms():
                index=n+j-d
                if index < j: continue
                fac=c*fall(index,j)
                if index == n+3: pivot+=fac
                elif index < len(aa): remainder+=fac*aa[index]
                else: assert fac == 0,(n,j,d,index)
        assert pivot != 0
        aa.append(s.cancel(-remainder/pivot))
    sol=sum(c*t**i for i,c in enumerate(aa))
    residual=s.Poly(s.expand(sum(p.as_expr()*s.diff(sol,t,j) for j,p in enumerate(pc))),t)
    checks['three_analytic_Frobenius_seed_'+str(seed)] = all(residual.nth(n)==0 for n in range(9))
    series.append([str(c) for c in aa])
facts['analytic_Frobenius_series_z8over9'] = series
facts['analytic_dimension_each_boundary'] = 3
facts['common_analytic_dimension_lower_bound_two_endpoints'] = 3+3-5

# Actual logarithmic density as a relative two-variable rational form.
hh=s.symbols('h',positive=True)
rr,ss=s.symbols('rr ss')
lr=lambda x: 1+b+(1-b)*x*x
for er in [-1,1]:
    for es in [-1,1]:
        rp=er*s.I*hh;sp=es*s.I*hh
        residue=(1-b*b)/(2*(1-b)*rp*2*(1-b)*sp)
        residue=s.cancel(residue.subs(b,(hh*hh-1)/(hh*hh+1)))
        checks[f'double_residue_{er}_{es}'] = zero(residue+s.Rational(er*es,4))
        facts[f'double_residue_{er}_{es}']=str(residue)
checks['density_derivative'] = zero(((1-b*b)/(lr(rr)*lr(ss))).subs(b,(hh*hh-1)/(hh*hh+1)) - (hh/(hh*hh+rr*rr))*(hh/(hh*hh+ss*ss)))

# Coefficients read from the independent phase-three decomposition.
yy,Z=s.symbols('yy Z')
c01=-s.sqrt(2)*Z*(9*Z*b*b-18*Z-8*b*b+8)/576
c21=-s.sqrt(2)*Z*(1-b)*(90*Z*yy*yy-45*Z-18*b-32*yy*yy+16)/288
c03=Z*Z*(1-b*b)/18
c23=-2*Z*Z*(1-b)*(2*yy*yy-1)/9
cut1=s.cancel((c01-c21*(1+b)/(1-b)/q**2)/4)
cut3=s.cancel((c03-c23*(1+b)/(1-b)/q**2)/4)
for deg in [1,2]:
    coeff=s.expand(cut1).coeff(Z,deg)
    facts['panel12_same_sign_cut_Z'+str(deg)]=str(s.factor(coeff))
    val=s.simplify(coeff.subs({b:0,yy:s.Rational(-1,2),q:s.Rational(1,2)}))
    checks['panel12_cut_exchange_fixed_nonzero_Z'+str(deg)]=val != 0
    facts['panel12_cut_exchange_fixed_Z'+str(deg)]=str(val)
facts['panel34_same_sign_cut_Z2']=str(s.factor(cut3/Z**2))
checks['panel34_cut_exchange_fixed_nonzero']=zero(cut3.subs(yy**2,s.Rational(1,2))-Z**2*(1-b*b)/72)

# Antisymmetry only cancels the symmetric segment [1,h]. It does not
# cancel the physical [0,min(1,h)] chain. h=1 supplies a simple counterexample.
mp.mp.dps=45
W=lambda x:mp.sqrt(1+(-1)*x*x+x**4)
g=lambda x:2/mp.pi*((mp.pi/2-mp.atan(x))**2-mp.atan(x)**2)/W(x)
moment=mp.quad(g,[0,1])
checks['involution_does_not_cancel_physical_N0']=moment > mp.mpf('.1')
facts['N0_h1_kminus1']=str(moment)

# Trigonometric coordinates used by the spectral evaluator: Y uncentered.
XX,YY=s.symbols('XX YY',real=True)
alpha=(XX+YY)/2;beta=(YY-XX)/2
checks['angular_boundary_identity']=s.trigsimp(1-s.cos(alpha)**2-s.cos(beta)**2+s.cos(XX)*s.cos(YY)) == 0
tt=s.symbols('tt')
checks['x1_unit_divisor_ramifies_order2']=s.limit((s.cos(tt)-1)/tt**2,tt,0)==-s.Rational(1,2)
facts['coordinate_warning']='X=alpha-beta, Y=alpha+beta. x1=cos(alpha),x2=cos(beta); x_i=+/-1 pull back quadratically, whereas a generic cos(X)=0 edge is simple.'

assert all(checks.values()),[n for n,ok in checks.items() if not ok]
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
if sys.platform != 'darwin': rss*=1024
result={'checks':checks,'facts':facts,'elapsed_seconds':time.time()-start,'peak_rss_bytes':rss,
        'conjecture_A':'OPEN','full_relative_GrW':'OPEN',
        'scope':'Exact measure/ladder/Frobenius/residue checks; no physical chain projection computed.'}
(out/'x21_verify.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'check_count':len(checks),'all_passed':True,'elapsed_seconds':result['elapsed_seconds'],
                  'peak_rss_bytes':rss,'facts':facts},indent=2))
