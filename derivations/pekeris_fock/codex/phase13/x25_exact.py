"""Independent exact local audit of the supplied polynomial L5."""
import json, time, resource, sys, io, pickle, zlib
from fractions import Fraction
from pathlib import Path
import sympy as S
start=time.time(); out=Path(__file__).resolve().parent; root=out.parents[1]
u,v,q,X,z,t,rho=S.symbols('u v q X z t rho')
checks={}; facts={}
P=q**4-1; D=q**4+1; A=1-v*v
B=-u*v*P+2*A*q*q; C=D*D-(u*P+2*v*q*q)**2
checks['discriminant_identity']=S.expand(B*B-A*C+P*P*(1-u*u-v*v))==0
F=-B
checks['circle_factorization']=S.rem(S.expand(F-u*v*(q*q-(u+1)/v)*(q*q-(u-1)/v)),S.Poly(u*u+v*v-1,u),u)==0
Y=S.symbols('Y'); a=S.symbols('a')
checks['central_normalization']=S.expand((a*(a*a+B))**2-a*a*(a*a+B)**2)==0
# Quotient (q,y)->(r=q^2,Z=q*y); pullback dr/Z=2 dq/y.
r=S.symbols('r')
Fr=u*v*r*r-2*(1-v*v)*r-u*v
checks['isogeny_equation']=S.expand(q*q*F-(r*Fr).subs(r,q*q))==0
lam=(u-1)/(u+1)
checks['lambda_ratio']=S.cancel(((u-1)/v)/((u+1)/v)-lam)==0
facts['residue_curve']='y^2=uv*q^4-2*(1-v^2)*q^2-uv; u^2+v^2=1'
facts['quotient_curve']='Z^2=r*(uv*r^2-2*(1-v^2)*r-uv), r=q^2,Z=q*y; dr/Z=2*dq/y'
# Binary-quartic j and quotient-cubic j need not agree.
I=12*(u*v)*(-u*v)+(-2*(1-v*v))**2
J=72*(u*v)*(-2*(1-v*v))*(-u*v)-2*(-2*(1-v*v))**3
jq=S.factor(6912*I**3/(4*I**3-J**2))
jr=256*(1-lam+lam*lam)**3/(lam*lam*(1-lam)**2)
facts['j_quartic_on_circle']=str(S.factor(jq.subs(v*v,1-u*u)))
facts['j_quotient']=str(S.factor(jr))
lp=S.diff(lam,u)
ga=(1-2*lam)*lp/(lam*(1-lam))-S.diff(lam,u,2)/lp
gb=-lp*lp/(4*lam*(1-lam))
g=(u*(u+1))**S.Rational(-1,2); gh=S.diff(g,u)/g
pa=(3*u*u-1)/(u*(u*u-1)); pb=(4*u*u+1)/(4*u*u*(u*u-1))
checks['tangential_PF_gauge_a']=S.simplify(ga-2*gh-pa)==0
checks['tangential_PF_gauge_b']=S.simplify(2*gh*gh-S.diff(g,u,2)/g-ga*gh+gb-pb)==0
facts['tangential_PF']='4*u^2*(u^2-1)*D_u^2+4*u*(3*u^2-1)*D_u+(4*u^2+1)'
checks['isogenous_not_generically_isomorphic']=S.simplify((jq-jr).subs({u:S.Rational(3,5),v:S.Rational(4,5)}))!=0
# Panel quartics at unit radial scale, before harmless algebraic scalings.
panel=[(1-u)*q**4-2*v*q*q+1+u,
       (1+u)*q**4-2*v*q*q+1-u,
       (1+v)*q**4+2*u*q*q+1-v,
       (1+v)*q**4-2*u*q*q+1-v]
for k,f in enumerate(panel,1):
    pp=S.Poly(f,q)
    disc=pp.nth(2)**2-4*pp.nth(4)*pp.nth(0)
    checks['panel_discriminant_'+str(k)]=S.expand(disc+4*(1-u*u-v*v))==0
    lead=pp.nth(4); middle=pp.nth(2)
    difference=S.expand(lead*f-(lead*q*q+middle/2)**2)
    checks['panel_square_on_circle_'+str(k)]=S.expand(difference-(1-u*u-v*v))==0
facts['panel_quartics']=[str(p) for p in panel]
# Read plain-text exported operator; coefficients are low to high.
cs=[sum(S.Integer(c)*u**i for i,c in enumerate(line.split()))
    for line in (root/'claude10/op_e2_A.txt').read_text().splitlines()]
assert len(cs)==6
# Sage's rational pickle stores base-32 integers. Decode only four known
# inert constructors; never import Sage or execute arbitrary pickle globals.
def rational32(s):
    parts=s.split('/')
    return Fraction(int(parts[0],32),int(parts[1],32) if len(parts)>1 else 1)
class ReadOnlySageDecoder(pickle.Unpickler):
    def find_class(self,module,name):
        known={
          ('sage.rings.rational','make_rational'):rational32,
          ('sage.rings.rational_field','RationalField'):lambda: 'QQ',
          ('sage.rings.polynomial.polynomial_ring_constructor','unpickle_PolynomialRing'):lambda *args:'QQ[l]',
          ('sage.rings.polynomial.polynomial_rational_flint','Polynomial_rational_flint'):lambda parent,coeff,*args:coeff}
        if (module,name) not in known: raise ValueError((module,name))
        return known[module,name]
decoded=ReadOnlySageDecoder(io.BytesIO(zlib.decompress((root/'claude10/op_e2_Q.sobj').read_bytes()))).load()
checks['sobj_matches_plain_text_exactly']=all(S.expand(sum(S.Rational(c.numerator,c.denominator)*u**i for i,c in enumerate(row))-p)==0 for row,p in zip(decoded,cs)) and len(decoded)==len(cs)
ders=[{0:S.Integer(1)}]
for j in range(5):
    dd={}
    for n,c in ders[-1].items():
        dd[n]=dd.get(n,0)+S.diff(c,u)
        dd[n+1]=dd.get(n+1,0)+2*u*c
    ders.append(dd)
zc=[S.expand(sum(cs[j]*ders[j].get(n,0) for j in range(6))/u).subs(u*u,z) for n in range(6)]
pc=[S.Poly(c.subs(z,t+S.Rational(8,9)),t) for c in zc]
fall=lambda x,j:S.prod(x-k for k in range(j))
ind=S.factor(sum(p.nth(j-3)*fall(rho,j) for j,p in enumerate(pc) if j>=3))
checks['indicial']=S.Poly(ind,rho).monic().as_expr()==S.expand(rho**2*(rho-1)**2*(rho-2))
facts['indicial']=str(ind)
seeds=S.symbols('a0:3'); aa=list(seeds)
for n in range(16):
    pivot=0; rest=0
    for j,p in enumerate(pc):
        for (d,),c in p.terms():
            ix=n+j-d
            if ix<j: continue
            fac=c*fall(ix,j)
            if ix==n+3: pivot+=fac
            elif ix<len(aa): rest+=fac*aa[ix]
            else: assert fac==0
    aa.append(S.factor(-rest/pivot))
facts['analytic_pivot']=str(S.factor(ind.subs(rho,S.symbols('n')+3)))
# L(h log t) = log(t)Lh + commutator; compute Laurent coefficients.
# D^j(t^n log t) nonlog coefficient is d/dn falling(n,j)*t^(n-j).
cc={}
nr=S.symbols('nr')
for j,p in enumerate(pc):
    df=S.diff(fall(nr,j),nr)
    for (d,),c in p.terms():
        for n,an in enumerate(aa):
            degree=d+n-j
            if degree>10: continue
            cc[degree]=cc.get(degree,0)+c*an*df.subs(nr,n)
neg={k:S.factor(val) for k,val in cc.items() if k<0 and S.expand(val)!=0}
facts['log_obstructions']={str(k):str(val) for k,val in sorted(neg.items())}
mat=S.linear_eq_to_matrix(list(neg.values()),seeds)[0]
checks['two_single_log_solutions']=mat.rank()==1
null=mat.nullspace()
facts['log_coefficient_seed_basis']=[[str(x) for x in vec] for vec in null]
facts['log_constraint']=[str(x) for x in mat.rowspace()[0]]
# All subsequent analytic inhomogeneous coefficients have the same nonzero
# pivot. Thus two independent h log(t)+g solutions + three analytic solutions
# form a basis: Jordan blocks 2,2,1, not 3,1,1.
checks['three_analytic_two_log_basis']=len(null)==2 and ind.subs(rho,3)!=0
facts['analytic_series']=[str(x) for x in aa[:8]]
facts['Jordan_blocks']=[2,2,1]
facts['Chebyshev_log10_rho']=str(S.N(S.log(3+S.sqrt(8),10),35))
facts['Chebyshev_log_model_slope_40_70']=str(S.N(S.log(3+S.sqrt(8),10)+S.log(S.Rational(70,40),10)/30,35))
assert all(checks.values()),[k for k,v in checks.items() if not v]
result={'checks':checks,'facts':facts,'elapsed_seconds':time.time()-start,
        'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)}
(out/'x25_exact.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
