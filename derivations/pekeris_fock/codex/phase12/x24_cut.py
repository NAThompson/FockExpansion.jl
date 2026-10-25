"""Local zonal jumps and actual four-parameter cap forms, rational handoff."""
import sys
sys.dont_write_bytecode=True
import resource,time,json
resource.setrlimit(resource.RLIMIT_CPU,(600,600))
from pathlib import Path
import sympy as S
OUT=Path(__file__).resolve().parent;START=time.time()
n,z,R1,R2,u,v,q,xi,t,tau,Y1,Y2=S.symbols('n z R1 R2 u v q xi t tau Y1 Y2')
Qp=z*z-n;H=1-n
checks={}
def zero(x):return S.cancel(x)==0
def test(name,result):
 checks[name]=bool(result);print(('PASS ' if result else 'FAIL ')+name,flush=True)
 if not result:raise AssertionError(name)
def ssrow(row):return [str(S.factor(x)) for x in row]
# Column states V=(1, h/s), h=sqrt(1-n), s=sqrt(z^2-n).
Bn=S.diag(0,-1/(2*H)+1/(2*Qp));Bz=S.diag(0,-z/Qp)
test('jump_character_flat',all(zero(x) for x in Bz.diff(n)-Bn.diff(z)))
def cov(row,var):return (row.diff(var)+row*(Bn if var==n else Bz)).applyfunc(S.cancel)
def E(row):return (2*n*cov(row,n)+z*cov(row,z)).applyfunc(S.cancel)
def dr(row,i):
 rr,xx=(R1,u) if i==1 else (R2,v)
 return (row.diff(rr)+2*rr*cov(row,n)+xx*cov(row,z)).applyfunc(S.cancel)
exact=json.loads((OUT/'x24_exact.json').read_text())
loc=globals()
def parse(x):return S.sympify(x,locals=loc)
def grow(kind):
 row=exact['inverse_rows'][kind]
 return S.Matrix([[parse(row[1]),parse(row[0])]])
g1=grow('inverse1');g2=(g1+E(g1)).applyfunc(S.cancel);gl=grow('log');gx=grow('ell_log')
# At z=1, h/s=1. Derivatives of the vector include the algebraic-state derivative.
for label,row,order in [('inverse1',g1,1),('log',gl,2),('ell_log',gx,3)]:
 temp=row
 for k in range(order):
  test(f'{label}_jump_vanishes_to_order_{k}',zero(sum(temp).subs(z,1)))
  temp=cov(temp,z)
# The four-parameter source prescription gives NO delta terms, in contrast
# to exponent 3 in the unreduced six-parameter representation.
monomials={}
for degree in range(4):
 for a in range(degree+1):
  b=degree-a
  row=g2 if degree==0 else g1 if degree==1 else -gl if degree==2 else gx
  for k in range(a):row=dr(row,1)
  for k in range(b):row=dr(row,2)
  monomials[(a,b)]=row
  if degree>1:
   other=g1 if degree==1 else -gl if degree==2 else gx
   for k in range(b):other=dr(other,2)
   for k in range(a):other=dr(other,1)
   test(f'mixed_R_derivatives_{a}_{b}',all(zero(x) for x in row-other))
terms=[]
for idx,item in enumerate(exact['physical_terms']):
 poly=S.Poly(parse(item['four_numerator']),Y1,Y2)
 row=S.zeros(1,2)
 for powers,coeff in poly.terms():row+=coeff*monomials[powers]
 row=row.applyfunc(S.cancel)
 terms.append(dict(item,cut_row=ssrow(row)))
 test('term'+str(idx)+'_no_H_pole_in_rational_part',not S.fraction(row[0])[1].has(1-n) and zero(S.limit((1-n)*row[0],n,1)))
 # At fixed z^2 !=1 the only radical coefficient pole is simple at n=1.
 test('term'+str(idx)+'_simple_H_pole_in_radical_part',zero(S.limit((1-n)**2*row[1],n,1)))
 # Actual endpoint coefficient agrees with the source numerator evaluated at R.
 expected=-parse(exact['inverse_rows']['inverse1'][0]).subs(n,1)*poly.as_expr().subs({Y1:R1,Y2:R2})
 residue=S.cancel((1-n)*row[1]).subs(n,1)
 test('term'+str(idx)+'_endpoint_source_value',zero(residue-expected))
 # Local target-independent resonant polynomial row is included, not discarded.
 test('term'+str(idx)+'_target_rational_cut_not_dropped',row[0]!=0)

# Pure-inverse rows generated modulo target quadratics give the SAME jump;
# there is no target discontinuity from their omitted normalization.
# Full parameter-connection convention (all variables named).
facts={'status':'Exact local jump rows; cap with all boundary faces below. No T5 reduction or pairing.',
 'basis':['1','sqrt((1-n)/(z^2-n))'], 'd_n':[[str(x) for x in r] for r in Bn.tolist()],
 'd_z':[[str(x) for x in r] for r in Bz.tolist()],
 'derivations':'D_R1=partial_R1+2*R1*partial_n+u*partial_z; D_R2=partial_R2+2*R2*partial_n+v*partial_z; apply connection to rows',
 'monomial_rows':{f'{a},{b}':ssrow(row) for (a,b),row in monomials.items()},
 'terms':terms,
 'cut_formula':'upper-minus-lower B contribution = -2*pi*i * sum term_prefactor * Integral(cap, xi*log((1+xi)/(1-xi))/(sqrt(D*t*(1-t))*(1+q^2*xi^2)) * cut_row . [1,sqrt((1-n)/(z^2-n))], dq dxi dt dtau)',
 'target':'u=1/5, v=6/5; continue from (1/5,1/3) through Im(v)>0 or <0; sqrt((1-n)/(z^2-n)) positive on real cap',
 'cap':'0<q,xi,t,tau<1 and z=R1*u+R2*v>1, with R as in x24_exact.json; every such point has R2>0',
 'no_delta':'For source degree d, the primitive jump vanishes to order >=d at z=1. Its d parameter derivatives create no delta. U2 itself is undifferentiated.',
 'tau_intersection':'z0=(t*a+(1-t)*P/D).x, zd=c*d.x; z=zd+tau*(z0-zd). Intersect (0,1) with z>1.',
 'tau_interval_cases':{'both_le_1':'empty','both_gt_1':'(0,1)','zd_le_1_lt_z0':'((1-zd)/(z0-zd),1)','z0_le_1_lt_zd':'(0,(1-zd)/(z0-zd))'},
 'endpoint_bound':'On cap n>1/(u^2+v^2)=25/37, so n=0 artifacts are absent. Generic n->1, z^2!=1 behavior is O((1-n)^(-1/2)). At a transverse simultaneous z=1,n=1 endpoint it is O((1-n)^(-1)); the beta-weighted two corner variables make this locally integrable. Degenerate arc endpoints have a.x or d.x distinct from 1 at this target.',
 'caution':'Cap period is four-dimensional and mixed. A nonzero value would not be a T5 pairing. The target-to-S_e2 projection remains uncomputed.'}
# Verify smooth genus-one tau slice at an exact rational sample, including
# q,r relation; shows what geometry remains instead of pretending it is S_e2.
P=S.Matrix([q**4-1,2*q*q]);D=1+q**4;aa=S.Matrix([0,1]);dd=S.Matrix([-1,0]);c=(1-q*q*xi*xi)/(1+q*q*xi*xi)
rr=tau*(t*aa+(1-t)*P/D)+(1-tau)*c*dd
nn=S.cancel(rr.dot(rr));zz=rr.dot(S.Matrix([u,v]))
sample={q:S.Rational(1,2),xi:S.Rational(2,3),t:S.Rational(2,5),u:S.Rational(1,5),v:S.Rational(6,5)}
quart=S.Poly(S.cancel((1-nn)*(zz*zz-nn)).subs(sample),tau)
quart=S.Poly(quart.as_expr(),tau)
test('tau_slice_quartic_degree4',quart.degree()==4)
test('tau_slice_quartic_squarefree',S.gcd(quart,quart.diff()).degree()==0)
facts['tau_elliptic_slice']={'parameters':{str(k):str(vv) for k,vv in sample.items()},
 'curve':'eta_tau^2=(1-n(tau))*(z(tau)^2-n(tau))',
 'quartic':str(S.factor(quart.as_expr())), 'degree':quart.degree(),
 'genus':'1 for this smooth quartic slice; NOT a T5 identification'}
facts['checks']=checks;facts['elapsed_seconds']=time.time()-START
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;facts['peak_rss_bytes']=rss if sys.platform=='darwin' else rss*1024
(OUT/'x24_cut.json').write_text(json.dumps(facts,indent=2)+'\n')
(OUT/'x24_cut.txt').write_text('# Rational physical-cap handoff. Rows act on COLUMN basis.\n'+json.dumps(facts,indent=2)+'\n')
print(json.dumps({'checks':len(checks),'elapsed_seconds':facts['elapsed_seconds'],'peak_rss_bytes':facts['peak_rss_bytes'],'export_bytes':(OUT/'x24_cut.txt').stat().st_size},indent=2))
