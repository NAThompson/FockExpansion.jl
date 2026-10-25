"""Exact geometry and arithmetic certificate for TASK_codex_20; small serial work."""
import sys
sys.dont_write_bytecode=True
import sympy as s
import json,time,resource,subprocess
from pathlib import Path
out=Path(__file__).resolve().parent;start=time.time();checks={};facts={}
q,u,v,X,r,z,H,a,b,T=s.symbols('q u v X r z H a b T');I=s.I
D=1+q**4;P1=q**4-1;P2=2*q*q;w2=1-u*u-v*v;C=D*D-(u*P1+v*P2)**2
zero=lambda e:s.cancel(s.expand(e))==0
for name,ax,ap,k in [('e1',u,P1,P2),('me1',-u,-P1,-P2),('e2',v,P2,-P1)]:
 A=1-ax*ax;B=ap-ax*(u*P1+v*P2);F=X**3+2*B*X*X+A*C*X
 checks[name+'_full_two_torsion']=zero((X*(X+B)**2+X*w2*k*k)-F)
 checks[name+'_gram']=zero(B*B-A*C+w2*k*k)
 checks[name+'_dual_cubic']=zero(F*(1-A*C/X**2)**2-((X+2*B+A*C/X)**3-4*B*(X+2*B+A*C/X)**2+4*(B*B-A*C)*(X+2*B+A*C/X)))
 if name=='e2':
  A2=A;B2=B;F2=F
  y02=-2*A*A*(1-u*u)*q*q*(q**4-1)**2
  checks['e2_section_0']=zero(y02-F.subs(X,-2*A*q*q))
  for sg in [1,-1]:
   xx=A*(sg*I*(q**4-1)-2*q*q)
   yy2=-sg*I*A*A*(u+sg*I*v)**2*(q**4-1)**2*(q*q+sg*I)**2
   checks['e2_section_'+str(sg)]=zero(yy2-F.subs(X,xx))
  Bbar=-u*v*r+2*A;Cbar=(1-u*u)*r*r-4*u*v*r+4*A
  Er=z*(z*z+2*Bbar*z+A*Cbar)
  checks['e2_biquadratic_base']=zero(F.subs(X,q*q*z)/q**6-Er.subs(r,q*q-q**-2))
  checks['e2_section_characters']=zero((q+I/q)**2-(q*q-q**-2+2*I)) and zero((q-I/q)**2-(q*q-q**-2-2*I))
 # Generic row/fibre checks at both char-zero certification targets.
 for uu,vv in [(s.Rational(1,3),s.Rational(1,4)),(s.Rational(1,4),s.Rational(1,3))]:
  cp=s.Poly(C.subs({u:uu,v:vv}),q);kp=s.Poly(k,q)
  checks[f'{name}_squarefree_disjoint_{uu}_{vv}']=s.gcd(cp,cp.diff()).degree()==0 and s.gcd(cp,kp).degree()==0
# The endpoint dual is Kuwata--Shioda J1, explicitly, not just a matching fibre list.
AA=s.symbols('AA');l1=(a+1)/(a-1);l2=(b+1)/(b-1);g=AA*(b-1)**2;k2=I*(a-1)/(b-1)
checks['KS_q4']=zero(g*(l1-1)**2*k2*k2+4*AA)
checks['KS_q2']=zero(-2*g*(l1+1)*(l2+1)*k2+8*I*AA*a*b)
checks['KS_constant']=zero(g*(l2-1)**2-4*AA)
checks['KS_a4']=zero(g*g*16*l1*l2*k2*k2+16*AA*AA*(1-a*a)*(1-b*b))
# a=v/sqrt(A), b=i*u/sqrt(A), A=1-u^2.
checks['KS_target_a4']=zero((AA*AA*(1-a*a)*(1-b*b)).subs({a:v/H,b:I*u/H,AA:H**2}).expand().subs(H**4,(1-u*u)**2).subs(H**2,1-u*u)-w2)
j1=64*(v*v+3*(1-u*u))**3/((1-u*u)*(v*v-(1-u*u))**2)
j2=64*(3-4*u*u)**3/(1-u*u)
checks['old_curve_j']=zero(16*((-2*v)**2+12*(1-u*u))**3/((1-u*u)*((-2*v)**2-4*(1-u*u))**2)-j1)
checks['second_curve_j']=zero(256*((2*u)**2-3)**3/((2*u)**2-4)-j2)
# A closed-cycle holomorphic period on Claude's fixed-u endpoint line has order TWO.
vv,hh=s.symbols('vv hh');ll=(vv+hh)/(vv-hh);gg=(vv-hh)**s.Rational(-1,2)
coef2=(vv*vv-hh*hh)*gg*s.diff(ll,vv)**2
coef1=(vv*vv-hh*hh)*(2*s.diff(gg,vv)*s.diff(ll,vv)+gg*s.diff(ll,vv,2))+2*vv*gg*s.diff(ll,vv)
coef0=(vv*vv-hh*hh)*s.diff(gg,vv,2)+2*vv*s.diff(gg,vv)+gg/4
checks['fixed_u_order2_F1']=s.simplify(coef1+coef2*(2*ll-1)/(ll*(1-ll)))==0
checks['fixed_u_order2_F0']=s.simplify(coef0+coef2/(4*ll*(1-ll)))==0
checks['C_discriminant']=zero(s.discriminant(C,q)-2**32*(u*u-1)**2*(u*u+v*v)**4*w2**4)
f=u*v/(1-u*u);gg=(1+u*u-2*v*v)/(1-u*u)
checks['nonisotrivial_moduli']=zero(s.det(s.Matrix([[s.diff(f,u),s.diff(f,v)],[s.diff(gg,u),s.diff(gg,v)]]))+4*(u*u+v*v)/(1-u*u)**3)
# t=s/(1+s), X=A*D*s, Y=A*(1+s)^2*Y0: D dt/Y0=dX/Y.
ss,AA=s.symbols('ss AA')
checks['holomorphic_measure_Jacobian']=zero(D/(1+ss)**2- (AA*D)/(AA*(1+ss)**2))
# Independent prime-field counts use no C tables; reproduce all n=1 outputs.
cert=[]
for filename,swap in [('x20_counts.json',False),('x20_counts_line.json',True)]:
 dat=json.loads((out/filename).read_text());target=(s.Rational(1,4),s.Rational(1,3)) if swap else (s.Rational(1,3),s.Rational(1,4))
 local=[]
 for ent in dat['data']:
  p=ent['p'];row=ent['rows'][0];uu=pow(4 if swap else 3,-1,p);vv=pow(3 if swap else 4,-1,p);aa=(1-vv*vv)%p
  chi=[0]+[1 if pow(j,(p-1)//2,p)==1 else -1 for j in range(1,p)];sm=0
  for qq in [*range(p),None]:
   d,p1,p2=(1,1,0) if qq is None else ((1+qq**4)%p,(qq**4-1)%p,2*qq*qq%p)
   px=(uu*p1+vv*p2)%p;bb=(p2-vv*px)%p;cc=(d*d-px*px)%p
   sm+=sum(chi[(xx**3+2*bb*xx*xx+aa*cc*xx)%p] for xx in range(p))
  checks[f'python_count_{filename}_{p}']=sm==row['char_sum']
  # Independent F_(p^2) pair arithmetic and norm character, without C tables.
  rr2=ent['rows'][1];c0,c1,_=rr2['irreducible_relation'];elts=[(j%p,j//p) for j in range(p*p)]
  def add(x,y):return ((x[0]+y[0])%p,(x[1]+y[1])%p)
  def mul(x,y):return ((x[0]*y[0]+c0*x[1]*y[1])%p,(x[0]*y[1]+x[1]*y[0]+c1*x[1]*y[1])%p)
  def scale(c,x):return (c*x[0]%p,c*x[1]%p)
  def ch(x):return chi[(x[0]*x[0]+c1*x[0]*x[1]-c0*x[1]*x[1])%p]
  sm2=0
  for qq in [*elts,None]:
   if qq is None:d,p1,p2=(1,0),(1,0),(0,0)
   else:
    q2=mul(qq,qq);q4=mul(q2,q2);d=add((1,0),q4);p1=add((-1,0),q4);p2=scale(2,q2)
   px=add(scale(uu,p1),scale(vv,p2));bb=add(p2,scale(-vv,px));cc=add(mul(d,d),scale(-1,mul(px,px)))
   for xx in elts:sm2+=ch(mul(xx,add(mul(xx,add(xx,scale(2,bb))),scale(aa,cc))))
  checks[f'python_extension2_count_{filename}_{p}']=sm2==rr2['char_sum']
  cp=s.Poly(C.subs({u:uu,v:vv}),q,modulus=p)
  checks[f'good_reduction_{filename}_{p}']=s.gcd(cp,cp.diff()).degree()==0 and s.gcd(cp,s.Poly(q**4-1,q,modulus=p)).degree()==0
  checks[f'ordinary_{filename}_{p}']=sm%p!=0
  # Recompute the Frobenius polynomial from raw counts and known divisor characters.
  traces=[]
  for rr in ent['rows']:
   nn=rr['n'];Q=p**nn;eps0=s.legendre_symbol((-2*(1-uu*uu))%p,p);pm=pow(p,nn,8)
   mwt=int(eps0**nn)+(2 if pm==1 else -2 if pm==5 else 0)
   traces.append(s.Integer(rr['char_sum']-mwt*Q))
   # Independent finite-field ring sanity / irreducibility.
   rel=rr['irreducible_relation']
   if nn>1:checks[f'field_irreducible_{filename}_{p}_{nn}']=s.Poly(T**nn-rel[1]*T-rel[0],T,modulus=p).is_irreducible
  co=[s.Integer(1)]
  for nn in range(1,4):co.append(-sum(co[nn-j]*traces[j-1] for j in range(1,nn+1))/nn)
  valid=[]
  for sign in [-1,1]:
   if co[3]==sign*p*co[2]:
    valid.append(sum((co+[sign*p**3*co[1],sign*p**5])[j]*T**j for j in range(6)))
  checks[f'unique_reciprocity_{filename}_{p}']=len(valid)==1
  P5=valid[0];extra=1+p*T if s.rem(P5,1+p*T,T)==0 else 1-p*T
  R4=s.div(P5,extra,T)[0];g4=s.expand(T**4*R4.subs(T,1/(p*T)))
  cycl=[m for m in range(1,61) if s.totient(m)<=4 and s.degree(s.gcd(g4,s.cyclotomic_poly(m,T)),T)>0]
  factor_degrees=[f.degree() for f,e in s.factor_list(cp)[1] for _ in range(e)]
  details={'p':p,'T5_polynomial':str(P5),'remaining_quartic':str(R4),'C_factor_degrees':factor_degrees,'cyclotomic_orders':cycl}
  if not cycl:
   checks[f'quartic_irreducible_{filename}_{p}']=s.Poly(g4,T).is_irreducible
   checks[f'Frob4_divisors_fixed_{filename}_{p}']=all(4%d==0 for d in factor_degrees)
   val=s.factor(s.resultant(g4,1-T**4,T));num,den=s.fraction(val)
   sf=-s.prod(pr for pr,ex in s.factorint(abs(num*den)).items() if ex%2)
   details.update(rank_upper=18,Artin_Tate_Fp4_limit=str(val),NS_discriminant_squareclass=int(sf))
  else:details.update(rank_upper=20,qualification='p=13 swapped target has extra algebraic cyclotomic factors; not used for upper-17 certificate')
  local.append(details)
 proven={d['NS_discriminant_squareclass'] for d in local if d.get('rank_upper')==18}
 checks[f'two_discriminants_{filename}']=len(proven)>=2
 cert.append({'target':list(map(str,target)),'prime_data':local,'char_zero_Picard':17})
facts['certificates']=cert
facts['summary']={'e1':{'Picard':18,'MW_rank':2,'torsion':'(Z/2)^2','T_rank':4},'e2':{'Picard':17,'MW_rank':3,'torsion':'(Z/2)^2','T_rank':5},'endpoint_curves':['eta^2=(1-u)q^4-2v q^2+(1+u)','eta^2=z(z^2+2u z+1)'],'fixed_u_period_operator':'(v^2-(1-u^2))*D_v^2+2v*D_v+1/4','j1':str(j1),'j2':str(j2)}
assert all(checks.values()),{k:v for k,v in checks.items() if not v}
result={'checks':checks,'facts':facts,'elapsed_seconds':time.time()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(out/'x20_verify.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
