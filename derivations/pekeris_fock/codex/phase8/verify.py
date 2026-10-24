"""Small independent checks against phase7 cover and actual integrals."""
import json,sys
from pathlib import Path
import sympy as S
import mpmath as m
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
sys.path.insert(0,str(OUT));import rational_panel as R
m.mp.dps=55
new=json.loads((OUT/'rational_connection.json').read_text())
old=json.loads((ROOT/'codex/phase7/connections.json').read_text());on=old['atoms'];nn=new['atoms']
def evaluator(data,variables):
 syms=S.symbols(variables);f=S.lambdify(syms,[S.sympify(e[2],locals=dict(zip(variables.split(),syms))) for e in data['entries']],'mpmath',cse=True)
 def ev(*args):return {(row[0],row[1]):v for row,v in zip(data['entries'],f(*args))}
 return ev
newf={a:evaluator(new['matrices'][a],'H k') for a in ['H','k']}
oldf={a:evaluator(old['cover']['plus_'+a],'h z') for a in ['h','z']}
def transform(H,k):
 h=m.sqrt(H);z=m.sqrt((-k+m.sqrt(k*k-4*H+0j))/2)
 zs=[z,-z,h/z,-h/z,1j,-1j,1j*h,-1j*h]
 scale=m.sqrt(2*H/(1+H));out={}
 for row,(j,fam,n) in enumerate(R.spec):
  jet='00' if n==0 else (str(n)+'0' if fam=='A' else '0'+str(n))
  fac=scale*(h**n if fam=='B' else 1)
  if j<2:out[row,on.index(('J' if j==0 else 'K')+'_'+jet)]=fac
  else:
   weights=[0]*8
   if j<6:
    for p in range(4):weights[p]=zs[p]**(j-2)/(4*zs[p]**3+2*k*zs[p])
   elif j in (6,7):weights[4:6]=[1/(2j),-1/(2j)] if j==6 else [m.mpf('.5')]*2
   else:weights[6:8]=[1/(2j*h),-1/(2j*h)] if j==8 else [m.mpf('.5')]*2
   for p,w in enumerate(weights):
    if w:out[row,on.index(f'J{p}_{jet}')]=fac*w
 out[nn.index('FP_B1'),on.index('FP_01')]=scale*h
 for fam,n in R.jets:
  jet='00' if n==0 else (str(n)+'0' if fam=='A' else '0'+str(n))
  out[R.endpoint(fam,n),on.index('E_'+jet)]=scale*(h**n if fam=='B' else 1)
 out[nn.index('F0'),on.index('F0_00')]=scale
 return out,h,z

def mul(A,B):
 out={};by={}
 for (i,j),v in B.items():by.setdefault(i,[]).append((j,v))
 for (i,j),v in A.items():
  for l,w in by.get(j,[]):out[i,l]=out.get((i,l),0)+v*w
 return out
worst=m.mpf(0);count=0
for H,k in [(m.mpf('1.7'),m.mpf('-.8')),(m.mpf('.6'),m.mpf('.4')),(m.mpf('2.3'),m.mpf('1.1'))]:
 B,h,z=transform(H,k)
 for axis in ['H','k']:
  Mh=oldf['h'](h,z);Mz=oldf['z'](h,z)
  hvel=1/(2*h) if axis=='H' else 0;zvel=-(1 if axis=='H' else z*z)/(4*z**3+2*k*z)
  Mo={ij:hvel*Mh.get(ij,0)+zvel*Mz.get(ij,0) for ij in set(Mh)|set(Mz)}
  lhs=mul(newf[axis](H,k),B);rhs=mul(B,Mo)
  for ij in set(lhs)|set(rhs)|set(B):
   db=m.diff(lambda t:transform(t,k)[0].get(ij,0),H) if axis=='H' else m.diff(lambda t:transform(H,t)[0].get(ij,0),k)
   err=abs(lhs.get(ij,0)-rhs.get(ij,0)-db);worst=max(worst,err);count+=1
assert worst<m.mpf('1e-47'),worst
print('PASS independent phase7 cover gauge identity:',count,'entries; worst',m.nstr(worst,5),flush=True)
# Check all rational Hermite certificates algebraically for the actual requests.
cert=0
# Rebuilding populates certificate cache. The cache keys are not public, so audit
# every top-level integrand explicitly, plus actual vector derivatives below.
for axis in [R.H,R.k]:
 l=(1/(2*R.H) if axis==R.H else 0)-S.diff(R.T,axis)/(2*R.T)
 for ker in R.kernels:
  r=R.C(S.diff(ker,axis)+ker*l);q,co=R.base_reduce(r)
  assert R.C(S.diff(q,R.s)-S.diff(R.T,R.s)*q/(2*R.T)+sum(v*R.kernels[j] for j,v in co.items())-r)==0
  cert+=1
print('PASS',cert,'exact rational Hermite certificates',flush=True)

def actual(H,k):
 h=m.sqrt(H);f=lambda s:h/m.sqrt(s**4+k*s*s+H)
 angle=lambda fam,s:m.atan(s) if fam=='A' else h*m.atan(s/h)
 out=[]
 for j,fam,n in R.spec:
  ker=lambda s:[1,s,s**0/(s**4+k*s*s+H),s/(s**4+k*s*s+H),s*s/(s**4+k*s*s+H),s**3/(s**4+k*s*s+H),1/(1+s*s),s/(1+s*s),1/(H+s*s),s/(H+s*s)][j]
  out.append(m.quad(lambda s:ker(s)*f(s)*angle(fam,s)**n/m.factorial(n),[0,1]))
 out.append(m.quad(lambda s:f(s)*angle('B',s)/s,[0,1]))
 out.extend(f(1)*angle(fam,1)**n/m.factorial(n) for fam,n in R.jets);out.append(m.mpf(1))
 return out
H=m.mpf('.83');k=m.mpf('-.41');X=actual(H,k);worst=0
# Symmetric finite differences give an independent actual-integral derivative.
step=m.mpf('1e-18')
for ax in ['H','k']:
 xp=actual(H+step,k) if ax=='H' else actual(H,k+step)
 xm=actual(H-step,k) if ax=='H' else actual(H,k-step)
 out=[m.mpf(0)]*55
 for (i,j),v in newf[ax](H,k).items():out[i]+=v*X[j]
 worst=max(worst,max(abs(out[i]-(xp[i]-xm[i])/(2*step)) for i in range(55)))
assert worst<m.mpf('1e-32'),worst
print('PASS 110 actual-integral derivatives; worst',m.nstr(worst,5),flush=True)

def moments(H,k,X):
 at=lambda name:X[nn.index(name)]
 def I2(fam,n):
  suf=('A' if n==0 else fam)+str(n)
  v=at('E_'+suf)-k*at('J_'+suf)-(2*H-k*k)*at('Q2_'+suf)+k*H*at('Q0_'+suf)
  if n:
   prev=('A' if n==1 else fam)+str(n-1)
   v-=at('K_'+prev)-at('U1_'+prev) if fam=='A' else H*at('K_'+prev)-H*H*at('V1_'+prev)
  return v
 F0=m.pi/2*at('J_A0')-2*at('J_A1')+4/m.pi*(at('J_A2')-at('J_B2')/H)
 Im=-2*at('Q2_B2')-k*at('Q0_B2')+at('FP_B1')-at('V1_B1')-at('E_B2')
 F2=4/m.pi*Im-m.pi/2*I2('A',0)+2*I2('A',1)-4/m.pi*I2('A',2)
 return F0,F2
F=moments(H,k,X);h=m.sqrt(H)
for j in (0,2):
 def fun(s):
  A=m.pi/2-m.atan(s);D=m.atan(s/h);f=h/m.sqrt(s**4+k*s*s+H)
  return 2/m.pi*((A*A-D*D) if j==0 else ((h*D/s)**2-(s*A)**2))*f
 ref=m.quad(fun,[0,1]);assert abs(ref-F[j//2])<m.mpf('1e-48')
print('PASS N0/N2 output functionals against direct integrals',flush=True)
# Check the old min-endpoint physical moments in both exchange orderings.
sys.path.insert(0,str(ROOT/'codex/phase3'))
import panel_connection as P
for xx,yy,cc in [(m.mpf(2),m.mpf(3),m.mpf(4)),(m.mpf(3),m.mpf(2),m.mpf(4))]:
 Q=xx*xx+yy*yy;ee=2*Q-cc*cc;hh=cc*cc/ee;kk=2*(xx*xx-yy*yy)/ee
 ff=moments(hh,kk,actual(hh,kk));b=(cc*cc-Q)/Q;Y=-xx/m.sqrt(Q)
 physical=P.direct_moments(b,Y)
 assert max(abs(m.sqrt(Q)/cc*ff[j]-physical[j]) for j in (0,1))<m.mpf('1e-46')
 co0=4*(1-b*b)/18;co2=-8*(1-b)*(2*Y*Y-1)/9
 direct=Q**m.mpf('1.5')*(co0*physical[0]+co2*physical[1])
 got=4*cc*ee/18*ff[0]-8*ee*(xx*xx-yy*yy)/(9*cc)*ff[1]
 assert abs(direct-got)<m.mpf('1e-45')
print('PASS physical panel3/panel4 moments and homogeneous psi30 output',flush=True)
# Quantify the boundary counterexample, x=1,y=2,c=1+epsilon.
for eps in [m.mpf('1e-4'),m.mpf('1e-6')]:
 x=1;y=2;c=1+eps;e=2*(x*x+y*y)-c*c;H=c*c/e;k=2*(x*x-y*y)/e;s0=m.sqrt(-k/2)
 val=m.quad(lambda s:m.sqrt(H)/(s**4+k*s*s+H)**m.mpf('1.5'),[0,s0,1])
 print('COUNTEREXAMPLE epsilon',m.nstr(eps,3),'epsilon*Q0_A0',m.nstr(eps*val,15),flush=True)
print('PASS all small serial checks',flush=True)
