"""Transport the actual regular zonal germ along specified target paths.

Small serial Taylor stepping; verifies physical branches at rational fibres.
This is not a transport of the full limiting relative chain to T5.
"""
import sys
sys.dont_write_bytecode=True
import resource,time,json
resource.setrlimit(resource.RLIMIT_CPU,(600,600))
from pathlib import Path
import sympy as S
import mpmath as m
m.mp.dps=75
OUT=Path(__file__).resolve().parent;START=time.time()
exact=json.loads((OUT/'x24_exact.json').read_text());checks={};facts={'fibres':[]}
def mpq(x):return m.mpf(str(S.N(x,85)))
def test(name,err,tol=m.mpf('2e-34')):
 checks[name]=bool(err<tol);print(('PASS ' if checks[name] else 'FAIL ')+name+' '+m.nstr(err,7),flush=True)
 if not checks[name]:raise AssertionError(name)
def actualJ(n,z):
 rho=m.sqrt(n);ga=m.acos(z/rho);theta=2*m.atan((1+rho)*m.tan(ga/2)/m.sqrt(1-n))
 return (ga-theta)/m.sqrt(n-z*z)
def taylor(n,zc,jc,lc,step,order):
 # (1-z)(n-z^2)J'=z(1-z)J+z-1+sqrt(1-n).
 D=[(1-zc)*(n-zc*zc),3*zc*zc-2*zc-n,3*zc-1,1]
 A=[zc*(1-zc),1-2*zc,-1]
 C=[zc-1+m.sqrt(1-n),1]
 coeff=[jc]
 for k in range(order):
  val=sum(A[i]*coeff[k-i] for i in range(min(2,k)+1))+(C[k] if k<2 else 0)
  val-=sum(D[i]*(k-i+1)*coeff[k-i+1] for i in range(1,min(3,k)+1))
  coeff.append(val/((k+1)*D[0]))
 jj=coeff[-1]
 for c in reversed(coeff[:-1]):jj=jj*step+c
 ratio=step/(1-zc)
 ll=lc-sum(ratio**k/k for k in range(1,order+1))
 return jj,ll

def transport(n,zbase,zend,R2,sign,order):
 # Target u=1/5 fixed; v path uses height 1/5, then endpoint v=6/5.
 h=R2/m.mpf(5)
 nodes=[zbase,zbase+sign*1j*h,zend+sign*1j*h,zend]
 jc=actualJ(n,zbase);lc=m.log(1-zbase);steps=0
 # Explicit path stepping with each local step inside 1/5 the nearest pole.
 zc=zbase
 for end in nodes[1:]:
  while abs(end-zc)>m.mpf('1e-65'):
   distance=min(abs(zc-1),abs(zc-m.sqrt(n)),abs(zc+m.sqrt(n)))
   delta=end-zc;sz=min(abs(delta),m.mpf('.18')*distance)
   step=delta*sz/abs(delta)
   jc,lc=taylor(n,zc,jc,lc,step,order);zc+=step;steps+=1
 return jc,lc,steps

samples=[];seen=set()
for item in exact['physical_terms']:
 # Same electron panel at the two source charges has identical geometry.
 geom=(item['panel'],tuple(item['a']))
 if geom in seen:continue
 seen.add(geom)
 if item['panel'] in (1,2):qq,xxi,tt,ta=S.Rational(3,4),S.Rational(1,2),S.Rational(4,5),S.Rational(9,10)
 else:qq,xxi,tt,ta=S.Rational(1,4),S.Rational(1,2),S.Rational(1,2),S.Rational(1,5)
 for inactive in ([False,True] if item['panel'] in (1,3) and item['a'][0]!='-1' else [False]):
  taus=S.Rational(1,5) if inactive and item['panel']==1 else S.Rational(9,10) if inactive else ta
  d=S.Matrix([S.Rational(x) for x in item['d']]);a=S.Matrix([S.Rational(x) for x in item['a']]);qS=S.Symbol('q')
  P=S.Matrix([S.sympify(x,locals={'q':qS}).subs(qS,qq) for x in item['P']]);D=1+qq**4
  cc=(1-qq*qq*xxi*xxi)/(1+qq*qq*xxi*xxi)
  R=taus*(tt*a+(1-tt)*P/D)+(1-taus)*cc*d
  nn=S.factor(R.dot(R));zb=S.factor(R[0]/5+R[1]/3);ze=S.factor(R[0]/5+6*R[1]/5)
  n,base,end,r2=map(mpq,[nn,zb,ze,R[1]])
  assert nn<1 and R[1]>0
  cap=bool(ze>1)
  results=[]
  for order in (50,65):
   jp,lp,sp=transport(n,base,end,r2,1,order)
   jm,lm,sm=transport(n,base,end,r2,-1,order)
   results.append((jp-jm,lp-lm,sp+sm))
  expectedj=-2j*m.pi/m.sqrt(end*end-n) if cap else 0
  expectedl=-2j*m.pi if cap else 0
  idx=len(samples)
  test(f'physical_fibre_{idx}_J_jump',abs(results[-1][0]-expectedj))
  test(f'physical_fibre_{idx}_log_jump',abs(results[-1][1]-expectedl))
  test(f'physical_fibre_{idx}_order_refinement',abs(results[-1][0]-results[0][0]))
  sample={'panel':item['panel'],'a':item['a'],'q':str(qq),'xi':str(xxi),'t':str(tt),'tau':str(taus),
          'R':[str(x) for x in R],'n':str(nn),'zbase':str(zb),'zend':str(ze),'cap':cap,
          'path_v':['1/3','1/3+I/5','6/5+I/5','6/5'],'other_path':'complex conjugate',
          'orders':[50,65],'steps_finest':results[-1][2],
          'J_jump_error':m.nstr(abs(results[-1][0]-expectedj),15),
          'refinement_error':m.nstr(abs(results[-1][0]-results[0][0]),15)}
  samples.append(sample);print(json.dumps(sample),flush=True)
facts['fibres']=samples;facts['checks']=checks
facts['scope']='Numerically verified based physical-germ transport at rational interior fibres. Exact cap incidence agrees. No enclosure for the full endpoint chain, no T5 pairing.'
facts['elapsed_seconds']=time.time()-START
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;facts['peak_rss_bytes']=rss if sys.platform=='darwin' else rss*1024
(OUT/'x24_transport.json').write_text(json.dumps(facts,indent=2)+'\n')
print(json.dumps({'checks':len(checks),'fibres':len(samples),'elapsed_seconds':facts['elapsed_seconds'],'peak_rss_bytes':facts['peak_rss_bytes']},indent=2))
