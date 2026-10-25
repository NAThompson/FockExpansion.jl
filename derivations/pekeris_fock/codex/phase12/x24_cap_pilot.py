"""Small serial floating-point screen of the actual four-dimensional cap.

This measures a candidate boundary-value variation, NOT a T5 pairing.
"""
import sys
sys.dont_write_bytecode=True
import resource,time,json
resource.setrlimit(resource.RLIMIT_CPU,(600,600))
from pathlib import Path
import sympy as S
import numpy as np
OUT=Path(__file__).resolve().parent;START=time.time()
cut=json.loads((OUT/'x24_cut.json').read_text())
n,z,R1,R2,u,v,q=S.symbols('n z R1 R2 u v q');loc={str(x):x for x in (n,z,R1,R2,u,v,q)}
compiled=[]
for item in cut['terms']:
 row=[S.sympify(x,locals=loc) for x in item['cut_row']]
 fn=S.lambdify((n,z,R1,R2,u,v,q),row,modules='numpy',cse=True)
 compiled.append((item,fn))
def gl(order,lo=0.,hi=1.):
 x,w=np.polynomial.legendre.leggauss(order)
 return lo+(hi-lo)*(x+1)/2,w*(hi-lo)/2

results=[]
for N in (8,12,18,26,36):
 begin=time.time();Q,WQ=gl(N)
 # xi=tanh(s), retaining the logarithm 2s exactly; tail is exponentially small.
 SX,WX=gl(N+12,0,16);XI=np.tanh(SX);DX=1/np.cosh(SX)**2
 TH,WTH=gl(N,0,np.pi/2);T=np.sin(TH)**2
 PH,WPH=gl(N+8,0,np.pi/2);SP=np.sin(PH)**2;JP=2*np.sin(PH)*np.cos(PH)
 T=T[:,None,None];XI=XI[None,:,None]
 baseweight=(2*WTH[:,None,None])*(WX*DX*2*SX)[None,:,None]*WPH[None,None,:]
 byterm=np.zeros(8);bysector={1:0.,2:0.,3:0.};cap_count=0
 maxvalue=0.;all_finite=True;minH=1.;minS=1.
 for qq,wq in zip(Q,WQ):
  D=1+qq**4;cc=(1-qq*qq*XI*XI)/(1+qq*qq*XI*XI)
  den=1+qq*qq*XI*XI
  PP={1:np.array([qq**4-1,2*qq*qq])/D,2:np.array([1-qq**4,2*qq*qq])/D,
      3:np.array([-2*qq*qq,1-qq**4])/D,4:np.array([2*qq*qq,1-qq**4])/D}
  for idx,(item,fn) in enumerate(compiled):
   aa=np.array([float(x) for x in item['a']]);dd=np.array([float(x) for x in item['d']]);pp=PP[item['panel']]
   r01=T*aa[0]+(1-T)*pp[0];r02=T*aa[1]+(1-T)*pp[1]
   zd=cc*(dd[0]*.2+dd[1]*1.2);z0=r01*.2+r02*1.2
   # Intersect the real tau interval with z(tau)>1, without sampling a Heaviside.
   dz=z0-zd
   tc=np.divide(1-zd,dz,out=np.zeros_like(z0+zd),where=dz!=0)
   lower=np.where(dz>0,np.maximum(0,tc),0)
   upper=np.where(dz<0,np.minimum(1,tc),1)
   bothbelow=(z0<=1)&(zd<=1);bothabove=(z0>1)&(zd>1)
   lower=np.where(bothabove,0,lower);upper=np.where(bothabove,1,upper)
   length=np.where(bothbelow,0,np.maximum(0,upper-lower))
   active=length>0
   tau=lower+length*SP[None,None,:]
   rr1=tau*r01+(1-tau)*cc*dd[0];rr2=tau*r02+(1-tau)*cc*dd[1]
   nn=rr1*rr1+rr2*rr2;zz=.2*rr1+1.2*rr2
   valid=np.broadcast_to(active,nn.shape)&(zz>1)&(nn<1)
   # Avoid evaluating rational formulas on the inactive set.
   nval=np.where(valid,nn,.8);zval=np.where(valid,zz,1.1)
   r1val=np.where(valid,rr1,.1);r2val=np.where(valid,rr2,.9)
   A,C=fn(nval,zval,r1val,r2val,.2,1.2,qq)
   val=A+C*np.sqrt((1-nval)/(zval*zval-nval))
   val=np.where(valid,val,0.)
   all_finite=all_finite and bool(np.all(np.isfinite(val)))
   maxvalue=max(maxvalue,float(np.max(np.abs(val))))
   if np.any(valid):
    minH=min(minH,float(np.min(1-nn[valid])));minS=min(minS,float(np.min(zz[valid]**2-nn[valid])))
   weight=baseweight*(length*JP[None,None,:])*XI/(np.sqrt(D)*den)
   bare=float(np.sum(val*weight))*wq
   pref=(-1 if item['psi40_sector'] in (1,2) else 1)*4*np.sqrt(2)/np.pi**2
   byterm[idx]+=pref*bare;cap_count+=int(np.sum(valid))
 for item,val in zip(cut['terms'],byterm):bysector[item['psi40_sector']]+=float(val)
 record={'order':N,'by_term':byterm.tolist(),'by_sector':bysector,'meaning':'Candidate (upper-minus-lower B contribution)/(-2*pi*i), not T5 pairing',
         'all_finite':all_finite,'cap_nodes':cap_count,'max_abs_local_cut':maxvalue,'min_1_minus_n':minH,'min_z2_minus_n':minS,
         'elapsed_seconds':time.time()-begin}
 results.append(record);print(json.dumps(record),flush=True)
 # This is a convergence screen. No hard success criterion on unvalidated values.
 if not all_finite:raise AssertionError('nonfinite cap value')
facts={'status':'NUMERICAL screen, not an enclosure, T5 projection, or proof of global chain transport',
       'target':{'u':'1/5','v':'6/5'},'quadrature':'q and t-angle Gauss; xi=tanh(s), s<16; cap tau exact interval, sin^2 angle',
       'results':results,'checks':{'all_sampled_values_finite':all(r['all_finite'] for r in results)},
       'elapsed_seconds':time.time()-START}
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;facts['peak_rss_bytes']=rss if sys.platform=='darwin' else rss*1024
(OUT/'x24_cap_pilot.json').write_text(json.dumps(facts,indent=2)+'\n')
print(json.dumps({'elapsed_seconds':facts['elapsed_seconds'],'peak_rss_bytes':facts['peak_rss_bytes']},indent=2))
