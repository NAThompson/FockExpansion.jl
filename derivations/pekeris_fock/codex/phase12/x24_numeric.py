"""Independent numeric checks of density reduction, Feynman chains and cuts."""
import sys
sys.dont_write_bytecode=True
import resource,time,json
resource.setrlimit(resource.RLIMIT_CPU,(600,600))
from pathlib import Path
import sympy as S
import mpmath as m
import numpy as np
OUT=Path(__file__).resolve().parent;START=time.time();m.mp.dps=55
facts={};checks={}
def test(name,value,tol):
 value=abs(value);checks[name]=bool(value<tol);facts[name]={'error':str(value),'tolerance':str(tol)}
 print(('PASS ' if checks[name] else 'FAIL ')+name+' '+str(value),flush=True)
 if not checks[name]:raise AssertionError(name)
# 1. New single-density identity, at four independently chosen interior points.
for idx,(bb,qq) in enumerate([('0.2','0.7'),('-0.6','0.4'),('0','1'),('0.85','0.9')]):
 b=m.mpf(bb);q=m.mpf(qq);h=m.sqrt((1+b)/(1-b))
 got=(1-b)*q*q*m.quad(lambda x:x*m.log((1+x)/(1-x))/(1+b+(1-b)*q*q*x*x),[0,m.mpf('.5'),1])
 test('density_single_parameter_'+str(idx),got-m.atan(q/h)**2,m.mpf('1e-50'))
# Original double integral separated exactly at fixed source: its 1D factor
# squared is a genuine independent check on both the normalization and h.
b=m.mpf('-.3');q=m.mpf('.65')
i=m.quad(lambda r:1/(1+b+(1-b)*r*r),[0,q])
new=(1-b)*q*q*m.quad(lambda x:x*m.log((1+x)/(1-x))/(1+b+(1-b)*q*q*x*x),[0,1])
test('double_vs_single_density',(1-b*b)*i*i-new,m.mpf('1e-50'))

# Density state derivatives, computed from actual atan values.
bb,qq=m.mpf('.23'),m.mpf('.41')
def density(b,q):
 h=m.sqrt((1+b)/(1-b));B=m.atan(q/h);return m.matrix([B*B,h*B,1])
exact=json.loads((OUT/'x24_exact.json').read_text())
bS,qS=S.symbols('b q')
for var,matrix in [('b',exact['density']['d_b']),('q',exact['density']['d_q'])]:
 fun=S.lambdify((bS,qS),S.Matrix([[S.sympify(x,locals={'b':bS,'q':qS}) for x in row] for row in matrix]),'mpmath')
 actual=m.diff((lambda b:density(b,qq)) if var=='b' else (lambda q:density(bb,q)),bb if var=='b' else qq)
 test('actual_density_connection_'+var,m.norm(actual-fun(bb,qq)*density(bb,qq),m.inf),m.mpf('1e-50'))

# 2. Zonal kernel, using an explicit elementary branch rather than matrices.
def J(n,z):
 s=m.sqrt(z*z-n);h=m.sqrt(1-n)
 return m.log((1-z)*(z+s)/(z-n+h*s))/s

def basis(n,z):return m.matrix([m.sqrt(1-n)*J(n,z),m.log(1-z),1,m.sqrt(1-n)])
nS,zS=S.symbols('n z')
ns,zs=m.mpf('.31'),m.mpf('.14')
for var,matrix in [('n',exact['zonal']['d_n']),('z',exact['zonal']['d_z'])]:
 fun=S.lambdify((nS,zS),S.Matrix([[S.sympify(x,locals={'n':nS,'z':zS}) for x in row] for row in matrix]),'mpmath')
 actual=m.diff((lambda n:basis(n,zs)) if var=='n' else (lambda z:basis(ns,z)),ns if var=='n' else zs)
 test('actual_zonal_connection_'+var,m.norm(actual-fun(ns,zs)*basis(ns,zs),m.inf),m.mpf('1e-49'))
# Same J from the actual angle branch at an independent real point.
gamma=m.acos(zs/m.sqrt(ns));theta=2*m.atan((1+m.sqrt(ns))*m.tan(gamma/2)/m.sqrt(1-ns))
test('J_log_vs_angles',J(ns,zs)-(gamma-theta)/m.sqrt(ns-zs*zs),m.mpf('1e-50'))

# Plain rows are evaluated modulo target degree-two polynomials; these
# polynomials have zero jump and do not affect the following comparisons.
rows={}
for key,rr in exact['inverse_rows'].items():
 rows[key]=S.lambdify((nS,zS),S.Matrix([[S.sympify(x,locals={'n':nS,'z':zS}) for x in rr]]),'mpmath')
def inv(kind,n,z):return (rows[kind](n,z)*basis(n,z))[0]
ns,zs=m.mpf('.37'),m.mpf('1.23');eps=m.mpf('1e-20')
for kind in ('inverse1','inverse2','inverse3','log','ell_log'):
 row=exact['inverse_rows'][kind]
 A=S.lambdify((nS,zS),S.sympify(row[1],locals={'n':nS,'z':zS}),'mpmath')
 C=S.lambdify((nS,zS),S.sympify(row[0],locals={'n':nS,'z':zS}),'mpmath')
 expected=A(ns,zs)+C(ns,zs)*m.sqrt((1-ns)/(zs*zs-ns))
 got=(inv(kind,ns,zs+1j*eps)-inv(kind,ns,zs-1j*eps))/(-2j*m.pi)
 test('actual_local_jump_'+kind,got-expected,m.mpf('2e-18'))

# 3. Feynman identities use independent finite quadrature before inversion.
def gl(order,lo=0.,hi=1.):
 x,w=np.polynomial.legendre.leggauss(order);return lo+(hi-lo)*(x+1)/2,w*(hi-lo)/2
q=.6;rs=.27;ss=.43;Y=np.array([.2,-.3]);D=1+q**4
P=np.array([q**4-1,2*q*q]);p=P/D;aa=np.array([0.,1.]);dd=np.array([-1.,0.])
b=-dd@Y;aa0=1-aa@Y;pp0=1-p@Y
cr,cs=1+rs*rs,1+ss*ss;kr,ks=(1-rs*rs)/cr,(1-ss*ss)/cs
lr=1+b+(1-b)*rs*rs;ls=1+b+(1-b)*ss*ss
th,wth=gl(24,0,np.pi/2);tt=np.sin(th)**2
TA,wa=gl(24);ET,we=gl(24)
T=tt[:,None,None];A=TA[None,:,None];E=ET[None,None,:]
R=A[...,None]*(T[...,None]*aa+(1-T[...,None])*p)+(1-A[...,None])*(E[...,None]*kr+(1-E[...,None])*ks)*dd
ell=1-np.sum(R*Y,axis=-1)
got=2/np.pi*np.sum((2*wth[:,None,None])*wa[None,:,None]*we[None,None,:]*(1-A)/ell**3)/(cr*cs)
expected=1/(np.sqrt(aa0*pp0)*lr*ls)
test('six_parameter_dirichlet_normalization',got-expected,2e-13)
# Four-parameter density-chain identity, fixed q and actual sector-one numerator.
# xi=tanh(s) turns the endpoint logarithm into 2s; truncate at s=20.
SX,wx=gl(56,0,20);XI=np.tanh(SX);dx=1/np.cosh(SX)**2
T=tt[:,None,None];A=TA[None,:,None];XI=XI[None,None,:]
cc=(1-q*q*XI*XI)/(1+q*q*XI*XI)
R=A[...,None]*(T[...,None]*aa+(1-T[...,None])*p)+(1-A[...,None])*cc[...,None]*dd
ell=1-np.sum(R*Y,axis=-1)
c0=(b*b-1)/72;c2=(1-b)*(9*b-8*Y[1])/144
F4=-q*q*c0*(1-b)+c2*(1+b)
weight=(2*wth[:,None,None])*wa[None,:,None]*(wx*dx*2*SX)[None,None,:]
got=F4/np.pi*np.sum(weight*XI/((1+q*q*XI*XI)*ell**2))
h=np.sqrt((1+b)/(1-b));K=-c0+c2*h*h/(q*q)
expected=K*np.arctan(q/h)**2/np.sqrt(aa0*pp0)
test('four_parameter_actual_source_chain',got-expected,3e-13)
facts['fixed_source_chain']={'Y':Y.tolist(),'q':q,'four_expected':float(expected),'four_got':float(got)}

# 4. Verify the exported actual cut rows by differentiating the elementary
# zonal primitives directly in R1,R2, including mixed derivatives.
cut=json.loads((OUT/'x24_cut.json').read_text())
R1S,R2S,uS,vS,qS=S.symbols('R1 R2 u v q');loc={'n':nS,'z':zS,'R1':R1S,'R2':R2S,'u':uS,'v':vS,'q':qS}
rv1,rv2,uv,vv,qq=m.mpf('.1'),m.mpf('.9'),m.mpf('.2'),m.mpf('1.2'),m.mpf('.5')
prim={0:'inverse2',1:'inverse1',2:'log',3:'ell_log'}
for idx in (0,1,4,7):
 item=cut['terms'][idx]
 pol=S.Poly(S.sympify(item['four_numerator'],locals={'q':qS,'Y1':S.Symbol('Y1'),'Y2':S.Symbol('Y2')}),S.Symbol('Y1'),S.Symbol('Y2'))
 def output(tv):
  total=m.mpc(0)
  for (a,b),co in pol.terms():
   d=a+b;fac=-1 if d==2 else 1
   def fn(r1,r2):return inv(prim[d],r1*r1+r2*r2,r1*uv+r2*tv)
   coeff=m.mpf(str(S.N(co.subs(qS,S.Rational(1,2)),60)))
   total+=coeff*fac*m.diff(fn,(rv1,rv2),(a,b))
  return total
 got=(output(vv+1j*eps)-output(vv-1j*eps))/(-2j*m.pi)
 row=S.lambdify((nS,zS,R1S,R2S,uS,vS,qS),S.Matrix([[S.sympify(x,locals=loc) for x in item['cut_row']]]),'mpmath')
 nn=rv1*rv1+rv2*rv2;zz=rv1*uv+rv2*vv
 pred=(row(nn,zz,rv1,rv2,uv,vv,qq)*m.matrix([1,m.sqrt((1-nn)/(zz*zz-nn))]))[0]
 test('actual_numerator_cut_term_'+str(idx),got-pred,m.mpf('2e-18'))

facts['checks']=checks;facts['elapsed_seconds']=time.time()-START
rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;facts['peak_rss_bytes']=rss if sys.platform=='darwin' else rss*1024
(OUT/'x24_numeric.json').write_text(json.dumps(facts,indent=2)+'\n')
print(json.dumps({'passed':len(checks),'elapsed_seconds':facts['elapsed_seconds'],'peak_rss_bytes':facts['peak_rss_bytes']},indent=2))
