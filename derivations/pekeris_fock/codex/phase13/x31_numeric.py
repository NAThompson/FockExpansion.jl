"""Independent nonsingular cap quadrature (mpmath), stdout only."""
import mpmath as m,resource
m.mp.dps=65
pi=m.pi
# s=b-X, t=b-Y, r=b-eta. This exact substitution removes both cap radicals.
def cap(s,t,f):
 D=m.sin(s)*m.sin(t);R=m.sqrt(1+D)
 center=2*m.atan2(m.sin((s+t)/2),m.cos((s-t)/2))
 def fun(th):
  r=center+2*m.asin(m.sqrt(D)/R*m.sin(th))
  return (1+4*D*m.cos(th)**2)*f(r)/pi
 return m.quad(fun,[-pi/2,0,pi/2])
modes={'sin X':lambda r:-m.cos(r),'1+sqrt2 cos X':lambda r:1-m.sqrt(2)*m.sin(r),'1-sqrt2 cos X':lambda r:1+m.sqrt(2)*m.sin(r)}
worst=m.mpf(0)
for s,t in [(m.mpf('.7'),m.mpf('1.1')),(m.mpf('2.4'),m.mpf('.8')),(m.mpf('1.3'),m.mpf('1.3')),(m.mpf('.9'),pi-m.mpf('.9'))]:
 for name,f in modes.items():
  err=abs(cap(s,t,f)-f(s)*f(t)/f(0));worst=max(worst,err)
  print('mode',name,'s,t=',m.nstr(s,6),m.nstr(t,6),'error',m.nstr(err,8))
 assert abs(cap(s,t,lambda r:1)-(1+2*m.sin(s)*m.sin(t)))<m.mpf('1e-60')
print('worst exact-mode error',m.nstr(worst,10));assert worst<m.mpf('1e-60')
f=lambda r:m.exp(r/7)+r*r/13
s=m.mpf('1.1')
for eps in map(m.mpf,['1e-3','1e-6','1e-9','1e-12']):
 print('endpoint t=',eps,'delta error',m.nstr(abs(cap(s,eps,f)-f(s)),12),'reflected error',m.nstr(abs(cap(s,pi-eps,f)-f(pi-s)),12))
# Weak boundary flux for a general analytic test, differentiating the regularized quadrature.
for eps in map(m.mpf,['1e-3','1e-6','1e-9']):
 flux=m.sin(eps)*m.diff(lambda t:cap(s,t,f),eps)
 print('weak endpoint sin(t)*d_t action',eps,m.nstr(flux,12))
# Distributional M_s=M_t, tested against a non-eigenfunction, including cap contacts.
def L(v):
 return lambda r:-m.sin(r)*m.diff(v,r,2)-m.cos(r)*m.diff(v,r)-2*m.sin(r)*v(r)
for s,t in [(m.mpf('.7'),m.mpf('1.1')),(m.mpf('1.3'),m.mpf('1.3')),(m.mpf('.9'),pi-m.mpf('.9'))]:
 Ms=-m.sin(s)*m.diff(lambda ss:cap(ss,t,f),s,2)-m.cos(s)*m.diff(lambda ss:cap(ss,t,f),s)-2*m.sin(s)*cap(s,t,f)
 err=abs(Ms-cap(s,t,L(f)));print('weak operator identity at',m.nstr(s,6),m.nstr(t,6),m.nstr(err,10));assert err<m.mpf('1e-58')
print('ru_maxrss bytes',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
print('ALL NUMERICAL CHECKS PASSED')
