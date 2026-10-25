"""Residue-cycle jets versus exact Frobenius space (numerical, not balls)."""
import json, time, resource, sys
from pathlib import Path
import mpmath as m
out=Path(__file__).resolve().parent; start=time.time(); m.mp.dps=95
u=m.sqrt(8)/3; v=m.mpf(1)/3
a=m.sqrt((u+1)/v); b=m.sqrt((1-u)/v)
roots=[a,-a,1j*b,-1j*b]
def integrate(c,R,N):
    sums=[m.mpc(0) for _ in range(4)]
    for j in range(N):
        th=2*m.pi*j/N; e=m.exp(1j*th); q=c+R*e
        y=m.sqrt(u*v)
        for rt in roots:
            if abs(rt-c)<R:
                y*=m.sqrt(R)*m.exp(1j*th/2)*m.sqrt(1-(rt-c)/(R*e))
            else:
                y*=m.sqrt(c-rt)*m.sqrt(1+R*e/(c-rt))
        P=q**4-1; F=y*y; k=v*P/F; l=P*P/(F*F)
        # u=sqrt(8/9+t); delta=w^2=-t.
        h1=1/(2*u); h2=-1/(8*u**3); h3=1/(16*u**5)
        # F(t)^(-1/2) + (3/16)P^2*t*F(t)^(-5/2)
        # +(105/1024)P^4*t^2*F(t)^(-9/2)
        # +(1155/16384)P^6*t^3*F(t)^(-13/2).
        g0=1
        g1=-k*h1/2+3*l/16
        g2=-k*h2/2+3*k*k*h1*h1/8-15*l*k*h1/32+105*l*l/1024
        g3=(-k*h3/2+3*k*k*h1*h2/4-5*k**3*h1**3/16
            +3*l/16*(-5*k*h2/2+35*k*k*h1*h1/8)
            -945*l*l*k*h1/2048+1155*l**3/16384)
        dq=1j*R*e
        for ix,g in enumerate([g0,g1,g2,g3]): sums[ix]+=dq*g/y
    return [s*2*m.pi/N*2*m.pi*1j for s in sums]
def fmt(x): return [m.nstr(x.real,85),m.nstr(x.imag,85)]
records={}; checks={}; jets=[]
circles=[('A',m.mpc(0),m.mpf(1),[256,512]),
         ('B',(a+1j*b)/2,None,[4096,8192])]
for name,c,R,ns in circles:
    if R is None:
        ds=sorted(abs(rt-c) for rt in roots)
        R=m.sqrt(ds[1]*ds[2])
    prev=None
    for N in ns:
        vals=integrate(c,R,N)
        rec={'N':N,'center':fmt(c),'radius':m.nstr(R,85),'jets':list(map(fmt,vals))}
        if prev:
            err=max(abs(x-y)/max(1,abs(y)) for x,y in zip(vals,prev))
            rec['relative_refinement_error']=str(err)
            checks[name+'_30_digit_refinement']=err<m.mpf('1e-35')
        constraint=23571*vals[0]+195840*vals[1]-65536*vals[2]
        rec['log_constraint_scaled_residual']=str(abs(constraint)/max(1,*map(abs,vals))/65536)
        predicted=(8505*vals[0]+373590*vals[1]+550912*vals[2])/131072
        rec['ODE_first_recurrence_scaled_residual']=str(abs(vals[3]-predicted)/max(1,abs(predicted)))
        records[name+'_'+str(N)]=rec; prev=vals
        print(name,N,'log constraint',rec['log_constraint_scaled_residual'],
              'ODE',rec['ODE_first_recurrence_scaled_residual'],flush=True)
    checks[name+'_log_constraint_30_digits']=abs(constraint)/max(1,*map(abs,vals))/65536<m.mpf('1e-35')
    checks[name+'_ODE_recurrence_30_digits']=abs(vals[3]-predicted)/max(1,abs(predicted))<m.mpf('1e-35')
    jets.append(vals)
det=jets[0][0]*jets[1][1]-jets[0][1]*jets[1][0]
checks['independent_normal_jets']=abs(det)>m.mpf('1e-5')
lam=(u-1)/(u+1)
legendre=[2*m.pi*1j*4*m.ellipk(lam)/(a*m.sqrt(u*v)),
          2*m.pi*1j*2j*m.ellipk(1-lam)/(a*m.sqrt(u*v))]
legendre_errors=[abs(jet[0]-value)/abs(value) for jet,value in zip(jets,legendre)]
for name,error in zip(['A','B'],legendre_errors):
    checks[name+'_Legendre_30_digits']=error<m.mpf('1e-35')
assert all(checks.values()),checks
res={'precision_digits':m.mp.dps,'u0':str(u),'v':str(v),'roots':list(map(fmt,roots)),
     'records':records,'two_by_two_jet_determinant':fmt(det),'checks':checks,
     'Legendre_values':list(map(fmt,legendre)),
     'Legendre_relative_errors':list(map(str,legendre_errors)),
     'scope':'Independent contour quadrature and first normal recurrence; not a rigorous enclosure or PF reconstruction.',
     'elapsed_seconds':time.time()-start,
     'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)}
(out/'x25_periods.json').write_text(json.dumps(res,indent=2)+'\n')
print('All',len(checks),'checks pass; seconds',res['elapsed_seconds'])
