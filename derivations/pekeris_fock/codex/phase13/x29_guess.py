"""Serial multipoint derivative QR calibration; exact source is x29_exact.json."""
import json,time,resource,sys,math
from pathlib import Path
from fractions import Fraction as F
import mpmath as m
import numpy as np
out=Path(__file__).resolve().parent; data=json.loads((out/'x29_exact.json').read_text())
m.mp.dps=210; start=time.time()
def poly(c,x):
    v=m.mpf(0)
    for a in reversed(c): v=v*x+m.mpf(a)
    return v
def state(x):
    z=(x-1)/(x+1); kk=m.ellipk(z); ee=m.ellipe(z); root=m.sqrt(-x*(x+1))
    p=kk/root
    pp=(ee/(2*z*(1-z))-kk/(2*z))*2/(x+1)**2/root-(2*x+1)/(2*x*(x+1))*p
    l=m.log(1+x)
    t=m.log(1-x*x)+m.log(-x)+l*l/4+m.polylog(2,-x)
    return [p,pp,m.mpf(1),l,t]
def derivatives(x,r):
    v=state(x)
    return [sum(poly(a,x)/poly(b,x)*z for (a,b),z in zip(row,v)) for row in data['jets'][:r+1]]
def qrtest(r,d):
    n=(r+1)*(d+1); N=n+60
    points=[-m.mpf(19)/20+m.mpf(9)/10*i/(N-1) for i in range(N)]
    vals=[derivatives(x,r) for x in points]
    cols=[np.array([((x+m.mpf(1)/2)/(m.mpf(9)/20))**j*v[k] for x,v in zip(points,vals)],dtype=object) for k in range(r+1) for j in range(d+1)]
    scales=[m.sqrt(np.dot(c,c)) for c in cols]; cols=[c/s for c,s in zip(cols,scales)]
    perm=list(range(n)); R=[[m.mpf(0) for j in range(n)] for i in range(n)]; diag=[]
    for k in range(n):
        j=max(range(k,n),key=lambda j:np.dot(cols[j],cols[j]))
        cols[k],cols[j]=cols[j],cols[k]; perm[k],perm[j]=perm[j],perm[k]
        for h in range(k): R[h][k],R[h][j]=R[h][j],R[h][k]
        norm=m.sqrt(np.dot(cols[k],cols[k])); R[k][k]=norm; diag.append(norm)
        q=cols[k]/norm
        for j in range(k+1,n):
            a=np.dot(q,cols[j]); cols[j]-=q*a
            b=np.dot(q,cols[j]); cols[j]-=q*b
            R[k][j]=a+b
    rec={'r':r,'d':d,'unknowns':n,'points':N,'last':str(diag[-1]),'previous':str(diag[-2]),'ratio':str(diag[-1]/diag[-2])}
    if (r,d)==(5,22):
        (out/'x29_samples.json').write_text(json.dumps({'dps':210,'rows':[{'u':str(x),'derivatives':list(map(str,v))} for x,v in zip(points,vals)]},indent=2)+'\n')
        z=[m.mpf(0)]*n; z[-1]=m.mpf(1)
        for i in reversed(range(n-1)): z[i]=-sum(R[i][j]*z[j] for j in range(i+1,n))/R[i][i]
        coeff=[m.mpf(0)]*n
        for i,j in enumerate(perm): coeff[j]=z[i]/scales[j]
        exact=[]
        for cs in data['coefficients']:
            cs=list(map(F,cs)); cc=[F(0)]*(d+1)
            for k,c in enumerate(cs):
                for j in range(k+1): cc[j]+=c*math.comb(k,j)*F(-1,2)**(k-j)*F(9,20)**j
            exact+=cc
        pivot=max(range(n),key=lambda j:abs(exact[j])); coeff=[c/coeff[pivot] for c in coeff]
        target=[q/exact[pivot] for q in exact]
        err=max(abs(c-m.mpf(t.numerator)/t.denominator) for c,t in zip(coeff,target))
        recovered=[F(str(c)).limit_denominator(10**40) for c in coeff]
        rec.update({'exact_nullvector_error':str(err),'rational_reconstruction_matches':recovered==target})
        assert recovered==target and err<m.mpf('1e-90')
    print(json.dumps(rec),flush=True)
    return rec
# Independent derivative checks by numerical differentiation, not connection rows.
def f(x):
    s=state(x);return (1+x)*s[0]+x*s[1]+s[4]
checks={}
for x in [m.mpf('-0.83'),m.mpf('-0.37')]:
    ds=derivatives(x,5)
    err=max(abs(ds[k]-m.diff(f,x,k))/max(1,abs(ds[k])) for k in range(6))
    checks[str(x)]=str(err); assert err<m.mpf('1e-190')
records=[qrtest(r,d) for r,d in [(4,22),(5,21),(5,22)]]
res={'dps':m.mp.dps,'records':records,'independent_derivative_errors':checks,'seconds':time.time()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)}
(out/'x29_guess.json').write_text(json.dumps(res,indent=2)+'\n');print('NUMERIC CHECKS PASS',res['seconds'],res['peak_rss_bytes'],flush=True)
