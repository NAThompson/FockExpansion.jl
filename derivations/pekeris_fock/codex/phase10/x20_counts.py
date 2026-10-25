"""Serial extension-field counts and reciprocal-polynomial reconstruction."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import sympy as s
import subprocess,json,time,resource,argparse
ap=argparse.ArgumentParser();ap.add_argument("--swap",action="store_true");args=ap.parse_args()
out=Path(__file__).resolve().parent;t=s.Symbol('t');T=s.Symbol('T');start=time.time();data=[]
for p in [11,13,19]:
 rows=[]
 for n in [1,2,3]:
  rel=[0,0,0]
  if n>1:
   for c0 in range(1,p):
    ok=False
    for c1 in range(p):
     f=s.Poly(t**n-c1*t-c0,t,modulus=p)
     if f.is_irreducible:
      rel=[c0,c1,0];ok=True;break
    if ok:break
  u=pow(4 if args.swap else 3,-1,p);v=pow(3 if args.swap else 4,-1,p)
  cmd=[str(out/'x20_count'),str(p),str(n),*map(str,rel),str(u),str(v)]
  row=json.loads(subprocess.check_output(cmd,text=True,timeout=240))
  # Three independent sections: quadratic character plus the two zeta8 characters.
  h=(-2*(1-u*u))%p
  eps0=1 if pow(h,(p-1)//2,p)==1 else -1
  pm=pow(p,n,8)
  tracepair=2 if pm==1 else -2 if pm==5 else 0
  row['MW_trace']=eps0**n+tracepair
  row['T5_trace']=row['char_sum']-row['MW_trace']*p**n
  row['irreducible_relation']=rel
  rows.append(row);print(json.dumps(row),flush=True)
 # L(T)=det(1-Frob*T) on the 5-dimensional orthogonal complement of known 17 divisors.
 a=[s.Integer(1)]
 for k in range(1,4):
  a.append(-sum(a[k-i]*rows[i-1]['T5_trace'] for i in range(1,k+1))/k)
 candidates=[]
 for sign in [-1,1]:
  # normalized degree5 polynomial has reciprocity a[5-k]=sign*p**(5-2k)*a[k].
  if a[3]!=sign*p*a[2]:continue
  aa=a+[sign*p**3*a[1],sign*p**5]
  pol=sum(aa[k]*T**k for k in range(6))
  candidates.append({'sign':sign,'polynomial':str(pol),'factorization':str(s.factor(pol))})
 result={'p':p,'rows':rows,'reciprocal_candidates':candidates};data.append(result)
 print('CANDIDATES',p,candidates,flush=True)
(out/('x20_counts_line.json' if args.swap else 'x20_counts.json')).write_text(json.dumps({'data':data,'elapsed_seconds':time.time()-start,'peak_child_rss_bytes':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,'peak_python_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2)+'\n')
