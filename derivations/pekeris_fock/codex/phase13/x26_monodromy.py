"""Independent based-loop experiment; NOT an identification of the X-strip path.
All coefficients read from plain-text rational operators. Serial 40-digit Taylor stepping.
"""
import os,json,time,resource
from pathlib import Path
import mpmath as m
m.mp.dps=int(os.environ.get('X26_DPS','40')); ORDER=int(os.environ.get('X26_ORDER','65'))
p=Path(__file__).resolve().parent;r=p.parents[1]; start=time.time()
def rd(x):
 a=x.split('/');return m.mpf(a[0])/(m.mpf(a[1]) if len(a)>1 else 1)
def load(n):return [[rd(x) for x in row.split()] for row in (r/'claude10'/n).read_text().splitlines()]
fall=[[m.mpf(m.factorial(n)/m.factorial(n-j)) if n>=j else 0 for j in range(5+1)] for n in range(ORDER+6)]
def shift(c,z):return [sum(c[i]*m.binomial(i,j)*z**(i-j) for i in range(j,len(c))) for j in range(len(c))]
def step(ops,z,h,Y):
 sh=[shift(c,z) for c in ops]
 aa=[[Y[j,k]/m.factorial(j) for k in range(5)] for j in range(5)]
 for n in range(ORDER-5):
  terms=[]
  for j,row in enumerate(sh):
   for d,co in enumerate(row):
    ix=n+j-d
    if ix<j or (j==5 and d==0):continue
    if co:terms.append((co*fall[ix][j],ix))
  den=sh[5][0]*fall[n+5][5]
  aa.append([-sum(co*aa[ix][k] for co,ix in terms)/den for k in range(5)])
 out=m.matrix(5,5)
 for j in range(5):
  weights=[fall[n][j]*h**(n-j) for n in range(j,ORDER)]
  for k in range(5):out[j,k]=sum(w*aa[n][k] for n,w in zip(range(j,ORDER),weights))
 return out
count=0
def transport(ops,path,sing):
 global count
 Y=m.eye(5);z=m.mpc(path[0])
 for goal in path[1:]:
  goal=m.mpc(goal)
  while abs(goal-z)>m.mpf('1e-35'):
   dist=min(abs(z-q) for q in sing);h=goal-z
   if abs(h)>m.mpf('.32')*dist:h*=m.mpf('.32')*dist/abs(h)
   Y=step(ops,z,h,Y);z+=h;count+=1
 return Y
def loop(ops,base,point,sing,radius,above=False):
 st=point+1j*radius
 stem=[base,base+1j*m.mpf('.25'),point+1j*m.mpf('.25'),st] if above else [base,point-radius if base<point else point+radius]
 st=stem[-1]
 circ=[point+(st-point)*m.exp(2j*m.pi*k/20) for k in range(1,21)]
 return transport(ops,stem+circ+list(reversed(stem[:-1])),sing)
def sv(A):return [m.nstr(x,8) for x in m.svd(A,compute_uv=False)]
A=load('op_e2_A.txt');C=load('op_eC.txt');b1=m.mpf(1)/5;b2=m.mpf(1)/3
la=m.sqrt(m.mpf(8)/9);lc=m.sqrt(m.mpf(24)/25)
sa=[la,-la,1,-1]
# Include all apparent roots, found numerically from exact leading polynomial.
# Known sextic roots from the independently supplied C operator.
sc=[lc,-lc,1,-1,m.mpf('1.54950727000600'),-m.mpf('1.54950727000600')]+[m.mpc(i*m.mpf('.965152178692833'),j*m.mpf('.121910872334630')) for i in (-1,1) for j in (-1,1)]
J={tuple(map(int,line.split()[:2])):[rd(x) for x in line.split()[2:]] for line in (r/'claude10/jets_Q.txt').read_text().splitlines()}
KC=m.matrix([J[(0,j)] for j in range(5)])
results={}; mats={}
for name,op,base,pt,sing,rad,above,K in [
 ('Q2',A,b1,-la,sa,m.mpf('.02'),False,m.eye(5)),
 ('R_u_minus',A,b1,-1,sa,m.mpf('.01'),True,m.eye(5)),
 ('U_v_plus',C,b2,1,sc,m.mpf('.004'),True,KC)]:
 T=K**-1*loop(op,base,pt,sing,rad,above)*K;mats[name]=T
 results[name]={'rank_singular_values':sv(T-m.eye(5)), 'involution_residual':m.nstr(m.norm(T*T-m.eye(5)),8),'nilpotent_square_residual':m.nstr(m.norm((T-m.eye(5))**2),8)}
 print(name,results[name],flush=True)
 (p/'x26_monodromy_partial.json').write_text(json.dumps(results,indent=2))
N=mats['Q2']-m.eye(5)
for name,G in [('R',mats['R_u_minus']),('U',mats['U_v_plus']),('RU',mats['R_u_minus']*mats['U_v_plus']),('UR',mats['U_v_plus']*mats['R_u_minus'])]:
 B=G*N*G**-1
 results['conjugation_'+name]={'stacked_images':sv(m.matrix([[N[i,j] for j in range(5)]+[B[i,j] for j in range(5)] for i in range(5)])),'product':sv(N*B),'triple':sv(B*N*B)}
 print(name,results['conjugation_'+name],flush=True)
results['matrices']={k:[[[m.nstr(m.re(T[i,j]),m.mp.dps),m.nstr(m.im(T[i,j]),m.mp.dps)] for j in range(5)] for i in range(5)] for k,T in mats.items()}
results['seconds']=time.time()-start;results['steps']=count;results['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
(p/'x26_monodromy.json').write_text(json.dumps(results,indent=2))
