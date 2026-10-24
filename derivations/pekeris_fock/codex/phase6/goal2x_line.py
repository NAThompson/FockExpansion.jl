import sys,json,time,argparse
sys.dont_write_bytecode=True
from goal2x_probe import *
from goal2x_rational import fit

def split(t,q):
 p=PRIME;ii=int(sp.sqrt_mod(-1,p));rt=int(sp.sqrt_mod(2,p))
 ca=(t+pow(t,-1,p))*pow(2,-1,p)%p;sa=(t-pow(t,-1,p))*pow(2*ii,-1,p)%p
 r=(ca-sa)*pow(rt,-1,p)%p;rr=(ca+sa)*pow(rt,-1,p)%p
 c=rt*(q+pow(q,-1,p))*pow(2,-1,p)%p;pp=-rt*(q-pow(q,-1,p))*pow(2*ii,-1,p)%p
 D=((1-c*c)**2-4*r*r*rr*rr)%p;disc=sp.sqrt_mod(D,p)
 if disc is None:return False
 try:
  v1=(1-c*c+int(disc))*pow(2*rr*rr,-1,p)%p
  v3=(rr*rr-r*r+int(disc))*pow(pp*pp,-1,p)%p
 except ValueError:return False
 return all(v and pow(v,(p-1)//2,p)==1 for v in (v1,v3))

def analyze(data,quick=False):
 samples=data['samples'];n=len(samples)-5
 if n<8:return
 xs=[x['u'] for x in samples];summary=[];hist={};passed=0;nonzero=0
 for a in range(10):
  row=[]
  for b in (range(45) if not quick else (0,1,5,44)):
   ys=[x['relations'][a][b] for x in samples]
   r=fit(xs[:n],ys[:n],xs[n:],ys[n:]);row.append(r)
   if any(ys):nonzero+=1;passed+=int(r['heldout']);key=str((r['degrees'][0],r['degrees'][1],r['heldout']));hist[key]=hist.get(key,0)+1
  summary.append(row)
 data['analysis']=dict(training=n,holdouts=5,quick=quick,nonzero=nonzero,nonzero_heldout_pass=passed,degree_histogram=hist,coefficients=summary)
 print('ANALYSIS',n,passed,'/',nonzero,hist,flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--count',type=int,default=37);ap.add_argument('--quick',action='store_true');args=ap.parse_args()
 path=OUT/'goal2x_line.json'
 data=json.loads(path.read_text()) if path.exists() else dict(prime=PRIME,pi=103,line='t=10+u, q=17+2u',samples=[],skips=[])
 u=data['samples'][-1]['u']+1 if data['samples'] else 0
 while len(data['samples'])<args.count:
  t=(10+u)%PRIME;q=(17+2*u)%PRIME
  if not split(t,q):u+=1;continue
  start=time.time()
  try:r=relations(t,q)
  except (ValueError,ZeroDivisionError) as e:
   data['skips'].append([u,str(e)]);u+=1;continue
  data['samples'].append(dict(u=u,point=[t,q],relations=r))
  path.write_text(json.dumps(data)+'\n')
  print('SAMPLE',len(data['samples']),u,'seconds',round(time.time()-start,3),flush=True)
  u+=1
 analyze(data,args.quick);path.write_text(json.dumps(data)+'\n')
