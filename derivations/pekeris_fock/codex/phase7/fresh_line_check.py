"""Fresh exact jets and independent elimination for saved per-panel line fits."""
from pathlib import Path
import json,sys
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(OUT.parent/'phase5'))
import goal2_combined_exact as G
A=G.A;A.N=7
A.IJ=[(i,k-i) for k in range(A.N+1) for i in range(k+1)];A.IX={v:k for k,v in enumerate(A.IJ)}
A.PAIRS=[[(A.IX[a,b],A.IX[i-a,j-b]) for a in range(i+1) for b in range(j+1)] for i,j in A.IJ]
P=G.P

def echelon(rows):
 out=[];n=len(rows)
 for k,row in enumerate(rows):
  v=row[:];rep=[int(j==k) for j in range(n)]
  for pivot,rr,cc in out:
   f=v[pivot]
   if f:
    v=[(x-f*y)%P for x,y in zip(v,rr)];rep=[(x-f*y)%P for x,y in zip(rep,cc)]
  pivot=next(i for i,x in enumerate(v) if x);scale=pow(v[pivot],-1,P)
  out.append((pivot,[x*scale%P for x in v],[x*scale%P for x in rep]))
 return out

def solve(B,row):
 v=row[:];out=[0]*len(B)
 for pivot,rr,cc in B:
  f=v[pivot]
  if f:
   v=[(x-f*y)%P for x,y in zip(v,rr)];out=[(x+f*y)%P for x,y in zip(out,cc)]
 assert not any(v),'target not in candidate span'
 return out

def ev(a,x):
 r=0
 for c in reversed(a):r=(r*x+c)%P
 return r

def rows_at(t,q,meta):
 _,C,T,pieces=G.build(t,q,True);out={}
 for piece in pieces:
  name=piece['name']
  if name=='classical':continue
  keys,targets=meta[name];wanted={tuple(k) for k in keys+targets};D=max(a+b for _,a,b in wanted);prev={(0,0):piece['rows']};saved={}
  for order in range(D+1):
   G.WORK=A.N-order;cur={}
   for a in range(order+1):
    b=order-a
    if not order:rows=piece['rows']
    else:
     axis=0 if a else 1;old=prev[a-1,b] if a else prev[a,b-1];rows=[]
     for row in old:
      o={i:v.diff(axis) for i,v in row.items()}
      for i,v in row.items():
       for j,w in C[axis].get(i,[]):o[j]=o.get(j,G.J())+v*w
      rows.append(o)
    cur[a,b]=rows
    for src,row in enumerate(rows):
     if (src,a,b) in wanted:saved[src,a,b]=G.project(row,T)
   prev=cur
  B=echelon([saved[tuple(k)] for k in keys]);out[name]=[solve(B,saved[tuple(k)]) for k in targets]
 return out

def run():
 data=json.loads((OUT/'line_inputs.json').read_text());fits=json.loads((OUT/'line_audit.json').read_text());results=[]
 used=set(data['xs']);checks=0
 for lam in range(2001,2100):
  if lam in used:continue
  try:out=rows_at(3+lam,4+2*lam,data['meta'])
  except (ValueError,ZeroDivisionError):continue
  for name,actual in out.items():
   block=fits['blocks'][name];pred=[[0]*len(data['meta'][name][0]) for _ in actual]
   for f in block['fits']:
    den=ev(f['denominator'],lam);assert den
    pred[f['target']][f['column']]=ev(f['numerator'],lam)*pow(den,-1,P)%P
   assert actual==pred,(lam,name)
   checks+=sum(len(row) for row in actual)
  results.append(lam);print('PASS fresh exact per-panel relations at lambda',lam,flush=True)
  if len(results)==2:break
 assert len(results)==2
 (OUT/'fresh_line_check.json').write_text(json.dumps(dict(prime=P,lambda_values=results,coefficients_checked=checks,method='own Taylor rows and independent modular elimination; every 240-state coordinate checked'),indent=2)+'\n')
 print('PASS',checks,'fresh coefficients',flush=True)
if __name__=='__main__':run()
