import json
from pathlib import Path
import goal2_combined_exact as B
A=B.A;J=B.J

def run(t,q,D=8):
 L,C,T,pieces=B.build(t,q,True);B.flatcheck(C);results=[]
 for piece in pieces:
  L=piece['rows'];prev={(0,0):L};bases=[[] for _ in L]+[[]];hist=[];selected=[[] for _ in bases]
  for order in range(D+1):
   B.WORK=A.N-order;current={}
   for a in range(order+1):
    b=order-a
    if not order:rows=L
    else:
     axis=0 if a else 1;old=prev[a-1,b] if a else prev[a,b-1];rows=[]
     for row in old:
      out={i:v.diff(axis) for i,v in row.items()}
      for i,v in row.items():
       for j,w in C[axis].get(i,[]):out[j]=out.get(j,J())+v*w
      rows.append(out)
    current[a,b]=rows
    for src,row in enumerate(rows):
     vec=B.project(row,T)
     if A.rankadd(bases[src],vec):selected[src].append([src,a,b])
     if A.rankadd(bases[-1],vec):selected[-1].append([src,a,b])
   prev=current;hist.append([len(bb) for bb in bases]);print(t,q,piece['name'],order,hist[-1],flush=True)
  results.append(dict(name=piece['name'],ranks=hist,selected=selected))
 return dict(point=[t,q],pieces=results)
if __name__=='__main__':
 prior=json.loads(Path(__file__).with_name('goal2_combined_exact.json').read_text());out=[]
 for result in prior['results']:
  out.append(run(*result['point']))
  Path(__file__).with_name('goal2_pieces_exact.json').write_text(json.dumps(dict(prime=B.P,pi=101,D=8,results=out),indent=2)+'\n')
