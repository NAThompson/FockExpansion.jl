import sys,json
sys.dont_write_bytecode=True
from goal2x_probe import *
from goal2x_rational import fit,ev

def main():
 d=json.loads((OUT/'goal2x_independent.json').read_text());result=[]
 for base in d['points']:
  point=base['point'];r=base['relations'];variations=[v for v in d['constant_results'] if v['point']==point]
  pis=[(103,r)]+[(v['pi'],v['relations']) for v in variations if 'pi' in v]
  xs=[v[0] for v in pis];fits=[];hist={}
  for a in range(10):
   row=[]
   for b in range(45):
    ys=[v[1][a][b] for v in pis]
    f=fit(xs[:4],ys[:4],xs[4:],ys[4:]);assert f['heldout'],(point,a,b,f)
    row.append(f);key=str(f['degrees']);hist[key]=hist.get(key,0)+1
   fits.append(row)
  parity=[]
  for v in variations:
   if 'si' not in v:continue
   rr=v['relations'];odd=[]
   for a in range(10):
    for b in range(45):
     if rr[a][b]!=r[a][b]:assert (rr[a][b]+r[a][b])%PRIME==0;odd.append([a,b])
   parity.append(dict(si=v['si'],sr=v['sr'],odd_entries=odd))
  result.append(dict(point=point,pi_fits=fits,pi_degree_histogram=hist,embedding_parity=parity))
  print('PASS constants',point,hist,flush=True)
 (OUT/'goal2x_constants.json').write_text(json.dumps(dict(prime=PRIME,results=result),indent=2)+'\n')
if __name__=='__main__':main()
