"""Independent Euclidean interpolation of frozen Claude per-block line samples."""
import json,sys,hashlib
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
sys.path.insert(0,str(OUT.parent/'phase6'))
import goal2x_rational as R
R.P=998244353

def run():
 snapshot=OUT/'line_inputs.json'
 if snapshot.exists():
  snap=json.loads(snapshot.read_text());xs=snap['xs'];meta=snap['meta'];combined=dict(zip(xs,snap['values']))
 else:
  combined={};meta=None;hashes={}
  for k in range(8):
   p=ROOT/'claude7'/f'bl_{k}.json';raw=p.read_bytes();hashes[str(p.relative_to(ROOT))]=hashlib.sha256(raw).hexdigest();d=json.loads(raw)
   if meta is None:meta=d['meta']
   assert d['meta']==meta
   for i,x in enumerate(d['xs']):
    values={n:d['vals'][n][i] for n in d['vals']}
    if x in combined:assert combined[x]==values
    combined[x]=values
  xs=sorted(combined);assert len(xs)>=160
  snap=dict(prime=R.P,line=['3+lambda','4+2*lambda'],pi=101,xs=xs,meta=meta,values=[combined[x] for x in xs],source_hashes=hashes)
  snapshot.write_text(json.dumps(snap)+'\n')
 result=dict(prime=R.P,training=128,holdouts=len(xs)-128,blocks={})
 for name in meta:
  keys,targets=meta[name];fits=[];hist={};bad=0;nonzero=0
  for ti in range(len(targets)):
   for j in range(len(keys)):
    ys=[combined[x][name][ti][j] for x in xs]
    if not any(ys):continue
    nonzero+=1;f=R.fit(xs[:128],ys[:128],xs[128:],ys[128:]);assert f
    if not f['heldout']:bad+=1
    else:hist[str(f['degrees'])]=hist.get(str(f['degrees']),0)+1
    fits.append(dict(target=ti,column=j,**f))
  result['blocks'][name]=dict(candidate_rank=len(keys),targets=len(targets),nonzero=nonzero,failed_holdouts=bad,degree_histogram=hist,fits=fits)
  print(name,'rank',len(keys),'nonzero',nonzero,'failed',bad,'degrees',hist,flush=True)
 (OUT/'line_audit.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':run()
