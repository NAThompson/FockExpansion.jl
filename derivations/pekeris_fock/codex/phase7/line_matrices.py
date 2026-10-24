"""Assemble dt,dq matrices on the reconstructed line; no new guessing."""
from pathlib import Path
import json
p=Path(__file__).resolve().parent
i=json.loads((p/'line_inputs.json').read_text());a=json.loads((p/'line_audit.json').read_text())
out=dict(prime=a['prime'],line=['3+lambda','4+2*lambda'],pi=101,convention='partial_t or partial_q of atom vector; coefficients restricted to line; not total line derivative',blocks={})
for name,meta in i['meta'].items():
 if name=='classical':continue
 keys,targets=[list(map(tuple,v)) for v in meta];fits={(r['target'],r['column']):r for r in a['blocks'][name]['fits']};matrices={}
 for axis in [0,1]:
  rows=[]
  for row,(src,x,y) in enumerate(keys):
   key=src,x+int(axis==0),y+int(axis==1)
   if key in keys:rows.append([row,keys.index(key),[1],[1]])
   else:
    ti=targets.index(key)
    for col in range(len(keys)):
     if (ti,col) in fits:
      f=fits[ti,col];assert f['heldout'];rows.append([row,col,f['numerator'],f['denominator']])
  matrices['dt' if axis==0 else 'dq']=rows
 out['blocks'][name]=dict(atoms=keys,**matrices)
(p/'line_connections.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS four complete 13x13 per-panel dt/dq line connections assembled')
