from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parents[2]
paths=['COLLAB.md','THEORY.md','claude11/chain.jl','claude11/lad_B1.jls','claude11/lad_B1_low.jls','claude11/op_A3.txt','claude11/fac_A3.sage','claude11/NOTES.md','claude11/claim_chain.md']
d={}
for p in paths:
 b=(root/p).read_bytes();d[p]={'size':len(b),'sha256':hashlib.sha256(b).hexdigest()}
(root/'codex/phase13/x27_baseline.json').write_text(json.dumps(d,indent=2)+'\n')
