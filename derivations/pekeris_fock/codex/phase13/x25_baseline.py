import hashlib, json
from pathlib import Path
root=Path(__file__).resolve().parents[2]
out=Path(__file__).resolve().parent
names=['COLLAB.md','THEORY.md','codex/TASK_codex_25.md',
       'claude10/op_e2_Q.sobj','claude10/op_e2_A.txt',
       'codex/phase9/x19_atoms.jl','codex/phase8/x16_connection.txt']
records={}
for name in names:
    p=root/name
    data=p.read_bytes()
    records[name]={'size':len(data),'sha256':hashlib.sha256(data).hexdigest()}
    if name in ['COLLAB.md','THEORY.md']:
        (out/('x25_prefix_'+name.replace('.md','.txt'))).write_bytes(data)
(out/'x25_baseline.json').write_text(json.dumps(records,indent=2)+'\n')
print('Recorded',len(records),'input hashes and two shared-log prefixes.')
