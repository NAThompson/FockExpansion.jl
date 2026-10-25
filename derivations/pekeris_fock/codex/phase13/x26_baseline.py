import hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent;r=p.parents[1]
d={str(f.relative_to(r)):hashlib.sha256(f.read_bytes()).hexdigest() for f in r.rglob('*') if f.is_file() and not str(f.relative_to(r)).startswith('codex/phase13/x26_') and str(f.relative_to(r)) not in ('COLLAB.md','THEORY.md')}
for n in ('COLLAB.md','THEORY.md'):
 b=(r/n).read_bytes();d[n]={'length':len(b),'sha256':hashlib.sha256(b).hexdigest()}
(p/'x26_baseline.json').write_text(json.dumps(d,indent=2))
