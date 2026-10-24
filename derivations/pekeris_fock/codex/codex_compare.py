"""Read-only comparison, run only after independent artifacts were produced."""
import json, hashlib, re, sympy as s
from pathlib import Path
from codex_derive import eps,E,Z,lam,l,m,n
base=Path(__file__).resolve().parent
ref=base.parent/'sym'/'recurrence.json'
ours=json.loads((base/'codex_recurrence.json').read_text()); theirs=json.loads(ref.read_text())
local={'lam':lam,'eps':eps,'epsilon':eps,'E':E,'Z':Z,'lambda':lam,'l':l,'m':m,'n':n}
results={}; diffs={}
for own,other,fixed in [('P','P_general',False),('P','P',True),('S2P','sigma2P',True),('C1','C1',True),('C2','C2',True)]:
    a={k:s.expand(s.sympify(re.sub(r'\blambda\b','lam',v),locals=local).subs(E,-eps**2) if fixed else s.sympify(re.sub(r'\blambda\b','lam',v),locals=local)) for k,v in ours[own].items()}
    a={k:v for k,v in a.items() if v!=0}
    b={k:s.sympify(re.sub(r'\blambda\b','lam',v),locals=local) for k,v in theirs[other].items()}
    delta={k:str(s.factor(a.get(k,0)-b.get(k,0))) for k in a.keys()|b.keys() if s.expand(a.get(k,0)-b.get(k,0))!=0}
    label=own+('_fixed' if fixed else '_general')
    results[label]={'ours_count':len(a),'reference_count':len(b),'union_count':len(a.keys()|b.keys()),'nonzero_differences':len(delta)}
    diffs[label]=delta
summary={'reference_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),'comparisons':results,'differences':diffs}
(base/'codex_comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
