"""Read-only protected-input/prefix audit plus saved-output checks."""
from pathlib import Path
import hashlib,json,re
D=Path(__file__).resolve().parent
R=D.parent.parent
before=json.loads((D/'x30_inputs.json').read_text())
checks={}
for name,b in before.items():
    raw=(R/name).read_bytes()
    candidate=raw[:b['size']] if name in ('COLLAB.md','THEORY.md') else raw
    checks[name]=hashlib.sha256(candidate).hexdigest()==b['sha256']
assert all(checks.values()),checks
control=(D/'x30_control.out').read_text()
assert float(re.search(r'max_abs=([^\n]+)',control)[1])<1e-62
fit=(D/'x30_fit.out').read_text()
assert 'rank=500' in fit and 'columns=604' in fit
for name in ('x30_eval.out','x30_fit.out','x30_rational.out','x30_denominator.out'):
    if (D/name).exists():
        text=(D/name).read_text()
        rss=list(map(int,re.findall(r'rss=(\d+)',text)))
        assert rss and max(rss)<1_500_000_000,(name,rss)
assert 'PASS exact cubic' in (D/'x30_exact.out').read_text()
artifacts={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(D.glob('x30_*')) if p.is_file() and p.name!='x30_audit_result.json'}
print(json.dumps(dict(status='PASS',protected_and_prefix_checks=checks,artifact_sha256=artifacts),indent=2))
