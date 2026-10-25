from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parents[2]
base=json.loads((root/'codex/phase13/x27_baseline.json').read_text())
out={}
for p,d in base.items():
 b=(root/p).read_bytes()
 if p in ('COLLAB.md','THEORY.md'):
  ok=len(b)>=d['size'] and hashlib.sha256(b[:d['size']]).hexdigest()==d['sha256']
  out[p]={'prefix_preserved':ok,'appended_bytes':len(b)-d['size']}
 else:
  ok=hashlib.sha256(b).hexdigest()==d['sha256'];out[p]={'unchanged':ok}
 assert ok,(p,out[p])
exact=json.loads((root/'codex/phase13/x27_exact.json').read_text())
assert len(exact['checks'])==7 and all(exact['checks'].values())
quad=json.loads((root/'codex/phase13/x27_quad.json').read_text())
assert float(quad['absolute_error'])<1e-60
saved=(root/'codex/phase13/x27_saved.out').read_text()
assert 'length 80 precision_first 256' in saved
rss=int(saved.split('maxrss_bytes ')[-1].strip())
assert rss<2_000_000_000
assert (root/'codex/phase13/x27_saved.err').stat().st_size==0
out['verification']={'exact_checks':7,'independent_test_quadrature':'PASS','saved_input_inspection':'PASS','julia_peak_rss_bytes':rss,'scientific_result':'OPEN; explicit missing-kernel/framing audit, not a B_A identification'}
(root/'codex/phase13/x27_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
