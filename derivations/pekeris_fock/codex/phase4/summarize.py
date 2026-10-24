from pathlib import Path
from decimal import Decimal as D
import re,json,hashlib
root=Path(__file__).resolve().parent
out={}
for h in (8,12):
 txt=(root/f'run_h{h}.out').read_text()
 rs=re.findall(r'RESULT N=(\d+) h=\d+ bits=\d+ dim=(\d+) E=([^ ]+) error=([^ ]+) a21=([^ ]+) residual=([^\n]+)',txt)
 assert len(rs)==3
 out[str(h)]={n:dict(dim=int(dim),E=E,error=err,a21=a,residual=res) for n,dim,E,err,a,res in rs}
 for n,row in out[str(h)].items():
  assert D(row['error'])>0 and D(row['residual'])<D('1e-45')
 assert all(D(out[str(h)][str(n)]['E'])>D(out[str(h)][str(n+1)]['E']) for n in (3,4))
for n in ('3','4','5'):
 de=abs(D(out['8'][n]['E'])-D(out['12'][n]['E']));da=abs(D(out['8'][n]['a21'])-D(out['12'][n]['a21']))
 assert de<D('1.2e-20') and da<D('1.1e-18')
 print('N',n,'energy refinement',de,'a21 refinement',da)
assert 'PASS seed' in (root/'check_seed.out').read_text()
assert 'PASS manufactured' in (root/'verify_jet.out').read_text()
for filename in ('corrected.out','global_state.out'):
 txt=(root/filename).read_text();assert 'ERROR:' not in txt
for p in (root/'input_states').glob('*.txt'):
 hashes=json.loads((root/'input_hashes.json').read_text());assert hashlib.sha256(p.read_bytes()).hexdigest()==hashes[str(p.relative_to(root.parent.parent))]
(root/'numerical_summary.json').write_text(json.dumps(out,indent=2)+'\n')
print('PASS nested Ritz energies, positive defects, residuals, refinement, seed and manufactured checks, input-state hashes')
for filename in ('review_state_clean.out','review_state_pole_removed.out'):
 txt=(root/filename).read_text();assert 'ERROR:' not in txt
 vals=re.findall(r'hinv=12 E-E0=([^ ]+) E-Eexport=([^\n]+)',txt)
 assert len(vals)==2
 assert all(abs(D(err))<D('9e-26') for _,err in vals)
print('PASS clean and pole-removed N4/N8 cutoff-integral comparisons')
def globals_from(filename):
 return {(int(n),int(q)):D(a) for n,q,a in re.findall(r'GLOBAL Claude N=(\d+) nq=(\d+) a21=([^ ]+)',(root/filename).read_text())}
a=globals_from('global_state.out');b=globals_from('global_state_pole_removed.out')
assert len(a)==8 and len(b)==6
for k in b:assert abs(a[k]-b[k])<D('1e-30')
for n in (4,8,10,11):assert abs(a[(n,24)]-a[(n,48)])<D('2.8e-15')
assert 'NONZERO polynomial' in (root/'pole_audit.out').read_text()
print('PASS global quadrature comparison, exact pole detection, and pole-removed global estimates unchanged within 1e-30')
