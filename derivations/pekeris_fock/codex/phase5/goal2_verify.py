"""Small, serial final audit; all output stays in the Codex allowlist."""
import json,re,hashlib
from pathlib import Path
root=Path(__file__).resolve().parent
for name,phrase in [('goal2_port_check.out','PASS: 312'),('goal2_classical_check.out','PASS 72'),('goal2_radicals.out','PASS z3^2'),('goal2_combined_check.out','PASS actual endpoint')]:
 text=(root/name).read_text();assert phrase in text and 'ERROR:' not in text and 'Traceback' not in text,name
 print('PASS saved verification',name)
classical=json.loads((root/'goal2_classical_exact.json').read_text())
assert len(classical['results'])==2
assert all(v['ranks']==[1,3,6,10,12,12,12,12,12] for v in classical['results'])
print('PASS exact classical two-point histories')
text=(root/'goal2_big_panel.out').read_text();rows=[l for l in text.splitlines() if 'order = 8,' in l]
assert len(rows)==6
for name,rank in [('N0',15),('N2',18),('joint',19)]:
 assert sum(f'output = "{name}"' in l and f'ranks = [{rank}, {rank}, {rank}]' in l for l in rows)==2
print('PASS high-precision panel two-point histories')
combined=json.loads((root/'goal2_combined_exact.json').read_text())
expected=[3,7,11,15,20,26,33,41]+[45]*7
assert len(combined['results'])==2
assert all(r['ranks']==expected and len(r['selected'])==45 for r in combined['results'])
assert combined['results'][0]['selected']==combined['results'][1]['selected']
assert (root/'goal2_combined_exact.out').read_text().count('PASS exact constant curvature')==2
print('PASS exact combined two-point stabilization through order 14 and common 45-row selection')
pieces=json.loads((root/'goal2_pieces_exact.json').read_text())
assert len(pieces['results'])==2
for point in pieces['results']:
 assert len(point['pieces'])==5
 for piece in point['pieces']:
  final=[12,12] if piece['name']=='classical' else [11,12,13]
  assert piece['ranks'][4:]==[final]*5
print('PASS all four panels and classical piece, two exact points, five stable orders')
(root/'goal2_candidate_generators.json').write_text(json.dumps(dict(
 qualification='Candidate in the endpoint quotient; actual-function minimality and generic closure unproved',
 coordinates=dict(t='(1-i)*(r1+i*r2)/(sqrt(2)*R)',q='(r12-i*P)/(sqrt(2)*R)'),
 sources=['1','chi modulo algebraic terms and a nonzero scale','R^3 psi30 modulo algebraic terms'],
 derivatives=combined['results'][0]['selected']),indent=2)+'\n')
inputs=[root.parent/'phase3'/n for n in ('panel_connection.py','classical_reference.py','structure.md','connection.md')]
inputs += [root.parent/'phase2'/n for n in ('order2.py','verify_atoms.py','minimal_rank.md')]
manifest={str(p.relative_to(root.parents[1])):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
(root/'goal2_input_hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('PASS input hashes saved; no external writes')
