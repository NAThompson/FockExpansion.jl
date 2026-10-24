"""Final artifact/inventory manifest; no writes outside phase8."""
from pathlib import Path
import json,hashlib
O=Path(__file__).resolve().parent;R=O.parents[1]
checks={'build.out':['H 225 [3, 4]','k 174 [3, 3]'],
'distances.out':['x 210 [7, 8]','y 230 [8, 9]','c 233 [7, 8]'],
'verify.out':['PASS all small serial checks','PASS physical panel3/panel4 moments'],
'audit.out':['PASS 673 reduced distance entries','PASS 1398 fresh exact distance pullback identities','PASS all 22 panel output coefficients']}
for name,need in checks.items():
 text=(O/name).read_text()
 for marker in need:assert marker in text,(name,marker)
 assert 'Traceback' not in text
files=list((R/'claude8').glob('panel1_*'))
assert not files,'Panel1 has appeared: perform conditional cross-check before finalizing'
for file in ['COLLAB.md','NOTES_codex.md','THEORY.md']:
 assert 'TASK_codex_15' in (R/file).read_text()
required=['rational_connection.json','distance_connection.json','outputs.json','row_clearing.json','report.md']
for file in required:assert (O/file).is_file()
inputs=['codex/phase7/connections.json','codex/phase3/panel_connection.py','codex/phase5/goal2_combined_exact.py','TASK_codex_15.md']
manifest={'inputs':{},'outputs':{},'status':{'rational_connection':'proved and verified','minimal_weight':'open; delta^4 sufficient','exact_weighted_polynomial_angular_jet':'open','claude_panel1_crosscheck':'pending; no panel1_* artifacts'}}
for path in inputs:manifest['inputs'][path]=hashlib.sha256((R/path).read_bytes()).hexdigest()
for p in sorted(O.iterdir()):
 if p.is_file() and p.name not in ['manifest.json','final_audit.out']:manifest['outputs'][p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
(O/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('PASS final saved-run checks, required exports/log entries and input/output hashes')
print('OPEN minimal weight and exact weighted polynomial angular-jet inclusion; PENDING Claude panel1 export')
