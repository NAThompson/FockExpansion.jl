"""Link the exact recurrence certificate to the displayed coefficient tables."""
from pathlib import Path
import subprocess,sys,pickle,json,hashlib
import sympy as s
HERE=Path(__file__).resolve().parent
from panel_operator import b,y,t
Z=s.Symbol('Z');table=json.loads((HERE/'inputs/final_table.json').read_text())
panels=pickle.load(open(HERE/'panels.pkl','rb'))
for name,kind,power,factor in [('angZ','angular',1,s.sqrt(2)),('angZ2','angular',2,s.sqrt(2)),('rotZ2','rotated+rr',2,1)]:
 q=sum(s.sympify(table[kind][key],locals={'b':b,'y':y,'Z':Z})*t**k for key,k in [('N0',0),('N2',2)])
 assert s.cancel(s.expand(q).coeff(Z,power)/factor-panels[name][0])==0
 print(name,'weight agrees exactly')
result=subprocess.run([sys.executable,str(HERE/'inputs/round8_classical/exact_check.py')],capture_output=True,text=True,check=True)
(HERE/'table_link.log').write_text(result.stdout+result.stderr)
assert 'ALL EXACT' in result.stdout and 'MISMATCH' not in result.stdout
print('Classical coefficient table and remainder agree exactly')
report=json.loads((HERE/'certificate.json').read_text())
report['table_link']='passed'
files=sorted(HERE.glob('*.py'))+sorted(p for p in (HERE/'inputs').rglob('*') if p.suffix in ('.py','.json','.pkl'))
report['sha256']={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(HERE/'certificate.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS: all input links checked')
