from pathlib import Path
import hashlib,json,sys
R=Path(__file__).resolve().parents[2]; O=Path(__file__).resolve().parent
manifest=O/'x31_initial.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def allowed(p):return p.parent==O and p.name.startswith('x31_')
if len(sys.argv)==1:
 files={str(p.relative_to(R)):sha(p) for p in R.rglob('*') if p.is_file() and not allowed(p) and '__pycache__' not in p.parts}
 logs={n:(R/n).stat().st_size for n in ['COLLAB.md','THEORY.md']}
 manifest.write_text(json.dumps({'files':files,'log_sizes':logs},indent=2)+'\n')
 print('initial protected snapshot',len(files))
else:
 d=json.loads(manifest.read_text());bad=[]
 for name,h in d['files'].items():
  p=R/name
  if name in d['log_sizes']:
   hh=hashlib.sha256(p.read_bytes()[:d['log_sizes'][name]]).hexdigest()
  else:hh=sha(p) if p.exists() else 'MISSING'
  if hh!=h:bad.append(name)
 new=[str(p.relative_to(R)) for p in R.rglob('*') if p.is_file() and not allowed(p) and '__pycache__' not in p.parts and str(p.relative_to(R)) not in d['files']]
 harness=[n for n in bad if n=='codex_31.log' and (R/n).read_text().startswith('Reading additional input from stdin...\nOpenAI Codex')]
 bad=[n for n in bad if n not in harness]
 print(json.dumps({'protected_changes':bad,'new_non_x31_files':new,'automatic_harness_transcript_changes':harness},indent=2));assert not bad and not new
