import sys,json,time,itertools
from pathlib import Path
sys.dont_write_bytecode=True
from goal2x_probe import *

def save(data):
 (OUT/'goal2x_independent.json').write_text(json.dumps(data,indent=2)+'\n')

def main():
 data=json.loads((OUT/'goal2x_independent.json').read_text()) if (OUT/'goal2x_independent.json').exists() else dict(prime=PRIME,pi_default=103,points=[],checks=[],constant_results=[])
 if data.get('complete'):
  print('Already complete: saved independent two-point tests');return
 for point in [(10,17),(13,22),(25,34),(46,55)]:
  if any(tuple(v["point"])==point for v in data["points"]):continue
  try:r=relations(*point)
  except (ValueError,ZeroDivisionError):continue
  data['points'].append(dict(point=point,relations=r));save(data)
  for signs in itertools.product((1,-1),repeat=3):
   if signs==(1,1,1):continue
   rr=relations(*point,signs=signs);assert rr==r
   data['checks'].append(dict(point=point,radical_signs=signs,equal=True));save(data)
   print('PASS radical',point,signs,flush=True)
  for pi in (107,109,113,127):
   rr=relations(*point,pi=pi)
   changes=[(a,b) for a in range(10) for b in range(45) if rr[a][b]!=r[a][b]]
   data['constant_results'].append(dict(point=point,pi=pi,changed_entries=changes,relations=rr));save(data)
   print('PI',point,pi,'changes',len(changes),flush=True)
  for si,sr in [(-1,1),(1,-1),(-1,-1)]:
   rr=relations(*point,si=si,sr=sr)
   changes=[(a,b) for a in range(10) for b in range(45) if rr[a][b]!=r[a][b]]
   data['constant_results'].append(dict(point=point,si=si,sr=sr,changed_entries=changes,relations=rr));save(data)
   print('EMBED',point,si,sr,'changes',len(changes),flush=True)
  if len(data['points'])==2:break
 assert len(data['points'])==2
 data['complete']=True;save(data)
 print('PASS independent two-point branch/constant checks',flush=True)
if __name__=='__main__':main()
