"""Exact [Z] panel 1/2 outputs on the existing x16 23-state modules."""
from pathlib import Path
import sympy as S,json,resource
O=Path(__file__).resolve().parent; R=O.parents[1]; I=R/'codex/phase8'
names=json.loads((I/'rational_connection.json').read_text())['atoms']; idx={n:i for i,n in enumerate(names)}
raw=json.loads((I/'outputs.json').read_text());H,k=S.symbols('H k');x,y,c=S.symbols('r1 r2 r12');Q=x*x+y*y
mom=[{i:S.sympify(v,locals={'H':H,'k':k}) for i,v in raw[n]} for n in ('F0','F2')]
summary={};lines=[]
for panel in (1,2):
 a,d=(x,y) if panel==1 else (y,x);hv=a*a/d**2;kv=(c*c-Q)/d**2
 one=idx['F0'];ep=a/c
 dep={idx['E_A0']:{one:ep},idx['E_A1']:{one:S.pi*ep/4},idx['E_A2']:{one:S.pi**2*ep/32}}
 for fam,n in [('A',0),('A',1),('B',1)]:
  en=idx['E_'+fam+str(n)];aa={en:-S.Integer(1)};bb={en:S.Integer(1)}
  if n==0:aa[one]=S.Integer(1)
  elif fam=='A':
   aa[idx['U0_A0']]=S.Integer(1);bb[idx['J_A0']]=-S.Integer(1);bb[idx['U0_A0']]=S.Integer(1)
  else:
   aa[idx['V0_A0']]=hv;bb[idx['J_A0']]=-hv;bb[idx['V0_A0']]=hv**2
  for power,terms in [(1,[(kv,aa),(-2,bb)]),(3,[(-2*hv,aa),(kv,bb)])]:
   rr={}
   for fac,row in terms:
    for j,v in row.items():
     for jj,w in dep.get(j,{j:S.Integer(1)}).items():rr[jj]=rr.get(jj,0)+fac*v*w/(kv*kv-4*hv)
   dep[idx[f'Q{power}_{fam}{n}']]={j:S.cancel(v) for j,v in rr.items() if S.cancel(v)!=0}
 co=[-a*d*d/18,d*d*(a*a-17*d*d+8*c*c)/(72*a)]
 out={}
 for t in (0,1):
  for j,v in mom[t].items():
   v=v.subs({H:hv,k:kv},simultaneous=True)*co[t]
   for jj,w in dep.get(j,{j:S.Integer(1)}).items():out[jj]=out.get(jj,0)+v*w
 out={j:S.cancel(v) for j,v in out.items() if S.cancel(v)!=0}
 mats={}
 for line in (I/'x16_connection.txt').read_text().splitlines():
  pn,ax,nn,mm,num,den=line.split('|')
  if int(pn)==panel:mats[(ax,idx[nn],idx[mm])]=S.sympify(num)/S.sympify(den)
 reach=set(out)
 while True:
  new=reach|{j for (_,i,j) in mats if i in reach}
  if new==reach:break
  reach=new
 assert len(reach)==23
 for j,v in sorted(out.items()):lines.append(f'{panel}|Z1_R3_output|{names[j]}|{v}')
 summary[panel]={'output_terms':len(out),'reachable': [names[j] for j in sorted(reach)],'connection_nnz':{ax:sum(i in reach for (aax,i,j) in mats if aax==ax) for ax in ('r1','r2','r12')},'co0':str(co[0]),'co2':str(co[1])}
 print('panel',panel,'output',len(out),'reachable',len(reach),'nnz',summary[panel]['connection_nnz'])
(O/'x31_source_output.txt').write_text('\n'.join(lines)+'\n');(O/'x31_source_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('ru_maxrss bytes',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
print('ALL SOURCE CHECKS PASSED')
