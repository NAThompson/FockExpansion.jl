"""Small serial exact four-panel exports; no imported script side effects."""
import json
from pathlib import Path
import sympy as S
import sys; sys.path.insert(0,'../codex/phase8'); from rational_panel import H,k,names,C
O=Path('../codex/phase8').resolve()
import sys
x,y,c=S.symbols('r1 r2 r12');Q=x*x+y*y;e=2*Q-c*c
D=json.loads((O/'rational_connection.json').read_text());F=json.loads((O/'outputs.json').read_text())
base=[{(i,j):S.sympify(v) for i,j,v,_ in D['matrices'][a]['entries']} for a in ['H','k']]
f=[{i:S.sympify(v) for i,v in F[a]} for a in ['F0','F2']]
summary={};lines=[];raw=[];outs=[];clears=[]
for panel in range(1,3):
 a,b=(x,y) if panel in (1,3) else (y,x)
 if panel<3:
  hv=a*a/(b*b);kv=(c*c-Q)/(b*b);bb=(a*a-b*b)/Q;Y2=c*c/(2*Q)
  co=[-Q**2*(-8*bb**2+8)/(576*a),-Q**2*(1-bb)*(-18*bb-32*Y2+16)/(288*a)]  # Z^1 part (Claude)
  ep=a/c
 else:
  hv=c*c/e;kv=2*(a*a-b*b)/e
  co=[4*c*e/18,-8*e*(a*a-b*b)/(9*c)];ep=c/(2*a)
 subs={H:hv,k:kv}
 output={i:C(sum(co[t]*f[t].get(i,0).subs(subs,simultaneous=True) if isinstance(f[t].get(i,0),S.Basic) else 0 for t in (0,1))) for i in range(55)}
 cache={v:C(v.subs(subs,simultaneous=True)) for M in base for v in M.values()}
 mats={}
 for axis in (x,y,c):
  dh=S.diff(hv,axis);dk=S.diff(kv,axis);M={}
  for ij in set(base[0])|set(base[1]):
   v=C(dh*cache.get(base[0].get(ij,0),0)+dk*cache.get(base[1].get(ij,0),0))
   if v:M[ij]=v
  mats[str(axis)]=M
  for (i,j),v in sorted(M.items()):
   num,den=S.fraction(v);raw.append(f'{panel}|{axis}|{names[i]}|{names[j]}|{num}|{den}')
 # Known actual endpoint dependencies; remove all three rational/constant endpoints.
 one=names.index('F0')
 dep={names.index('E_A0'):{one:ep},names.index('E_A1'):{one:S.pi*ep/4},names.index('E_A2'):{one:S.pi**2*ep/32}}
 # Two integration-by-parts identities remove every reachable odd quartic pair.
 for fam,n in [('A',0),('A',1),('B',1)]:
  suf=fam+str(n);en=names.index('E_'+suf)
  aa={en:-S.Integer(1)};bbrow={en:S.Integer(1)}
  if n==0:aa[one]=S.Integer(1)
  elif fam=='A':
   aa[names.index('U0_A0')]=S.Integer(1)
   bbrow[names.index('J_A0')]=-S.Integer(1);bbrow[names.index('U0_A0')]=S.Integer(1)
  else:
   aa[names.index('V0_A0')]=hv
   bbrow[names.index('J_A0')]=-hv;bbrow[names.index('V0_A0')]=hv**2
  for power,terms in [(1,[(kv,aa),(-2,bbrow)]),(3,[(-2*hv,aa),(kv,bbrow)])]:
   rr={}
   for fac,row in terms:
    for j,v in row.items():
     for jj,w in dep.get(j,{j:S.Integer(1)}).items():rr[jj]=rr.get(jj,0)+fac*v*w/(kv**2-4*hv)
   dep[names.index('Q'+str(power)+'_'+suf)]={j:C(v) for j,v in rr.items() if C(v)!=0}
 def cols(row):
  out={}
  for j,v in row.items():
   for jj,w in dep.get(j,{j:S.Integer(1)}).items():out[jj]=out.get(jj,0)+v*w
  return {j:C(v) for j,v in out.items() if C(v)!=0}
 output=cols(output)
 mats={ax:{(i,j):v for i in range(55) if i not in dep for j,v in cols({j:v for (ii,j),v in M.items() if ii==i}).items()} for ax,M in mats.items()}
 reach=set(output)
 while True:
  new=reach|{j for M in mats.values() for i,j in M if i in reach}
  if new==reach:break
  reach=new
 joint={i:S.Integer(1) for i in reach}
 for ax,M in mats.items():
  for (i,j),v in sorted(M.items()):
   if i not in reach:continue
   num,den=S.fraction(v);lines.append(f'{panel}|{ax}|{names[i]}|{names[j]}|{num}|{den}');joint[i]=S.lcm(joint[i],den)
 for i,v in sorted(output.items()):outs.append(f'{panel}|output|{names[i]}|{v}')
 for i,v in sorted(joint.items()):clears.append(f'{panel}|clearing|{names[i]}|{S.factor(v)}')
 summary[panel]=dict(H=str(hv),k=str(kv),reachable=[names[i] for i in sorted(reach)],dropped=[n for i,n in enumerate(names) if i not in reach],nnz={ax:sum(i in reach for i,j in M) for ax,M in mats.items()},output_terms=len(output),max_clearing_degree=max(S.Poly(v,x,y,c).total_degree() for v in joint.values()))
 print('panel',panel,'reachable',len(reach),'nnz',summary[panel]['nnz'],'clear degree',summary[panel]['max_clearing_degree'],flush=True)
Path('z1_panel_output.txt').write_text('\n'.join(outs)+'\n')
