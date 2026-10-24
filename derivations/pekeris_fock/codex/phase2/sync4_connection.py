#!/usr/bin/env python3
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import sympy as S
import mpmath as mp
import order2 as O
from verify_atoms import atoms,matrix_from_json


sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'claude2'))
import order2_atoms as C
mp.mp.dps=45
G=C.gamma_table()
a,b,d=C.RV
# Transformation from Claude rows to Codex columns in physical distances.
T=S.zeros(8,6); names=C.ATOMS
T[names.index('one'),0]=1;T[names.index('la'),1]=1;T[names.index('al'),2]=1
T[names.index('beP'),3]=1/d;T[names.index('P'),3]=S.pi/(2*d)
T[names.index('Om'),4]=1;T[names.index('lb'),4]=S.pi;T[names.index('lR'),5]=-1
CC=O.connection((a,b,d));res=[]
for k,z in enumerate((a,b,d)):
 GC=S.zeros(8)
 for source in names:
  for target,f in G[source][k].items():GC[names.index(target),names.index(source)]=f
 diff=T.diff(z)+GC*T-T*CC[k]
 for i in range(8):
  for j in range(6):
   f=S.cancel(diff[i,j]);assert f==0,(k,i,j,f)
 print('PASS exact transformed connection',z,48,'entries',flush=True)
# Also verify actual function normalizations at points on both folded/order charts.
for pt in [('0.62','0.41','0.53'),('0.41','0.62','0.53')]:
 aa,bb,dd=map(mp.mpf,pt); v=C.atom_values(aa,bb,dd)
 # Inverse of a=2v+w, b=2u+w, d=2(u+v).
 u=(dd-aa+bb)/4;vv=(dd+aa-bb)/4;w=(aa+bb-dd)/2
 got=atoms(u,vv,w)
 mapped=[v['one'],v['la'],v['al'],(v['beP']+mp.pi*v['P']/2)/dd,v['Om']+mp.pi*v['lb'],-v['lR']]
 err=max(abs(x-y) for x,y in zip(got,mapped));print('Actual basis identity',pt,mp.nstr(err,8),flush=True);assert err<mp.mpf('1e-35')
print('PASS 144 exact connection entries and 12 actual basis identities')
