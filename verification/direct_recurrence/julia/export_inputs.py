"""Export formula inputs, never differentiated residuals, to a portable expression DAG.
Maintenance tool only; the Julia verifier does not invoke Python.
Run after ../verify_all.py using its uv environment.
"""
from pathlib import Path
import sys,json
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import exact_residual as e
import sympy as s
nodes=[];seen={}
def emit(v):
 v=s.sympify(v)
 if v in seen:return seen[v]
 constants={s.pi:'P',s.log(2):'L',s.sqrt(2):'R',s.I:'J'}
 if v in constants:node=['symbol',constants[v]]
 elif v.is_Rational:node=['rational',str(v.p),str(v.q)]
 elif v.is_Symbol:node=['symbol',str(v)]
 elif v.is_Add:node=['add',*[emit(x) for x in v.args]]
 elif v.is_Mul:node=['mul',*[emit(x) for x in v.args]]
 elif v.is_Pow and v.exp.is_Integer:node=['pow',emit(v.base),int(v.exp)]
 else:raise ValueError(v)
 idx=len(nodes)+1;nodes.append(node);seen[v]=idx;return idx
roots={}
for tag,zp,ep,apow in [('Z0',0,0,0),('Z1',1,0,0),('Z2',2,0,0),('Z3',3,0,0),('E',0,1,0),('EZ',1,1,0),('a21',0,0,1),('a21Z',1,0,1)]:
 def coeff(v):return s.expand(v).coeff(e.Z,zp).coeff(e.E,ep).coeff(e.A21,apow)
 roots['formula_'+tag]=emit(coeff(e.C+e.state))
 roots['source_'+tag]=emit(coeff(e.src))
for name,values in e.paneldata.items():
 for label,v in zip(('Q','R','W','S'),values):roots[name+'_'+label]=emit(v)
 H=s.cancel(values[3]/e.P);roots[name+'_h0']=emit(s.expand(H).coeff(e.t,1));roots[name+'_h2']=emit(s.expand(H).coeff(e.t,3))
 rat=s.factor(values[2]-H*e.h[e.t]/e.P);D=1+e.b+(1-e.b)*e.t**2
 cs=s.symbols('c1:4');eq=s.Poly(s.cancel((rat-cs[0]*e.t/(1+e.t**2)-cs[1]*e.t/D-cs[2]*e.t/D**2)*(1+e.t**2)*D**2),e.t)
 sol=s.solve(eq.all_coeffs(),cs,dict=True)[0]
 for i,c in enumerate(cs):roots[name+'_c'+str(i+1)]=emit(sol[c])
# Elementary conic evaluations in each geometric chart, before any panel weights.
for i,pa in enumerate(e.panels+e.rot):
 bv,yv,ab,be,sinab,sy,s1mb,tanab2=pa
 uf=e.sc.universal(ab,be,s1mb,sy,-yv,tanab2,bv)
 for label,v in [('b',bv),('y',yv),('ab',e.sc.ang(*ab)),('be',e.sc.ang(*be)),('sinab',sinab),('sy',sy),('Kinf',uf['Kinf']),('K1',uf['K1'])]:roots[f'chart{i+1}_{label}']=emit(v)
# Keep both sides of each classical coefficient-table identity, not its residual.
links=[]
def capture(name,old,new):
 mapping={'r1':e.r1,'r2':e.r2,'xi':e.xi,'eta':e.eta,'delta':e.dl,'sigma':e.sg,'G':e.G,'log2':s.log(2)}
 def convert(v):
  v=s.sympify(v)
  return v.subs({x:mapping[str(x)] for x in v.free_symbols if str(x) in mapping},simultaneous=True)
 j=len(links)+1;links.append(name)
 roots[f'link{j}_left']=emit(convert(old));roots[f'link{j}_right']=emit(convert(new))
text=(HERE.parent/'inputs/round8_classical/exact_check.py').read_text()
a=text.index('def chk(');b=text.index('pairs =',a)
text=text[:a]+text[b:]
sys.path.insert(0,str(HERE.parent/'inputs/round8_classical'))
exec(compile(text,'coefficient_table_export','exec'),{'chk':capture})
(HERE/'formula_inputs.json').write_text(json.dumps({'description':'Undifferentiated formula, lower-order source, panel primitives, and conic evaluations. No residuals.','nodes':nodes,'roots':roots,'table_links':links},separators=(',',':'))+'\n')
print(len(nodes),'nodes exported')
