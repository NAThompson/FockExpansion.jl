"""Fresh actual-state comparison of saved independent connections, serial."""
import ast,json,math
from pathlib import Path
import mpmath as m
import sympy as S
import rational_panel as R
O=Path(__file__).resolve().parent;root=O.parents[1];m.mp.dps=45
# Read only function definitions, avoiding either source script's top-level jobs.
def loadfunc(path,name,ns):
 tree=ast.parse(path.read_text());node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name);exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),'exec'),ns);return ns[name]
states=loadfunc(root/'claude8/panel1_check.py','states',{'m':m})
actual=loadfunc(O/'verify.py','actual',{'m':m,'R':R})
D=json.loads((O/'rational_connection.json').read_text());F=json.loads((O/'outputs.json').read_text());CD=json.loads((root/'claude8/panel1_connection.json').read_text())
H,k=S.symbols('H k');mf={ax:S.lambdify((H,k),[S.sympify(v) for i,j,v,_ in D['matrices'][ax]['entries']],'mpmath',cse=True) for ax in ['H','k']}
ff=[S.lambdify((H,k),[S.sympify(v) for i,v in F[a]],'mpmath') for a in ['F0','F2']]
def monos(d):return [(i,t-i) for t in range(d+1) for i in range(t+1)]
def cv(r,h,Y):
 def poly(a,d):return sum(m.mpf(v)*h**i*Y**j for v,(i,j) in zip(a,monos(d)))
 return poly(r['num'],r['dn'])/poly(r['den'],r['dd'])
def moment(h,Y,X):
 # Undo symmetric quartet: sum z^j Ji = S_j. The Q identity only needs S1,S3.
 b=(h*h-1)/(h*h+1);c4=1-b;c2=4*Y-2;c0=1+b
 def q(fam):
  q=[X[fam,11,n] for n in range(3)]
  # quartic exponent -1/2; sum z_i=0, sum z_i^2=-2*c2/c4
  for n in range(3):q[n]+=( -2*c2/c4*X[fam,0,n]+X[fam,5,n])/2
  zs=[1j,-1j,1j*h,-1j*h]
  for p,z in enumerate(zs):
   a=(1/(2j) if p%2==0 else -1/(2j)) if ((fam=='A' and p<2) or (fam=='D' and p>=2)) else 0
   for n in (1,2):q[n]-=a*(z*X[fam,1,n-1]+z*z*X[fam,0,n-1]+z**3*X[fam,6+p,n-1])
  return q
 qa=q('A');qd=q('D')
 n0=m.pi/2*X['A',0,0]-2*X['A',0,1]+4/m.pi*(X['A',0,2]-X['D',0,2])
 hd=(X['D',8,1]+X['D',9,1])/2
 n2=-m.pi/2*qa[0]+2*qa[1]-4/m.pi*(qa[2]-qd[2])+4*h/m.pi*(X['D',1,1]+X['D',10,1]+4*(Y-1)/(1-b)*hd)-4*(c0+c2+c4)*X['D',11,2]/(m.pi*c4)
 return [n0,n2]
worst=m.mpf(0);checks=0
for hs,ys in [('1.31','.37'),('.79','.44')]:
 h=m.mpf(hs);Y=m.mpf(ys);hv=h*h;kv=(h*h+1)*(2*Y-1);scale=m.sqrt(2*h*h/(1+h*h));X=actual(hv,kv);CX=states(h,Y)
 def ours(hh,yy,XX):
  HH=hh*hh;kk=(hh*hh+1)*(2*yy-1);sc=m.sqrt(2*HH/(1+HH))
  return [sum(v*XX[i] for (i,_),v in zip(F[a],ff[t](HH,kk)))/sc for t,a in enumerate(['F0','F2'])]
 for got,ref in zip(ours(h,Y,X),moment(h,Y,CX)):
  err=abs(got-ref);worst=max(worst,err);assert err<m.mpf('1e-35'),err;checks+=1
 for ax in (0,1):
  dh=2*h if ax==0 else 0;dk=2*h*(2*Y-1) if ax==0 else 2*(h*h+1)
  DX=[m.mpf(0)]*55
  for a,vel in [('H',dh),('k',dk)]:
   for (i,j,_,_),v in zip(D['matrices'][a]['entries'],mf[a](hv,kv)):DX[i]+=vel*v*X[j]
  DC={key:0 for key in CX}
  for key,r in CD.items():
   aa,part,i,j=key.split('|')
   if int(aa[0])!=ax:continue
   fam=aa[1];i=int(i);j=int(j);DC[fam,i//3,i%3]+=(1 if part=='0' else 1j)*cv(r,h,Y)*CX[fam,j//3,j%3]
  om=ours(h,Y,DX);cm=moment(h,Y,DC)
  for t in (0,1):
   od=m.diff(lambda z:ours(z,Y,X)[t],h) if ax==0 else m.diff(lambda z:ours(h,z,X)[t],Y)
   cd=m.diff(lambda z:moment(z,Y,CX)[t],h) if ax==0 else m.diff(lambda z:moment(h,z,CX)[t],Y)
   err=abs(om[t]+od-cm[t]-cd);worst=max(worst,err);assert err<m.mpf('1e-32'),(h,Y,ax,t,err);checks+=1
 print('PASS fresh moments and h,Y2 derivatives',hs,ys,flush=True)
print('PASS',checks,'comparisons; worst absolute error',m.nstr(worst,8),flush=True)
(O/'x16_compare.json').write_text(json.dumps(dict(checks=checks,worst=str(worst),points=[['1.31','.37'],['.79','.44']]),indent=2)+'\n')
