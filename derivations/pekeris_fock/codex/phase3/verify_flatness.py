"""Check the full padded connection, independent of chosen initial data."""
import mpmath as m,json
from pathlib import Path
from panel_connection import rhs
m.mp.dps=42

def matrix(b,y,axis):
 C=m.matrix(65)
 for col in range(65):
  i,j=divmod(col,5)
  xa=[[m.mpf(0)]*3 for _ in range(13)];xd=[[m.mpf(0)]*3 for _ in range(13)]
  if j==0:xa[i][0]=xd[i][0]=1
  elif j<3:xa[i][j]=1
  else:xd[i][j-2]=1
  a=rhs(b,y,axis,'A',xa);d=rhs(b,y,axis,'D',xd)
  for k in range(13):
   for l,v in enumerate([a[k][0],a[k][1],a[k][2],d[k][1],d[k][2]]):C[5*k+l,col]=v
 return C
b,y=m.mpf('.3'),m.mpf('-.6')
cb,cy=matrix(b,y,0),matrix(b,y,1)
errors=[]
for step in (m.mpf('1e-8'),m.mpf('5e-9')):
 dyb=(matrix(b+step,y,1)-matrix(b-step,y,1))/(2*step)
 dby=(matrix(b,y+step,0)-matrix(b,y-step,0))/(2*step)
 curvature=dyb-dby+cy*cb-cb*cy
 err=max(abs(x) for x in curvature)
 errors.append(err)
 print('CURVATURE step',step,'max',m.nstr(err,10),flush=True)
assert errors[-1]<m.mpf('1e-13')
assert errors[-1]<errors[0]/3
counts=[sum(abs(x)>m.mpf('1e-35') for x in c) for c in (cb,cy)]
print('65x65 padded parameter connection numerical nonzero counts at generic point:',counts)
print('These are NOT Laguerre stencil counts or a symbolic sparsity certificate.')
Path(__file__).with_name('flatness_checks.json').write_text(json.dumps(dict(errors=list(map(str,errors)),numerical_parameter_counts=counts),indent=2)+'\n')
print('PASS')
