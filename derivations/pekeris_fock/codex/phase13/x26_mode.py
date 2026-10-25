# Independent Sturm flux computation, using only explicit psi10 and psi21.
import mpmath as m,json,time,resource
from pathlib import Path
m.mp.dps=45
p=m.pi
# Z coefficient of -wh/16, wh=w*8psi21-2wV*psi10.
def r(x,y):
 a=(x+y)/2;b=(y-x)/2
 w=2*m.sin(a)*m.sin(b)
 e=2*m.sqrt(2)*m.sin(a)*m.cos(b/2)
 n=-4*m.sin(b)*(m.cos(a/2)+m.sin(a/2))
 f0=m.sin(b/2)/m.sqrt(2);f1=-(m.cos(a/2)+m.sin(a/2))
 p21=-(p-2)/(3*p)*m.cos(b)
 return -(8*w*p21-2*(e*f1+n*f0))/16
out={}
for tag,h in [('cos',lambda x:m.cos(x/2)),('sin',lambda x:m.sin(x/2))]:
 # Split quadratures use entire trigonometric sources, no spectral continuation.
 integ=m.quad(lambda x:h(x)*m.quad(lambda y:h(y)*r(x,y),[p/2,p,3*p/2])/(p/2),[p/2,p,3*p/2])
 out[tag+'_positive_log_coefficient']=m.nstr(integ/h(3*p/2),43)
 out[tag+'_code_negative_log_coefficient']=m.nstr(-integ/h(3*p/2),43)
print(json.dumps(out,indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
