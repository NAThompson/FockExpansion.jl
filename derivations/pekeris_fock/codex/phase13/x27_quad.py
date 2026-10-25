from pathlib import Path
import mpmath as m,json,resource
m.mp.dps=65
p=m.pi
r=lambda x,y:m.sqrt(2)/4*m.sin((x+y)/2)**2*m.cos((y-x)/4)
q=m.quad(lambda x:m.sin(x)*m.quad(lambda y:m.sin(y)*r(x,y),[p/2,p,3*p/2]),[p/2,p,3*p/2])*(-2/p)
e=16*(184-191*m.sqrt(2))/(1575*p)
assert abs(q-e)<m.mpf('1e-60')
out={'test_input':'phi3=sin(alpha), NOT the physical ladder', 'digits':65,'quadrature':str(q),'exact_value':str(e),'absolute_error':str(abs(q-e)),'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
