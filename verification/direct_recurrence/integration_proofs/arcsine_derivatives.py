import sympy as s
t,b,y=s.symbols('t b y')
N=(1+t*t)*(1+b+(1-b)*t*t)
P=1+b+(2-4*y*y)*t*t+(1-b)*t**4
assert s.expand(N-4*y*y*t*t-P)==0
v=2*y*t
expected={t:2*y*(1+b-(1-b)*t**4)/N,y:2*t,b:-y*t*(1-t*t)/(1+b+(1-b)*t*t)}
for x,e in expected.items():
    assert s.cancel(s.diff(v,x)-v*s.diff(N,x)/(2*N)-e)==0
    print(x,'exact PASS')
