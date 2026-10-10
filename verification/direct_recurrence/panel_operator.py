"""Exact differential action on the two universal integral panels."""
import sympy as s, pickle
from pathlib import Path
b,y,t=s.symbols('b y t'); P=1+b+(2-4*y*y)*t*t+(1-b)*t**4
N=(1+t*t)*(1+b+(1-b)*t*t)
# A=1+2/pi asin(2yt/sqrt(N)). A_i = h_i/(pi sqrt(P)).
h={b:-2*y*t*(1-t*t)/(1+b+(1-b)*t*t),y:4*t,t:4*y*(1+b-(1-b)*t**4)/N}
metric={(b,b):-4*(1-b*b),(y,y):-(1-y*y),(b,y):-2*b*(1-2*y*y)/y}
drift={b:12*b,y:-(2-5*y*y)/y}
def action(Q):
    # operator on A Q/sqrt(P) = A R/sqrt(P) + W/pi
    D=lambda q,x:s.diff(q,x)-q*s.diff(P,x)/(2*P)
    R=-21*Q; W=0
    for (i,j),g in metric.items():
        R+=g*D(D(Q,i),j)
        W+=g*(h[i]*D(Q,j)+h[j]*D(Q,i)+Q*(s.diff(h[i],j)-h[i]*s.diff(P,j)/(2*P)))/P
    for i,g in drift.items():
        R+=g*D(Q,i); W+=g*h[i]*Q/P
    return s.factor(R),s.factor(W)
if __name__=='__main__':
    out={}
    for name,Q in [('angZ',(8*b*b-8)/576+(b-1)*(-18*b-32*y*y+16)*t*t/288),('angZ2',-(9*b*b-18)/576+(b-1)*(90*y*y-45)*t*t/288),('rotZ2',(1-b*b)/18+2*(b-1)*(2*y*y-1)*t*t/9)]:
        R,W=action(Q); print(name,'R=',R,flush=True)
        # seek R/sqrtP = d_t(S/P^(3/2)).
        cs=s.symbols('c0:8'); S=sum(cs[i]*t**i for i in range(8))
        num=s.cancel(R*P**2-(s.diff(S,t)*P-s.Rational(3,2)*S*s.diff(P,t)))
        sol=s.solve(s.Poly(num,t).all_coeffs(),cs,dict=True)
        print('primitive solutions',len(sol),flush=True)
        if sol:
            S=s.factor(S.subs(sol[0])); print('S=',S,flush=True)
            assert s.cancel(R*P**2-s.diff(S,t)*P+s.Rational(3,2)*S*s.diff(P,t))==0
            out[name]=(Q,R,W,S)
    assert len(out)==3
    pickle.dump(out,open(Path(__file__).with_name('panels.pkl'),'wb'))
