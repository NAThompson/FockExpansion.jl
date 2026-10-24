"""Exact distance pullback audit; no external source imports."""
import sympy as s
from pathlib import Path
r,t,c=s.symbols('r t c', nonzero=True)
Q=r*r+t*t
D=(Q-c*c)**2-4*r*r*t*t
maps=[((r*r-t*t)/Q,c*c/(2*Q)), (-(r*r-t*t)/Q,c*c/(2*Q)), ((c*c-Q)/Q,r*r/Q), ((c*c-Q)/Q,t*t/Q)]
lines=[]
for j,(b,y2) in enumerate(maps,1):
    discr=s.factor((4*y2-2)**2-4*(1-b*b))
    assert s.cancel(discr-4*D/Q**2)==0
    h2=s.factor((1+b)/(1-b))
    # Root collisions with the four angle poles, via T(i), T(i*h).
    Ti=s.factor((1-b)-(4*y2-2)+1+b)
    Tih=s.factor((1-b)*h2*h2-(4*y2-2)*h2+1+b)
    lines.append(f'panel {j}: Delta^2={discr}; h^2={h2}; T(i)={Ti}; T(ih)={Tih}')
assert s.factor(D-((r+t)**2-c*c)*((r-t)**2-c*c))==0
lines += ['Shared discriminant numerator D has degree 4.', 'R^2 and P^2 have degree 2; D degree 4.', 'PASS: 4 discriminant pullbacks, Heron factorization, endpoint/angle collision polynomials computed exactly.']
output='\n'.join(lines)+'\n'
print(output,end='')
Path(__file__).with_name('verify_radicals6.out').write_text(output)
