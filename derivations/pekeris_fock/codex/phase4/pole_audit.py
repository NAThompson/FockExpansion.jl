from pathlib import Path
from decimal import Decimal as D,getcontext
getcontext().prec=100
root=Path(__file__).resolve().parent
for p in sorted((root/'input_states').glob('*.txt')):
 terms=[]
 for line in p.read_text().splitlines()[1:]:
  n,k,a,b,g,c=line.split()
  if (k,a,b,g)==('0','0','0','V'):
   terms.append((int(n),D(c)))
 if not terms:
  print(p.name,"no bare V term; no such pole")
  continue
 assert any(c!=0 for n,c in terms)
 print(p.name, 'pole coefficients n,c:',terms)
 print('At s=1 coefficient of pi*exp(-eps*s)/r12 =',sum((c for n,c in terms),D(0)))
 print('NONZERO polynomial in s => unscreened 1/r12 pole, infinite kinetic form')
print('Analytic: x=t*y, |t|<1: P=sqrt(1+x^2-y^2)->1, beta=asin(1-2y^2/(1+x^2))->pi/2, y*V->pi.')
print('Integral |grad(1/r12)|^2 d^3(r12) ~ integral_0^delta dr12/r12^2 diverges.')
