# Local exponents and factorisation of the reconstructed telescoper.
import sys
from ore_algebra import *
tag = sys.argv[1]
cs = load('op_%s_Q.sobj' % tag)
r = len(cs) - 1
lead = cs[-1]
print(tag, 'order', r, 'degree', max(c.degree() for c in cs), 'leading:', lead.factor())
for fac, mult in lead.factor():
    a = fac.roots(QQbar)[0][0]
    K = QQbar
    RK.<z> = K[]
    ops = [RK(c.change_ring(K)(a + z)) for c in cs]
    vals = [(c.valuation(z) - m if c != 0 else oo) for m, c in enumerate(ops)]
    mn = min(vals)
    s = polygen(K, 's')
    ind = sum(ops[m][m + mn] * prod(s - j for j in range(m)) for m in range(r + 1) if vals[m] == mn)
    print('  root of', fac, ': indicial roots', [ (CC(x), e) for x, e in ind.roots()])
# infinity: lam = 1/w
w = polygen(QQ, 'w')
A = OreAlgebra(PolynomialRing(QQ, 'l'), 'Dl')
L = A(cs)
try:
    print('right factor attempt:', L.right_factor())
except Exception as e:
    print('right_factor error:', type(e).__name__, str(e)[:200])
try:
    print('local basis exponents at infinity:', L.annihilator_of_composition(1/A.base_ring().gen()).indicial_polynomial(A.base_ring().gen()).roots(QQbar))
except Exception as e:
    print('inf error:', type(e).__name__, str(e)[:150])
