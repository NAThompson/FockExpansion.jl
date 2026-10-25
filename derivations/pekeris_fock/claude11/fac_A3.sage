from ore_algebra import *
R.<t> = QQ[]; A.<Dt> = OreAlgebra(R)
rows = [l.split() for l in open('op_A3.txt')]
c = {}
for k, j, v in rows: c[(int(k), int(j))] = QQ(round(2*RR(v)))
L = sum(sum(c[(k, j)] * t^j for j in range(13)) * Dt^k for k in range(3))
print(L)
print('leading factored:', factor(L.leading_coefficient()))
print('rational sols:', L.rational_solutions())
try:
    print('factor:', L.factor())
except Exception as e:
    print('factor failed', type(e), e)
for p in [0, 1, -1]:
    print(p, L.local_basis_monomials(p))
