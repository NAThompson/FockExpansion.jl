# Rigorous (ball-arithmetic) version of the phase-11 no-common-invariant test, using ore_algebra's certified transition matrices.
from ore_algebra import *
import itertools, pickle
R.<l> = QQ[]; A.<Dl> = OreAlgebra(R)
def op(fn): return sum(R(c) * Dl^m for m, c in enumerate(load(fn)))
LA, LC = op('op_e2_Q.sobj'), op('op_eC_Q.sobj')
r = QQ(1)/150
def loop(L, base, a):
    s = sign(a - base)   # approach from the base side
    p0 = a - s*r
    return L.numerical_transition_matrix([base, p0, a - s*r*I, a + s*r, a + s*r*I, p0, base], 1e-45)
aA = QQbar(sqrt(QQ(8)/9)); aC = QQbar(sqrt(QQ(24)/25))
TA1 = loop(LA, QQ(1)/5, aA); TA2 = loop(LA, QQ(1)/5, -aA)
TC1 = loop(LC, QQ(1)/3, aC); TC4 = loop(LC, QQ(1)/3, -aC)
CBF = TA1.base_ring()
Dg = diagonal_matrix(CBF, [factorial(k) for k in range(5)])
conv = lambda T: Dg * T * Dg.inverse()        # ore basis (f^(k)/k!) -> derivative jets
# jets: rows (d1^i d2^j Pi) in terms of d1^k Pi
J = {tuple(map(int, s.split()[:2])): [QQ(x) for x in s.split()[2:]] for s in open('jets_Q.txt')}
KC = matrix(QQ, [J[(0, j)] for j in range(5)])
TA1d, TA2d = conv(TA1), conv(TA2)
TC1d = KC.change_ring(CBF).inverse() * conv(TC1) * KC.change_ring(CBF)
TC4d = KC.change_ring(CBF).inverse() * conv(TC4) * KC.change_ring(CBF)
diff = max(abs(x) for x in (TA1d - TC1d).list())
print('Q1 from line A vs line C: max |difference| (ball):', diff)
# compare with the 50-digit Python result
I5 = identity_matrix(CBF, 5)
St = block_matrix([[TA1d - I5], [TA2d - I5], [TC4d - I5]])
best = None
for rows in itertools.combinations(range(15), 5):
    d = St.matrix_from_rows(rows).det()
    if not d.contains_zero():
        best = (rows, d); break
print('certified nonzero 5x5 minor of stacked (T-1):', best)
