"""Find a constant-coefficient basis of the 5-dim K-span of weight-1 functions (if one exists):
intersect the row spaces at several generic points."""
import sympy as sp
import mpmath as mp
from atom_closure import allrows, W1, r1, r2, r12, R, P

mp.mp.dps = 40
pts = [(sp.Rational(7, 10), sp.Rational(3, 10), sp.Rational(1, 2)), (sp.Rational(5, 13), sp.Rational(9, 17), sp.Rational(3, 11)),
       (sp.Rational(2, 3), sp.Rational(1, 7), sp.Rational(5, 7)), (sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 3))]
bases = []
for a, b, c in pts:
    sub = {r1: a, r2: b, r12: c}; sub[R] = sp.sqrt(a**2 + b**2); sub[P] = sp.sqrt(2 * (a**2 + b**2) - c**2)
    M = mp.matrix([[mp.mpf(str(sp.N(e.subs(sub), 45))) for e in row] for row in allrows])
    # orthonormal basis of row space via SVD
    U, S, V = mp.svd_r(M)
    k = sum(1 for x in S if x > 1e-20)
    bases.append(mp.matrix([[V[i, j] for j in range(len(W1))] for i in range(k)]))
# intersection: vectors n with (I - Q_i) n = 0 for all i
n = len(W1)
blocks = []
for B in bases:
    Q = B.T * B
    blocks.append(mp.eye(n) - Q)
Big = mp.matrix(n * len(blocks), n)
for t, Bl in enumerate(blocks):
    for i in range(n):
        for j in range(n):
            Big[t * n + i, j] = Bl[i, j]
U, S, V = mp.svd_r(Big)
null = [i for i in range(n) if S[i] < 1e-20]
print('dimension of constant-coefficient subspace:', len(null), ' singular values:', [mp.nstr(S[i], 3) for i in range(n)])
# rational basis via reduced row echelon of null vectors
N = sp.Matrix([[sp.nsimplify(float(V[i, j]), tolerance=1e-12, rational=True) for j in range(n)] for i in null])
print('W1 order:', W1)
for row in N.rref()[0].tolist():
    if any(x != 0 for x in row):
        print({W1[j]: x for j, x in enumerate(row) if x != 0})
