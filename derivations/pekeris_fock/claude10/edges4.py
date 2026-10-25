# Inertias of T5 at the four physical collinear arcs, in one basis (line-A jet at b = (1/5, 1/3)); common invariant subspace.
from monodromy import *
from mpmath import mpf, mpc, sqrt, nstr, matrix, eye, svd, binomial, inverse, norm
from fractions import Fraction
def rd(s):
    return mpf(Fraction(s).numerator) / mpf(Fraction(s).denominator)
def load_ops(fn): return [[rd(x) for x in l.split()] for l in open(fn)]
J = {tuple(map(int, l.split()[:2])): [rd(x) for x in l.split()[2:]] for l in open('jets_Q.txt')}
def Kmat(dirvec):
    # rows: (a d1 + c d2)^j Pi in terms of d1^k Pi
    a, c = dirvec; K = matrix(5, 5)
    for j in range(5):
        for i in range(j + 1):
            w = binomial(j, i) * a**i * c**(j - i)
            for k in range(5): K[j, k] += w * J[(i, j - i)][k]
    return K
def inert(ops, base, pt, sing, K):
    T = monodromy_matrix(ops, base, pt, mpf(1) / 150, sing)
    return inverse(K) * T * K
b1, b2 = mpf(1) / 5, mpf(1) / 3
res = {}
# line A: x2 = 1/3, param x1
opsA = load_ops('op_e2_A.txt'); l0 = sqrt(mpf(8) / 9)
singA = [l0, -l0, mpf(1), mpf(-1)]
res['A+ (Q1)'] = inert(opsA, b1, l0, singA, eye(5)); res['A- (Q2)'] = inert(opsA, b1, -l0, singA, eye(5))
# line C: x1 = 1/5, param x2
opsC = load_ops('op_eC.txt'); m0 = sqrt(mpf(24) / 25)
singC = [m0, -m0, mpf(1), mpf(-1), mpf('1.54950727000600'), mpf('-1.54950727000600')] + [mpc(s1 * mpf('0.965152178692833'), s2 * mpf('0.121910872334630')) for s1 in (1, -1) for s2 in (1, -1)]
KC = Kmat((0, 1))
res['C+ (Q1)'] = inert(opsC, b2, m0, singC, KC); res['C- (Q4)'] = inert(opsC, b2, -m0, singC, KC)
import os
if os.path.exists('op_eD.txt'):
    # line D: x1 = l, x2 = l + 2/15
    opsD = load_ops('op_eD.txt')
    c = mpf(2) / 15; disc = sqrt(2 - c**2)   # 2 l^2 + 2 c l + c^2 - 1 = 0
    lp, lm = (-c + disc) / 2, (-c - disc) / 2
    singD = [lp, lm] + [mpc(x) for x in __import__('json').load(open('singD.json'))]
    KD = Kmat((1, 1))
    res['D+ (Q1)'] = inert(opsD, b1, lp, singD, KD); res['D- (Q3)'] = inert(opsD, b1, lm, singD, KD)
for k, T in res.items():
    print(k, 'rank(T-1) sv:', [nstr(x, 3) for x in svd(T - eye(5), compute_uv=False)])
print('consistency |T_A+ - T_C+| =', nstr(norm(res['A+ (Q1)'] - res['C+ (Q1)']), 3))
St = matrix(20, 5)
keys = [k for k in ['A+ (Q1)', 'A- (Q2)', 'C- (Q4)', 'D- (Q3)'] if k in res]
St = matrix(5 * len(keys), 5)
for n, key in enumerate(keys):
    for i in range(5):
        for j in range(5): St[5 * n + i, j] = (res[key] - eye(5))[i, j]
print('common invariants of Q1..Q4 inertias, singular values:', [nstr(x, 3) for x in svd(St, compute_uv=False)])
