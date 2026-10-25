# 2-D jets of the S_e2 period at a point b = (x1, x2): reduce d1^i d2^j F^{-1/2} (i+j <= 4) modulo exact forms to the basis
# {d1^k F^{-1/2}, k = 0..4}, over GF(p).  Output: rows "i j c0 c1 c2 c3 c4".
import sys, os
p = Integer(os.environ.get("PRIME", 2^31 - 1)); k = GF(p)
x10, x20 = QQ(sys.argv[1]), QQ(sys.argv[2]); out = sys.argv[3]; EXTRA = int(sys.argv[4]) if len(sys.argv) > 4 else 0
R.<X,q,a,b> = k[]
D = 1 + q^4; P1, P2 = q^4 - 1, 2*q^2
x1, x2 = a, b
Px = P1*x1 + P2*x2
A = 1 - x2^2; B = P2 - x2*Px; C = D^2 - Px^2
F = X^3 + 2*B*X^2 + A*C*X
r = 4
# numerators: d1^i d2^j F^{-1/2} = N_ij / F^{i+j+1/2}
N = {(0, 0): R(1)}
for s in range(1, r + 1):
    for i in range(s + 1):
        j = s - i
        if i > 0:
            Np = N[(i - 1, j)]; m = s - 1
            N[(i, j)] = Np.derivative(a)*F - (m + k(1)/2)*Np*F.derivative(a)
        else:
            Np = N[(0, j - 1)]; m = s - 1
            N[(0, j)] = Np.derivative(b)*F - (m + k(1)/2)*Np*F.derivative(b)
S.<X,q> = k[]
phi = R.hom([X, q, k(x10), k(x20)], S)
Fs = phi(F); FX, Fq = Fs.derivative(X), Fs.derivative(q)
def lift(key):   # to common denominator F^{r+1/2}
    s = key[0] + key[1]
    return phi(N[key]) * Fs^(r - s)
sgn = k(r) - k(1)/2
dUX, dUq = 3*r + EXTRA, 8*r + 2*EXTRA
mons = [(i, j) for i in range(dUX + 1) for j in range(dUq + 1)]
exact = []
for (i, j) in mons:
    u = X^i*q^j
    exact.append(u.derivative(X)*Fs - sgn*u*FX)
    exact.append(u.derivative(q)*Fs - sgn*u*Fq)
basis = [lift((t, 0)) for t in range(5)]
cols0 = basis + exact
allm = sorted(set(sum([c.monomials() for c in cols0 + [lift(key) for key in N]], [])))
idx = {mm: t for t, mm in enumerate(allm)}
def colvec(c):
    v = [0]*len(allm)
    for co, mm in zip(c.coefficients(), c.monomials()): v[idx[mm]] = co
    return v
M0 = matrix(k, [colvec(c) for c in cols0]).transpose()
print('system', M0.dimensions(), 'rank', M0.rank()); sys.stdout.flush()
with open(out, 'w') as f:
    for key in sorted(N):
        rhs = vector(k, colvec(lift(key)))
        try:
            sol = M0.solve_right(rhs)
            f.write('%d %d %s\n' % (key[0], key[1], ' '.join(str(x) for x in sol[:5])))
        except ValueError:
            f.write('%d %d FAIL\n' % key)
        f.flush()
