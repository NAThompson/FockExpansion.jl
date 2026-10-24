# Verify the closed form of psi41^(2) against (i) spectral 60-digit data, (ii) the independent Feynman-integral data (data1.txt).
import pickle, sympy as S
from mpmath import mp, mpf, sqrt, atan, asin, log, clsin, pi as mpi, sin, cos
mp.dps = 70
T = pickle.load(open('closed41.pkl', 'rb'))
r1, r2, r12, pi, a_h, b_h = S.symbols('r1 r2 r12 pi a_h b_h')
r1, r2, r12, pi = [S.Symbol(n, positive=True) for n in ('r1', 'r2', 'r12', 'pi')]
f = {k: S.lambdify((r1, r2, r12, pi, a_h, b_h), v, 'mpmath') for k, v in T.items()}
def atoms(x1, x2, x12):
    R2 = x1**2 + x2**2; P = sqrt(2 * R2 - x12**2); d, s = x1 - x2, x1 + x2; L5 = x12 * P + s * d
    al = 2 * atan(x2 / x1); be = asin((R2 - x12**2) / R2)
    cl = lambda z: clsin(2, z)
    DL = cl(al - be) + cl(mpi + be - al) - cl(mpi - al - be) - cl(al + be)
    lg = be * log((2 * x1 * x2)**2 / L5**2) + (be + al) * log((s**2 - x12**2) / (2 * x1 * x2)) + (be - al) * log((x12**2 - d**2) / (2 * x1 * x2))
    Om = lg - DL
    return {'1': 1, 'alpha': al, 'L': log((s + x12) / sqrt(R2)), 'U': Om + mpi * log((d + x12)**2 / L5), 'V': P * (be + mpi / 2)}
def F(a, t, ah, bh):
    x1, x2 = cos(a / 2), sin(a / 2)
    x12 = sqrt(2 * sin(mpi / 4 - a / 2)**2 + 2 * sin(a) * sin(t / 2)**2)
    A = atoms(x1, x2, x12)
    return sum(f[k](x1, x2, x12, mpi, ah, bh) * A[k] for k in f)
c = {l.split('|')[0].strip(): mpf(l.split('|')[1]) for l in open('coef5.txt')}
bh = c['1 r1^0 r2^0 r12^4'] / mpf(1.5); ah = (c['1 r1^4 r2^0 r12^0'] - 2 * bh) / 16
L = [list(map(mpf, l.split())) for l in open('spec41.txt')][2000:2200]
print('spectral data, 200 held-out pts, max |closed - data| =', mp.nstr(max(abs(F(a, t, ah, bh) - v) for a, t, v in L), 5))
# Feynman data = pure psi41^(2) + shift; compare up to a harmonic 2-parameter fit (shift constants are Float64 in repo)
D = [list(map(mpf, l.split())) for l in open('data1.txt')][:200]
import mpmath
rows = []; rhs = []
for a, t, v in D:
    x1, x2 = cos(a / 2), sin(a / 2)
    rows.append([4 * (x1**2 - x2**2)**2 - 1, sin(a)**2 * (3 * cos(t)**2 - 1) / 2]); rhs.append(v - F(a, t, 0, 0) + 0)
M = mpmath.matrix(rows); y = mpmath.matrix(rhs)
sol = mpmath.lu_solve(M.T * M, M.T * y)
res = max(abs(rhs[i] - sol[0] * rows[i][0] - sol[1] * rows[i][1]) for i in range(len(rows)))
print('Feynman data (independent method), 200 pts, after 2 harmonic constants: max resid =', mp.nstr(res, 5))
print('harmonic constants of Feynman-normalized form:', mp.nstr(sol[0], 25), mp.nstr(sol[1], 25))
rs = sorted([(abs(rhs[i] - sol[0] * rows[i][0] - sol[1] * rows[i][1]), D[i][0], D[i][1]) for i in range(len(rows))], reverse=True)
for r, a, t in rs[:6]: print(mp.nstr(r, 3), mp.nstr(a, 6), mp.nstr(t, 6))
print('median', mp.nstr(rs[len(rs)//2][0], 3))
