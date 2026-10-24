"""Direct-operator check of the classical-block recurrence (classical_stencils_*.json): random Laguerre arrays in all 21 sector
atoms x log levels 0,1; LHS = sigma^2 * [-(2v+w)(2u+w)(u+v)] e^{sigma/2} (H - E) Psi by mpmath differentiation of the actual atoms;
RHS = sum_{k,b'} log(sigma)^k b'(pt) G_{k,b'}(pt) / D_{b'}(pt). eps = 3/2, E = -eps^2, Z = 2, lambda = 1, pi = actual pi."""
import json, random, glob
import mpmath as m
import sympy as sp
from classical_check import atoms as atom_vals
m.mp.dps = 60
EPS = m.mpf(3) / 2
ops = json.load(open('classical_ops.json')); ATOMS = ops['atoms']
st = {}
for f in glob.glob('classical_stencils_*.json'): st.update(json.load(open(f)))
l, n_, mm, eps, PI = sp.symbols('l n m epsilon PI')
L_, M_, N_ = sp.symbols('l m n')
def comp(rec):
    return [(eval(o), sp.lambdify((L_, M_, N_), sp.sympify(cf, locals={'epsilon': eps, 'PI': PI}).subs({eps: sp.Rational(3, 2), PI: sp.pi}), 'mpmath')) for o, cf in rec.items()]
BLK = {}
for k, v in st.items():
    name, bp, b = eval(k); BLK[(name, bp, b)] = comp(v['rec'])
u_, v_, w_ = sp.symbols('u v w')
DEN = {a: sp.lambdify((u_, v_, w_), sp.sympify(ops['D'][a], locals={'epsilon': eps}).subs(eps, sp.Rational(3, 2)), 'mpmath') for a in ATOMS}
def sector_val(name, x, y, z):
    A = atom_vals(x, y, z); Rr = m.sqrt(x * x + y * y); Pp = m.sqrt(2 * Rr * Rr - z * z)
    i = int(name.lstrip('P').lstrip('g'))
    return A[i] * (Pp if name.startswith('P') else 1)
lag = lambda n, x: m.laguerre(n, 0, x)
def psi(x, y, z, A):
    U, V, W = EPS * (y + z - x), EPS * (x + z - y), 2 * EPS * (x + y - z); sig = U + V + W
    tot = 0
    for (k, b), arr in A.items():
        F = sum(cf * lag(i, U) * lag(j, V) * lag(q, W) for (i, j, q), cf in arr.items())
        tot += m.log(sig)**k * sector_val(b, x, y, z) * F
    return m.exp(-sig / 2) * tot
def HmE(f, x, y, z):
    E = -EPS**2; d = lambda *o: m.diff(f, (x, y, z), o)
    lap = (d(2, 0, 0) + 2 / x * d(1, 0, 0) + d(0, 2, 0) + 2 / y * d(0, 1, 0) + 2 * d(0, 0, 2) + 4 / z * d(0, 0, 1)
           + (x*x - y*y + z*z) / (x * z) * d(1, 0, 1) + (y*y - x*x + z*z) / (y * z) * d(0, 1, 1))
    return -lap / 2 + (-2 / x - 2 / y + 1 / z + EPS**2) * f(x, y, z)
def rhs(x, y, z, A):
    U, V, W = EPS * (y + z - x), EPS * (x + z - y), 2 * EPS * (x + y - z); sig = U + V + W
    tot = 0
    for k in (0, 1):
        for bp in ATOMS:
            G = 0
            for name, kk, w8 in (('sigma2P', k, 1), ('C1', k + 1, k + 1), ('C2', k + 2, (k + 1) * (k + 2))):
                for (nm, t, b), stl in BLK.items():
                    if nm != name or t != bp or (kk, b) not in A: continue
                    arr = A[(kk, b)]
                    for e in {(i - o[0], j - o[1], q - o[2]) for (i, j, q) in arr for o, _ in stl}:
                        if min(e) < 0: continue
                        val = sum(f(*e) * arr[(e[0] + o[0], e[1] + o[1], e[2] + o[2])] for o, f in stl if (e[0] + o[0], e[1] + o[1], e[2] + o[2]) in arr)
                        if val != 0: G += w8 * val * lag(e[0], U) * lag(e[1], V) * lag(e[2], W)
            tot += m.log(sig)**k * sector_val(bp, x, y, z) * G / DEN[bp](U, V, W)
    return tot
random.seed(5)
A = {(k, b): {(i, j, q): m.mpf(random.randint(-5, 5)) / 7 for i in range(2) for j in range(2) for q in range(2) if random.random() < 0.3}
     for k in (0, 1) for b in ATOMS}
for pt in [('0.62', '0.41', '0.53'), ('0.35', '0.8', '0.6')]:
    x, y, z = map(m.mpf, pt)
    U, V, W = EPS * (y + z - x), EPS * (x + z - y), 2 * EPS * (x + y - z); sig = U + V + W
    lhs = sig**2 * (-(2 * V + W) * (2 * U + W) * (U + V)) * m.exp(sig / 2) * HmE(lambda a, b, c_: psi(a, b, c_, A), x, y, z)
    r = rhs(x, y, z, A)
    print(pt, 'LHS', m.nstr(lhs, 18), 'RHS', m.nstr(r, 18), 'rel', m.nstr(abs(lhs - r) / abs(lhs), 3), flush=True)
