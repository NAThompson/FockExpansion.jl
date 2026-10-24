"""Direct-operator check of one panel's Julia-generated recurrence (panel{j}_rec.txt) with ACTUAL panel atoms (Codex's 23 reduced
rational-family atoms, codex/phase8/report.md + x16_summary.json): random Laguerre arrays in all sectors x log levels 0,1;
LHS = sigma^2 [-(2v+w)(2u+w)(u+v)] e^{sigma/2}(H-E)Psi via mpmath differentiation; RHS from stencils / row factors. eps=3/2, E=-eps^2."""
import sys, json, random
import mpmath as m
import sympy as sp
j = sys.argv[1]
m.mp.dps = int(sys.argv[2]) if len(sys.argv) > 2 else 40
EPS = m.mpf(3) / 2
summ = json.load(open('../codex/phase8/x16_summary.json'))[j]
xs, ys, cs = sp.symbols('r1 r2 r12')
Hf = sp.lambdify((xs, ys, cs), sp.sympify(summ['H']), 'mpmath'); kf = sp.lambdify((xs, ys, cs), sp.sympify(summ['k']), 'mpmath')
KER = {'J': lambda s, H, T: 1, 'K': lambda s, H, T: s, 'Q0': lambda s, H, T: 1 / T, 'Q1': lambda s, H, T: s / T, 'Q2': lambda s, H, T: s**2 / T,
       'Q3': lambda s, H, T: s**3 / T, 'U0': lambda s, H, T: 1 / (s * s + 1), 'U1': lambda s, H, T: s / (s * s + 1),
       'V0': lambda s, H, T: 1 / (s * s + H), 'V1': lambda s, H, T: s / (s * s + H), 'FP': lambda s, H, T: 1 / s}
_cache = {}
def atoms(x, y, c):
    key = (x, y, c)
    if key in _cache: return _cache[key]
    H, k = Hf(x, y, c), kf(x, y, c); sH = m.sqrt(H)
    T = lambda s: s**4 + k * s * s + H
    f = lambda s: sH / m.sqrt(T(s))
    jet = {'A0': lambda s: f(s), 'A1': lambda s: f(s) * m.atan(s), 'A2': lambda s: f(s) * m.atan(s)**2 / 2,
           'B1': lambda s: f(s) * sH * m.atan(s / sH), 'B2': lambda s: f(s) * (sH * m.atan(s / sH))**2 / 2}
    out = {}
    for name in summ['reachable']:
        if name == 'F0': out[name] = m.mpf(1); continue
        kern, jt = name.split('_')
        if kern == 'E': out[name] = jet[jt](m.mpf(1)); continue
        out[name] = m.quad(lambda s: KER[kern](s, H, T(s)) * jet[jt](s), [0, 1])
    _cache[key] = out
    return out
ATOMS = sorted(summ['reachable'])
rowf = json.load(open(f'panel{j}_rowfactors.json'))
u_, v_, w_, eps_ = sp.symbols('u v w epsilon')
DEN = {a: sp.lambdify((u_, v_, w_), sp.sympify(rowf[a]).subs(eps_, sp.Rational(3, 2)), 'mpmath') for a in ATOMS}
lag = lambda nn, t: m.laguerre(nn, 0, t)
def psi(x, y, z, A):
    U, V, W = EPS * (y + z - x), EPS * (x + z - y), 2 * EPS * (x + y - z); sig = U + V + W
    at = atoms(x, y, z); tot = 0
    for (k, b), arr in A.items():
        F = sum(cf * lag(i, U) * lag(jj, V) * lag(q, W) for (i, jj, q), cf in arr.items())
        tot += m.log(sig)**k * at[b] * F
    return m.exp(-sig / 2) * tot
def HmE(f, x, y, z):
    d = lambda *o: m.diff(f, (x, y, z), o)
    lap = (d(2, 0, 0) + 2 / x * d(1, 0, 0) + d(0, 2, 0) + 2 / y * d(0, 1, 0) + 2 * d(0, 0, 2) + 4 / z * d(0, 0, 1)
           + (x*x - y*y + z*z) / (x * z) * d(1, 0, 1) + (y*y - x*x + z*z) / (y * z) * d(0, 1, 1))
    return -lap / 2 + (-2 / x - 2 / y + 1 / z + EPS**2) * f(x, y, z)
random.seed(int(j))
A = {(k, b): {(i, jj, q): m.mpf(random.randint(-5, 5)) / 7 for i in range(2) for jj in range(2) for q in range(2) if random.random() < 0.25}
     for k in (0, 1) for b in ATOMS}
import subprocess
from fractions import Fraction
PTS = [('0.62', '0.41', '0.53'), ('0.45', '0.7', '0.6')]
with open(f'arrays_{j}.txt', 'w') as f:
    for (k, b), arr in A.items():
        for (i, jj, q), cf in arr.items():
            fr = Fraction(int(cf * 7), 7)
            f.write(f"{k}|{b}|{i} {jj} {q}|{fr.numerator}|{fr.denominator}\n")
with open(f'points_{j}.txt', 'w') as f:
    for p in PTS: f.write(' '.join(p) + '\n')
subprocess.run(['julia', '-t', '2', 'rhs_eval.jl', f'panel{j}_terms.txt', f'arrays_{j}.txt', f'points_{j}.txt', f'G_{j}.txt'], check=True)
Gv = {}
for line in open(f'G_{j}.txt'):
    p, k, t, g = line.rstrip('\n').split('|'); Gv[(int(p), int(k), t)] = m.mpf(g)
for pi_, pt in enumerate(PTS, start=1):
    x, y, z = map(m.mpf, pt)
    U, V, W = EPS * (y + z - x), EPS * (x + z - y), 2 * EPS * (x + y - z); sig = U + V + W
    lhs = sig**2 * (-(2 * V + W) * (2 * U + W) * (U + V)) * m.exp(sig / 2) * HmE(lambda a, b, c_: psi(a, b, c_, A), x, y, z)
    at = atoms(x, y, z)
    r = sum(m.log(sig)**k * at[t] * g / DEN[t](U, V, W) for (p, k, t), g in Gv.items() if p == pi_)
    print('panel', j, pt, 'LHS', m.nstr(lhs, 15), 'RHS', m.nstr(r, 15), 'rel', m.nstr(abs(lhs - r) / abs(lhs), 3), flush=True)
