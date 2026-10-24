"""Differential closure of the order-2 atom set.
Field K = Q(r1,r2,r12)(R,P). Weight-1 candidates t_j (alpha, beta, logs of letters); weight-2 atom
Omega = logs - DL (constant coefficients, from chi = c*Omega/(6 pi) + weight<=1 terms).
Computes grad Omega = sum_j C_kj t_j + K-terms symbolically, checks it numerically against finite differences,
and finds the dimension of the smallest K-span of weight-1 functions closed under differentiation."""
import random
import mpmath as mp
import sympy as sp

r1, r2, r12, R, P = sp.symbols('r1 r2 r12 R P', positive=True)
RV = (r1, r2, r12)
sh, d = r1 + r2, r1 - r2
Rsq = r1**2 + r2**2
Psq = 2 * Rsq - r12**2
v = (R**2 - r12**2) / R**2
g = r12 * P / R**2
L1, L2, L3, L4, L5 = sh * P + d * r12, sh * r12 - d * P, sh * P - d * r12, sh * r12 + d * P, r12 * P + sh * d
# weight-1 atoms: name -> (gradient as explicit K expressions)
dR = [r1 / R, r2 / R, 0]
dP = [2 * r1 / P, 2 * r2 / P, -r12 / P]


def tot(expr, k):
    """total derivative wrt RV[k] of an expression in r's, R, P (no atoms)"""
    return sp.diff(expr, RV[k]) + sp.diff(expr, R) * dR[k] + sp.diff(expr, P) * dP[k]


al_grad = [-2 * r2 / R**2, 2 * r1 / R**2, 0]
vgrad = [tot(v, k) for k in range(3)]
be_grad = [vg / g for vg in vgrad]
letters = {'lr1': r1, 'lr2': r2, 'lR': R, 'lsp': sh + r12, 'lsm': sh - r12, 'lvp': r12 + d, 'lum': r12 - d,
           'l5': L5, 'l1': L1, 'l2': L2, 'l3': L3, 'l4': L4}
W1 = ['al', 'be'] + list(letters)
grad = {'al': al_grad, 'be': be_grad}
for k, L in letters.items():
    grad[k] = [tot(L, j) / L for j in range(3)]
T = {n: sp.Symbol(n) for n in W1}
# Omega = logs - DL ; logs in letter form, dDL = (l1 - l2) d(al-be) + (l3 - l4) d(al+be)
ln2 = sp.log(2)
logs = (T['be'] * (2 * (ln2 + T['lr1'] + T['lr2']) - 2 * T['l5'])
        + (T['be'] + T['al']) * (T['lsm'] + T['lsp'] - ln2 - T['lr1'] - T['lr2'])
        + (T['be'] - T['al']) * (T['lum'] + T['lvp'] - ln2 - T['lr1'] - T['lr2']))


def grad_omega(k):
    g_logs = sum(sp.diff(logs, T[n]) * grad[n][k] for n in W1)
    g_dl = (T['l1'] - T['l2']) * (grad['al'][k] - grad['be'][k]) + (T['l3'] - T['l4']) * (grad['al'][k] + grad['be'][k])
    return sp.expand(g_logs - g_dl)


GO = [grad_omega(k) for k in range(3)]
# ---------- numerical check of grad Omega against finite differences ----------
mp.mp.dps = 50


def num_atoms(p):
    a, b, c = p
    Rn = mp.sqrt(a**2 + b**2); Pn = mp.sqrt(2 * Rn**2 - c**2)
    sub = {r1: a, r2: b, r12: c, R: Rn, P: Pn}
    vals = {'al': 2 * mp.atan(b / a), 'be': mp.asin((Rn**2 - c**2) / Rn**2)}
    for k, L in letters.items():
        vals[k] = mp.log(sp.N(L.subs(sub), 50))
    return sub, vals


def omega_num(p):
    sub, vals = num_atoms(p)
    cl2 = lambda t: mp.clsin(2, t)
    al, be = vals['al'], vals['be']
    DL = cl2(al - be) + cl2(mp.pi + be - al) - cl2(mp.pi - al - be) - cl2(al + be)
    lg = sp.N(logs.subs({T[n]: vals[n] for n in W1}), 50)
    return mp.mpf(str(lg)) - DL


p0 = (mp.mpf('0.7'), mp.mpf('0.4'), mp.mpf('0.5'))
for k in range(3):
    h = mp.mpf('1e-15')
    pp = list(p0); pm = list(p0); pp[k] += h; pm[k] -= h
    fd = (omega_num(pp) - omega_num(pm)) / (2 * h)
    sub, vals = num_atoms(p0)
    an = sp.N(GO[k].subs(sub).subs({T[n]: vals[n] for n in W1}), 40)
    print(f'dOmega/d{RV[k]}: analytic {sp.N(an, 20)}  finite-diff {mp.nstr(fd, 20)}')
# ---------- closure dimension ----------
# weight-1 content: coefficient rows C_k (over W1); iterate rows -> derivative rows until K-rank (at a generic point) stabilizes.
def coeff_row(expr):
    return [sp.expand(expr).coeff(T[n]) for n in W1]


rows = [coeff_row(GO[k]) for k in range(3)]


def deriv_row(row, k):
    # d/dr_k (sum c_j t_j) = sum (dc_j) t_j + K ; keep the weight-1 part
    return [tot(c, k) for c in row]


pt = {r1: sp.Rational(7, 10), r2: sp.Rational(3, 10), r12: sp.Rational(1, 2)}
pt[R] = sp.sqrt(pt[r1]**2 + pt[r2]**2); pt[P] = sp.sqrt(2 * pt[R]**2 - pt[r12]**2)
pt2 = {r1: sp.Rational(5, 13), r2: sp.Rational(9, 17), r12: sp.Rational(3, 11)}
pt2[R] = sp.sqrt(pt2[r1]**2 + pt2[r2]**2); pt2[P] = sp.sqrt(2 * pt2[R]**2 - pt2[r12]**2)


def rank(rs, pnt):
    M = sp.Matrix([[sp.N(c.subs(pnt), 40) for c in r] for r in rs])
    return M.rank(iszerofunc=lambda x: abs(x) < 1e-25)


level = rows
allrows = list(rows)
for it in range(3):
    print(f'derivative level {it}: rows {len(allrows)}, generic K-rank {rank(allrows, pt)} / {rank(allrows, pt2)} (of {len(W1)} candidates)')
    level = [deriv_row(r, k) for r in level for k in range(3)]
    allrows += level
# which candidates are actually used
used = [n for j, n in enumerate(W1) if any(r[j] != 0 for r in rows)]
print('weight-1 candidates appearing in grad Omega:', used)
