"""Symbolic version of the classical-part atomization (exact coefficients).
Algebraic quantities are rational in T = tan(al/4), U = tan(ps/4) (and sqrt2); atoms as in atoms.py with exact constants.
Output: total classical part (per Z power) as sum of coefficient * atom-monomial, coefficients simplified."""
import os as _os
_HERE_ = _os.path.dirname(_os.path.abspath(__file__))
import sys, json, pickle
sys.path.insert(0, '..')
import sympy as sp
from fractions import Fraction as Fr
T, U, Z, AL, PS, Gc = sp.symbols('T U Z AL PS G')
pi, log2, r2 = sp.pi, sp.log(2), sp.sqrt(2)
ca, sa = (1-T**2)/(1+T**2), 2*T/(1+T**2)          # cos(al/2), sin(al/2)
cp, sp_ = (1-U**2)/(1+U**2), 2*U/(1+U**2)          # cos(ps/2), sin(ps/2)
REG = {}
def lam(m1, m2, k):
    m1, m2, k = int(m1), int(m2), int(k)
    if m1 == 0 and m2 == 0: return sp.log(sp.Abs(2*sp.sin(k*pi/8)))
    if m1 % 2 == 0 and m2 % 2 == 0 and k % 2 == 0: return lam(m1//2, m2//2, k//2)+lam(m1//2, m2//2, k//2+4)
    if (m1, m2) < (0, 0): m1, m2, k = -m1, -m2, -k
    k %= 8; s = sp.Symbol('L_%d_%d_%d' % (m1, m2, k)); REG[s] = ('lam', (m1, m2, k)); return s
def cl(m1, m2, k):
    m1, m2, k = int(m1), int(m2), int(k)
    if m1 == 0 and m2 == 0: return sp.Symbol('Cl2c_%d' % (k % 16))     # constant Clausen value (kept symbolic)
    if m1 % 2 == 0 and m2 % 2 == 0 and k % 2 == 0: return 2*cl(m1//2, m2//2, k//2)+2*cl(m1//2, m2//2, k//2+8)
    sg = 1
    if (m1, m2) < (0, 0): m1, m2, k, sg = -m1, -m2, -k, -1
    k %= 16; s = sp.Symbol('C_%d_%d_%d' % (m1, m2, k)); REG[s] = ('cl', (m1, m2, k)); return sg*s
Rat = lambda x: sp.Rational(Fr(x).numerator, Fr(x).denominator)
def ang(ca_, cp_, cpi): return Rat(ca_)*AL+Rat(cp_)*PS+pi*Rat(cpi)
def form(t, shift=0): return (Fr(t[0])*4, Fr(t[1])*4, (Fr(t[2])+Fr(shift))*8)
def lam_a(t, shift=0): return lam(*form(t, shift))
def cl_a(t, shift=0): return cl(*form(t, shift))
half = lambda t: tuple(Fr(x)/2 for x in t)
def add_(t1, t2, s1=1, s2=1): return tuple(s1*Fr(x)+s2*Fr(y) for x, y in zip(t1, t2))
b_, y_, Z_ = sp.symbols('b y Z')
TAB = json.load(open(_os.path.join(_HERE_, '../final_table.json'))); CONST = json.load(open(_os.path.join(_HERE_, '../constants_simplified.json')))
def universal(ab, bet, s1mb, sy, absy, tanab2, bval):
    out = {}; abx, betx = ang(*ab), ang(*bet)
    g = [add_((0, 0, Fr(1, 4)), half(ab)), add_((0, 0, Fr(1, 4)), half(ab), 1, -1), add_((0, 0, -Fr(1, 4)), half(ab)), add_((0, 0, -Fr(1, 4)), half(ab), 1, -1)]
    bmg = [add_(bet, gk, 1, -1) for gk in g]
    ls = [lam_a(half(x))-log2 for x in bmg]; lc = [lam_a(half(x), Fr(1, 2))-log2 for x in bmg]
    Cm = [cl_a(x) for x in bmg]; Cp = [cl_a(x, 1) for x in bmg]
    s = [-1, -1, 1, 1]; u = [1, -1, -1, 1]
    v = sum(abx/4*s[k]*(ls[k]-lc[k])+betx/2*u[k]*(ls[k]-lc[k])+u[k]*(Cm[k]-Cp[k])/2 for k in range(4))
    v += pi/8*(-ls[0]-3*lc[0]+ls[1]+3*lc[1]-ls[2]+5*lc[2]+ls[3]-5*lc[3])
    out['Kinf'] = v/(pi*s1mb)
    s = [-1, 1, -1, 1]
    lcb = lam_a(bet, Fr(1, 2))-log2
    v = sum(abx/4*s[k]*(ls[k]+lc[k])+betx/2*(ls[k]+lc[k])+(Cm[k]+Cp[k])/2 for k in range(4))
    v += -2*betx*lcb-cl_a(add_(bet, bet), 1)+pi/8*(-ls[0]+3*lc[0]-ls[1]+3*lc[1]+ls[2]+5*lc[2]+ls[3]+5*lc[3])-pi*lcb+log2*betx+pi*log2/2
    out['K1'] = 2*v/(pi*2*sy)
    s = [1, -1, -1, 1]; u = [-1, -1, 1, 1]
    v = sum(abx/4*s[k]*(ls[k]+lc[k])+betx/2*u[k]*(ls[k]+lc[k])+u[k]*(Cm[k]+Cp[k])/2 for k in range(4))
    cs = lam_a(add_(half(ab), (0, 0, Fr(1, 4))))-log2/2
    if bval != 0: v += (abx/2-pi/4)*(lam_a(add_(half(ab), (0, 0, Fr(3, 4))))-log2/2)
    v += -(abx/2+pi/4)*cs-cl_a(add_(ab, (0, 0, Fr(1, 2))))/2+cl_a(add_(ab, (0, 0, -Fr(1, 2))))/2
    v += pi/8*(ls[0]-3*lc[0]+ls[1]-3*lc[1]+ls[2]+5*lc[2]+ls[3]+5*lc[3])+3*pi*log2/4
    out['Km1'] = v/(pi*2*absy)
    out['Lg2'] = log2; out['R1'] = pi/4; out['L1'] = pi**2/16
    out['Rrho'] = abx/2*tanab2; out['Lrho'] = abx**2/4; out['JL'] = pi**2/8-abx**2/2
    out['Lreg'] = (pi**2/4-abx**2)/(4*bval) if bval != 0 else None
    uang = add_((0, 0, Fr(1, 4)), half(ab), 1, -1); ux = ang(*uang)
    ti2 = ux*(lam_a(uang)-lam_a(uang, Fr(1, 2)))+(cl_a(add_(uang, uang))+cl_a(add_((0, 0, 1), add_(uang, uang), 1, -1)))/2
    out['J1'] = pi/4*(2*lam_a(add_(half(ab), (0, 0, Fr(1, 4))))-2*log2)+ti2
    return out
def panel(name, bval, yval, s1mb, s1pb, sy, absy, tanab2, ab, bet):
    Uf = universal(ab, bet, s1mb, sy, absy, tanab2, bval)
    sub = {b_: bval, y_: yval}
    tot = 0
    Ah1 = 1+2*ang(*bet)/pi
    for key, val in CONST[name].items():
        k = sp.sympify(key, locals={'sb': sp.Symbol('sb'), 'sy': sp.Symbol('sy'), 'log2': sp.Symbol('log2')}).subs({sp.Symbol('sb'): s1pb, sp.Symbol('sy'): sy, sp.Symbol('log2'): log2})
        vv = sp.sympify(val, locals={'b': b_, 'y': y_, 'Z': Z, 'Ah1': sp.Symbol('Ah1'), 'log2': sp.Symbol('log2')}).subs(sub).subs({sp.Symbol('log2'): log2, sp.Symbol('Ah1'): Ah1})
        tot += k*vv
    for k, v in TAB[name].items():
        if k in ('C', 'N0', 'N2'): continue
        if Uf[k] is None: continue
        tot += sp.sympify(v, locals={'b': b_, 'y': y_, 'Z': Z}).subs(sub)*Uf[k]
    return tot
def all_panels():
    a = ca**2-sa**2; w = cp**2-sp_**2
    tot = panel('angular', a, -sp_, r2*sa, r2*ca, cp, sp_, sa/ca, (1, 0, 0), (0, -Fr(1, 2), 0))
    tot += panel('angular', -a, -sp_, r2*ca, r2*sa, cp, sp_, ca/sa, (-1, 0, 1), (0, -Fr(1, 2), 0))
    for yv, syv, betv in ((-ca, sa, (Fr(1, 2), 0, -Fr(1, 2))), (-sa, ca, (-Fr(1, 2), 0, 0))):
        tot += panel('rotated+rr', -w, yv, r2*cp, r2*sp_, syv, -yv, cp/sp_, (0, -1, 1), betv)
        tot += Z**2/36*panel('rr0', 0, yv, 1, 1, syv, -yv, 1, (0, 0, Fr(1, 2)), betv)
    return tot
if __name__ == '__main__':
    import time; t0 = time.time()
    e = all_panels()
    pickle.dump(e, open(_os.path.join(_HERE_, 'sym_panels.pkl'), 'wb'))
    print('panels built', time.time()-t0, 'atoms', len([s for s in e.free_symbols if s in REG]))
