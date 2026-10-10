"""Common loader: exact classical coefficients (round2) as sympy expressions in alpha, psi."""
import os, sys, pickle
import sympy as sp
H = os.path.dirname(os.path.abspath(__file__)); R2 = os.path.join(H, '..', 'round2')
al, ps = sp.symbols('alpha psi', real=True); Z = sp.Symbol('Z', positive=True)
_C = pickle.load(open(os.path.join(R2, 'classical_formula_coeffs.pkl'), 'rb'))
def coef(key):
    d = _C[key]; return sp.expand(sum(c*(Z if str(z) == 'Z' else Z**2 if str(z) == 'Z**2' else 1) for z, c in d.items()))
def cpm(s, k): return coef('4Lam((al%sps)/4+%dpi/8)' % ('+1' if s > 0 else '-1', k))
KEYS = list(_C.keys())

# Fock variables: r1 = cos(alpha/2), r2 = sin(alpha/2), xi = sqrt2 sin(psi/2) (= r12/R), eta = sqrt2 cos(psi/2)
r1, r2, xi, eta = sp.symbols('r1 r2 xi eta')
s, dl = sp.symbols('sigma delta')   # sigma = r1 + r2, delta = r1 - r2
_A2, _P2 = sp.symbols('A2 P2')
def to_fock(e):
    e = sp.expand(sp.expand_trig(sp.sympify(e).subs({al: 2*_A2, ps: 2*_P2})))
    e = sp.expand(sp.expand_trig(e))
    e = e.subs({sp.cos(_A2): r1, sp.sin(_A2): r2, sp.cos(_P2): eta/sp.sqrt(2), sp.sin(_P2): xi/sp.sqrt(2)})
    return reduce_fock(sp.expand(e))
_G = None
def reduce_fock(e, gens=(r1, r2, xi, eta)):
    """normal form modulo r1^2 + r2^2 = 1, xi^2 + eta^2 = 2 (eliminate r1^2 -> 1 - r2^2? no: eliminate eta^2 and r1^2)"""
    e = sp.expand(e)
    p = sp.Poly(e, r1, r2, xi, eta)
    out = 0
    for (i, j, k, l), c in p.terms():
        t = c*r2**j*xi**k
        t *= (1 - r2**2)**(i//2)*r1**(i % 2)*(2 - xi**2)**(l//2)*eta**(l % 2)
        out += t
    return sp.expand(out)
def nterms(e):
    e = sp.expand(e)
    return len(sp.Add.make_args(e)) if e != 0 else 0

# centred angles a = alpha - pi/2, p = psi - pi (fixed points of r1<->r2 and psi -> 2pi - psi)
a, p = sp.symbols('a p', real=True)
from sympy.simplify.fu import TR8, TR10i
def centre(e): return sp.sympify(e).subs({al: a + sp.pi/2, ps: p + sp.pi}, simultaneous=True)
def sumform(e):
    """linear combination of cos/sin(n a/2 + m p/2)"""
    return sp.expand(TR8(sp.expand(sp.expand_trig(sp.expand(e)))))
def prodform(e):
    """sum of products f(n a/2) g(m p/2) (f, g in cos, sin)"""
    A2, P2 = sp.symbols('A2c P2c')
    x = sp.expand(sp.expand_trig(sp.expand(e)))
    # expand everything in powers of cos/sin of the base angles, then go back per variable
    v = list(x.free_symbols - {Z})
    s_ = sumform(e)
    out = 0
    for t in sp.Add.make_args(s_):
        out += sp.expand(sp.expand_trig(t))
    return sp.expand(out)
def ncount(e):
    s1 = sumform(e)
    return nterms(s1)

def prod_terms(e, u=None, v=None):
    """separable form sum_i c_i f_i(n u) g_i(m v): returns dict {(fa, fb): coeff} with fa, fb trig factors"""
    if u is None: u, v = a, p
    s_ = sumform(e); out = {}
    for t in sp.Add.make_args(s_):
        if t == 0: continue
        c, f = sp.S(1), None
        for fac in sp.Mul.make_args(t):
            if isinstance(fac, (sp.sin, sp.cos)): f = fac
            else: c *= fac
        if f is None:
            out[(1, 1)] = out.get((1, 1), 0) + c; continue
        arg = sp.expand(f.args[0]); au = sp.expand(arg.coeff(u)*u); av = sp.expand(arg - au)
        if isinstance(f, sp.sin): parts = [(sp.sin(au), sp.cos(av), 1), (sp.cos(au), sp.sin(av), 1)]
        else: parts = [(sp.cos(au), sp.cos(av), 1), (sp.sin(au), sp.sin(av), -1)]
        for fa, fb, sg in parts:
            if fa == 0 or fb == 0: continue
            # normalise signs: sin(-x) -> -sin(x) handled by sympy automatically
            cc = c*sg; fa_, fb_ = fa, fb
            if fa_.could_extract_minus_sign(): fa_, cc = -fa_, -cc
            if fb_.could_extract_minus_sign(): fb_, cc = -fb_, -cc
            out[(fa_, fb_)] = out.get((fa_, fb_), 0) + cc
    return {k: sp.nsimplify(sp.simplify(c)) for k, c in out.items() if sp.simplify(c) != 0}
def prod_str(e, u=None, v=None):
    d = prod_terms(e, u, v)
    return len(d), ' + '.join('(%s)*%s*%s' % (c, fa, fb) for (fa, fb), c in d.items())

def E0_fock():
    E = pickle.load(open(os.path.join(H, '..', 'round3', 'E0_polynomial.pkl'), 'rb'))
    sym = {str(x): x for x in E.free_symbols}
    rep = {sym.get('r1', r1): r1, sym.get('r2', r2): r2}
    if 'd' in sym: rep[sym['d']] = xi
    if 'z' in sym: rep[sym['z']] = eta
    if 'Z' in sym: rep[sym['Z']] = Z
    E = E.subs(rep)
    return E, sym

def to_sigma(e):
    """express an r1<->r2 symmetric polynomial (mod r1^2+r2^2=1) in sigma = r1+r2 and xi; returns None if not symmetric"""
    e = sp.expand(e)
    # substitute r1 = (s + dl)/2, r2 = (s - dl)/2 with dl^2 = 2 - s^2
    f = sp.expand(e.subs({r1: (s + dl)/2, r2: (s - dl)/2}, simultaneous=True))
    P_ = sp.Poly(f, dl); out = 0
    for (k,), c in P_.terms():
        out += c*(2 - s**2)**(k//2)*dl**(k % 2)
    return sp.expand(out)
