"""Pekeris perimetric operator, derived from the Hylleraas S-state Hamiltonian.

u = eps(r2+r12-r1), v = eps(r1+r12-r2), w = 2 eps(r1+r2-r12),
Psi = exp(-(u+v+w)/2) F(u,v,w).  Returns the polynomial-coefficient operator
P F = mult * exp(+(u+v+w)/2) (H - E)[exp(-(u+v+w)/2) F]
as a dict {(a,b,c): coeff} meaning coeff * d^a_u d^b_v d^c_w F.
"""
import sympy as sp

u, v, w, eps, E, Z, lam = sp.symbols('u v w epsilon E Z lambda')
R1, R2, R12 = sp.symbols('r1 r2 r12', positive=True)

r1 = (2 * v + w) / (4 * eps)
r2 = (2 * u + w) / (4 * eps)
r12 = (u + v) / (2 * eps)

# d(u,v,w)/d(r1,r2,r12)
du = {R1: (-eps, eps, 2 * eps), R2: (eps, -eps, 2 * eps), R12: (eps, eps, -2 * eps)}
X = (u, v, w)
Fs = sp.Function('F')(u, v, w)
sigma = u + v + w


def dr(expr, r):
    return sum(c * sp.diff(expr, x) for c, x in zip(du[r], X))


def hamiltonian_minus_E(psi):
    t = (dr(dr(psi, R1), R1) + 2 / r1 * dr(psi, R1)
         + dr(dr(psi, R2), R2) + 2 / r2 * dr(psi, R2)
         + 2 * dr(dr(psi, R12), R12) + 4 / r12 * dr(psi, R12)
         + (r1**2 - r2**2 + r12**2) / (r1 * r12) * dr(dr(psi, R1), R12)
         + (r2**2 - r1**2 + r12**2) / (r2 * r12) * dr(dr(psi, R2), R12))
    V = -Z / r1 - Z / r2 + lam / r12
    return -t / 2 + (V - E) * psi


# multiplier: -2 * (2v+w)(2u+w)(u+v) / (8 eps^3) * 4 eps^2 ... choose to clear all denominators
MULT = -sp.Rational(1, 1) * (2 * v + w) * (2 * u + w) * (u + v) / eps


def operator_dict(expr_of_F, F=Fs):
    """Collect expr (linear in F and its derivatives) into {(a,b,c): coeff}."""
    expr = sp.expand(expr_of_F)
    out = {}
    derivs = [d for d in expr.atoms(sp.Derivative) if d.expr == F] + [F]
    for d in derivs:
        if d == F:
            key = (0, 0, 0)
        else:
            cnt = dict(d.variable_count)
            key = (cnt.get(u, 0), cnt.get(v, 0), cnt.get(w, 0))
        c = sp.expand(expr.coeff(d))
        if c != 0:
            out[key] = sp.factor(c)
    # sanity: nothing left
    rest = sp.expand(expr - sum(c * (F if k == (0, 0, 0) else sp.Derivative(F, *[(x, n) for x, n in zip(X, k) if n])) for k, c in out.items()))
    assert sp.simplify(rest) == 0, rest
    return out


def pekeris_operator():
    g = sp.exp(-sigma / 2)
    expr = sp.simplify(MULT * hamiltonian_minus_E(g * Fs) / g)
    return operator_dict(sp.expand(sp.simplify(expr)))


if __name__ == '__main__':
    P = pekeris_operator()
    for k in sorted(P):
        print(k, sp.factor(P[k]))
