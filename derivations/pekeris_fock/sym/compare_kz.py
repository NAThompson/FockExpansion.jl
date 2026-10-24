"""Spot-check our 33-term recurrence against Koutschan-Zeilberger's machine
reconstruction of Pekeris's Eq. (22) (oPEKERIS1, shift-operator form, lambda=1, e=eps)."""
import json
import sympy as sp

l, m, n, e, Z, lam, eps = sp.symbols('l m n e Z lambda epsilon')
rec = json.load(open('recurrence.json'))['P']
rec = {eval(k): sp.sympify(v.replace('lambda', 'lam'), locals={'epsilon': eps, 'Z': Z, 'lam': lam, 'l': l, 'm': m, 'n': n}) for k, v in rec.items()}
kz = {
    (-2, 0, 0): 4 * l * (l - 1) * (-e * n - e * m - e + Z),
    (-2, 0, 1): -4 * e * l * (n + 1) * (l - 1),
    (0, 0, -2): -n * (n - 1),
    (0, 0, 2): -(n + 2) * (n + 1),
    (-1, -1, 0): 4 * m * l * (-e * l - e * m + 2 * Z - 1),
    (0, -2, 1): -4 * e * m * (n + 1) * (m - 1),
    (2, 0, -1): -4 * e * n * (l + 2) * (l + 1),
    (-1, 0, -1): 2 * l * n * (-e * n - 1 + 2 * Z - e - 2 * e * m),
    (1, 1, 0): 4 * (m + 1) * (l + 1) * (-e * l - e * m - 1 + 2 * Z - 2 * e),
}
print('number of offsets:', len(rec))
for off, val in kz.items():
    ours = sp.expand(rec[off].subs(lam, 1) / (2 * eps)).subs(eps, e)
    print(off, 'match' if sp.expand(ours - val) == 0 else f'MISMATCH ours={sp.factor(ours)} kz={sp.factor(val)}')
