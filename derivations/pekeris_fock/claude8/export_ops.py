"""Export cleared operators (claude7/classical_ops.json format) as exact monomial term lists for the Julia stencil generator:
{block_key: [[deriv(3), [i,j,k,a,b] (u^i v^j w^k eps^a PI^b), num, den], ...]}."""
import json, sys
import sympy as sp
u, v, w, eps, PI = sp.symbols('u v w epsilon PI')
src, dst = sys.argv[1], sys.argv[2]
d = json.load(open(src))
out = {}
for key, op in d['ops'].items():
    terms = []
    for dk, s in op.items():
        poly = sp.Poly(sp.sympify(s, locals={'epsilon': eps, 'PI': PI, 'u': u, 'v': v, 'w': w}), u, v, w, eps, PI)
        for mono, cf in poly.terms():
            cf = sp.Rational(cf)
            terms.append([list(eval(dk)), list(mono), str(cf.p), str(cf.q)])
    out[key] = terms
json.dump(out, open(dst, 'w'))
print(len(out), 'blocks exported')
