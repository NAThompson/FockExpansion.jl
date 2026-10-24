"""Convert sympy stencil shards (claude7/classical_stencils_*.json) to the same monomial text format as stencil2.jl (one shard per process)."""
import json, sys
import sympy as sp
l, m, n, eps, PI = sp.symbols('l m n epsilon PI')
k = int(sys.argv[1])
d = json.load(open(f'../claude7/classical_stencils_{k}.json'))
with open(f'classical_rec_sympy_{k}.txt', 'w') as f:
    for key, v in d.items():
        for off, cf in v['rec'].items():
            poly = sp.Poly(sp.sympify(cf, locals={'epsilon': eps, 'PI': PI, 'l': l, 'm': m, 'n': n}), l, m, n, eps, PI)
            o = ' '.join(map(str, eval(off)))
            for mono, c in poly.terms():
                c = sp.Rational(c)
                f.write(f"{key}|{o}|{' '.join(map(str, mono))}|{c.p}|{c.q}\n")
