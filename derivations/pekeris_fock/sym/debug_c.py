import itertools, sympy as sp
from verify_recurrence import *
from recurrence import op_scaled, compose_log_ops
params = {eps: sp.Rational(3, 2), Z: 2, lam: 1, E: -sp.Rational(7, 3)}
P = op_scaled(); S2P, C1, C2 = compose_log_ops(P)
for name, Op in [('C2', C2), ('C1', C1)]:
    rec = stencil_to_recurrence(operator_to_stencil(Op))
    for idx in itertools.product(range(3), repeat=3):
        A = {idx: sp.Integer(1)}
        Fp = sp.laguerre(idx[0], u) * sp.laguerre(idx[1], v) * sp.laguerre(idx[2], w)
        expr = 0
        for (a, b, c), coef in Op.items():
            d = Fp
            for x, k in zip(X, (a, b, c)):
                if k: d = sp.diff(d, x, k)
            expr += coef.subs(params) * d
        direct = to_laguerre(expr)
        pred = apply_rec(rec, A, params)
        bad = {k: (direct.get(k, 0), pred.get(k, 0)) for k in set(direct) | set(pred) if sp.simplify(direct.get(k, 0) - pred.get(k, 0)) != 0}
        if bad: print(name, idx, bad)
print('done')
