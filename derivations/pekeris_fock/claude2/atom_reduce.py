"""Reduce the weight-1 rows modulo the letter identities (s_h^2-r12^2) L4 = L5 L3 and (r12^2-d^2) L1 = L5 L2,
then find the log-part direction."""
import sympy as sp
from atom_closure import allrows, W1, r1, r2, r12, R, P

idx = {n: i for i, n in enumerate(W1)}


def reduce(row):
    row = [sp.sympify(x) for x in row]
    c3 = row[idx['l3']]; row[idx['l3']] = 0  # l3 = lsp + lsm + l4 - l5
    for n, s in (('lsp', 1), ('lsm', 1), ('l4', 1), ('l5', -1)):
        row[idx[n]] += s * c3
    c2 = row[idx['l2']]; row[idx['l2']] = 0  # l2 = lvp + lum + l1 - l5
    for n, s in (('lvp', 1), ('lum', 1), ('l1', 1), ('l5', -1)):
        row[idx[n]] += s * c2
    return row


logcols = [n for n in W1 if n not in ('al', 'be', 'l2', 'l3')]
for pt in [(sp.Rational(7, 10), sp.Rational(3, 10), sp.Rational(1, 2)), (sp.Rational(5, 13), sp.Rational(9, 17), sp.Rational(3, 11))]:
    sub = {r1: pt[0], r2: pt[1], r12: pt[2]}; sub[R] = sp.sqrt(pt[0]**2 + pt[1]**2); sub[P] = sp.sqrt(2 * (pt[0]**2 + pt[1]**2) - pt[2]**2)
    M = sp.Matrix([[sp.N(sp.sympify(reduce(r)[idx[n]]).subs(sub), 30) for n in W1] for r in allrows])
    full = M.rank(iszerofunc=lambda x: abs(x) < 1e-20)
    Mlog = sp.Matrix([[sp.N(sp.sympify(reduce(r)[idx[n]]).subs(sub), 30) for n in logcols] for r in allrows])
    rl = Mlog.rank(iszerofunc=lambda x: abs(x) < 1e-20)
    # direction of the log part: first nonzero row normalized
    for i in range(Mlog.rows):
        rowv = Mlog.row(i)
        if max(abs(x) for x in rowv) > 1e-15:
            j0 = max(range(len(logcols)), key=lambda j: abs(rowv[j]))
            print('rank after reduction', full, ' log-part rank', rl, ' direction', {logcols[j]: sp.N(rowv[j] / rowv[j0], 8) for j in range(len(logcols)) if abs(rowv[j]) > 1e-15})
            break
    print('point', pt, ': rank after reduction', full, ', log-part rank', rl)
# symbolic: reduced grad Omega rows (level 0) -> coefficients of al, be only?
for k, r in enumerate(allrows[:3]):
    rr = reduce(r)
    print(f'dOmega/d{(r1, r2, r12)[k]}:', {n: sp.simplify(rr[idx[n]].subs({R: sp.sqrt(r1**2 + r2**2), P: sp.sqrt(2 * (r1**2 + r2**2) - r12**2)})) for n in W1 if sp.sympify(rr[idx[n]]) != 0})
