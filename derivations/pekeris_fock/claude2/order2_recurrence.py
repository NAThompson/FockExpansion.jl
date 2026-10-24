"""Recurrence for Psi = exp(-sigma/2) sum_k log(sigma)^k sum_b b F_{k,b}, b in the order-2 closed atom set
(order2_atoms.ATOMS). Fast rational-function arithmetic in the field Q(u,v,w,eps); Z=2, lambda=1, E=-eps^2.

Block (target b') <- (source b) of the log-free operator:
  O_{b'b} = delta P + 2 sum_kl M_kl Gamma^k_{b b'} d_l
            + [sum_kl M_kl (d_l Gamma^k_{bb'} + sum_m Gamma^k_{bm} Gamma^l_{mb'}) + sum_k bvec_k Gamma^k_{bb'}],
with d_k b = sum_b' Gamma^k_{bb'} b'. Log tower: (k,b') <- (k,b): sigma^2 O; (k+1,b): (k+1)[delta C1 + 2 sigma sum_l m_l Gamma^l_{bb'}];
(k+2,b): (k+1)(k+2) delta C2. Each target row (b') is multiplied by its own clearing polynomial D_{b'}."""
import json
import sys
import sympy as sp
sys.path.insert(0, '../sym')
from pekeris_operator import u, v, w, eps, E, Z, lam
from recurrence import op_scaled, split_operator, operator_to_stencil, stencil_to_recurrence
import order2_atoms as oa

FLD, fu, fv, fw, fe = sp.field([u, v, w, eps], sp.QQ)
FX = (fu, fv, fw)
SX = (u, v, w)
rsub = {oa.r1: (2 * v + w) / (4 * eps), oa.r2: (2 * u + w) / (4 * eps), oa.r12: (u + v) / (2 * eps)}
JAC = [[sp.diff(rsub[r], x) for r in oa.RV] for x in SX]
sigma = fu + fv + fw


def to_f(expr):
    return FLD.from_expr(sp.sympify(expr).subs({E: -eps**2, Z: 2, lam: 1}))


def gamma_uvw():
    G = oa.gamma_table()
    out = {a: [dict() for _ in range(3)] for a in oa.ATOMS}
    for a in oa.ATOMS:
        for xi in range(3):
            acc = {}
            for k in range(3):
                for t, cf in G[a][k].items():
                    acc[t] = acc.get(t, 0) + JAC[xi][k] * cf.subs(rsub, simultaneous=True)
            out[a][xi] = {t: to_f(sp.cancel(c)) for t, c in acc.items() if sp.cancel(c) != 0}
    return out


def fdiff(f, i):
    return f.diff(FX[i])


def build(Gm):
    P = op_scaled()
    M, bvec, c = split_operator(P)
    Mf = [[to_f(M[i, j]) for j in range(3)] for i in range(3)]
    bf = [to_f(x) for x in bvec]
    Pf = {k: to_f(cf) for k, cf in P.items()}
    msum = [sum((Mf[i][j] for i in range(3)), FLD(0)) for j in range(3)]
    trM = sum((Mf[i][j] for i in range(3) for j in range(3)), FLD(0))
    A = oa.ATOMS
    O, C1, C2 = {}, {}, {}
    for b in A:
        for bp in A:
            op = {}
            if b == bp:
                op = dict(Pf)
            z0 = FLD(0)
            for k in range(3):
                gk = Gm[b][k].get(bp, FLD(0))
                for l in range(3):
                    if gk != 0 and Mf[k][l] != 0:
                        key = tuple(1 if t == l else 0 for t in range(3))
                        op[key] = op.get(key, FLD(0)) + 2 * Mf[k][l] * gk
                    if Mf[k][l] != 0:
                        z0 += Mf[k][l] * (fdiff(gk, l) + sum((Gm[b][k].get(m, FLD(0)) * Gm[m][l].get(bp, FLD(0)) for m in A), FLD(0)))
                z0 += bf[k] * gk
            if z0 != 0:
                op[(0, 0, 0)] = op.get((0, 0, 0), FLD(0)) + z0
            op = {kk: cf for kk, cf in op.items() if cf != 0}
            if op:
                O[(bp, b)] = op
            # C1
            op1 = {}
            if b == bp:
                for j in range(3):
                    op1[tuple(1 if t == j else 0 for t in range(3))] = 2 * sigma * msum[j]
                op1[(0, 0, 0)] = sigma * sum(bf, FLD(0)) - trM
            extra = sum((2 * sigma * msum[l] * Gm[b][l].get(bp, FLD(0)) for l in range(3)), FLD(0))
            if extra != 0:
                op1[(0, 0, 0)] = op1.get((0, 0, 0), FLD(0)) + extra
            op1 = {kk: cf for kk, cf in op1.items() if cf != 0}
            if op1:
                C1[(bp, b)] = op1
            if b == bp:
                C2[(bp, b)] = {(0, 0, 0): trM}
    return O, C1, C2


def row_denominators(O, C1):
    """D_{b'} clears every block in row b' (sigma^2 O and C1); returned as field element (polynomial)."""
    D = {}
    for blk, mult in ((O, sigma**2), (C1, FLD(1))):
        for (bp, b), op in blk.items():
            for cf in op.values():
                den = (mult * cf).denom
                D[bp] = D.get(bp, den.ring.one).lcm(den)
    # Laguerre finiteness: every coefficient of d/dx_l must be divisible by x_l (x L' is finite, bare L' is not).
    # Multiply the row by x_l where needed.
    for bp in list(D):
        for _ in range(3):
            changed = False
            for blk, mult in ((O, sigma**2), (C1, FLD(1))):
                for (t, b), op in blk.items():
                    if t != bp:
                        continue
                    for kk, cf in op.items():
                        if sum(kk) != 1:
                            continue
                        l = kk.index(1)
                        q = FLD.from_expr(D[bp].as_expr()) * mult * cf
                        assert q.denom == q.denom.ring.one or q.denom.as_expr().free_symbols <= {eps}
                        num = q.numer
                        if num.rem(FX[l].numer) != 0:
                            D[bp] = D[bp] * FX[l].numer
                            changed = True
            if not changed:
                break
    return D


def to_sym(fe_):
    num = fe_.numer.as_expr()
    den = fe_.denom.as_expr()
    q = sp.cancel(num / den)
    assert sp.fraction(q)[1].free_symbols <= {eps}, sp.fraction(q)[1]
    return q


if __name__ == '__main__':
    Gm = gamma_uvw()
    O, C1, C2 = build(Gm)
    D = row_denominators(O, C1)
    result = {'atoms': oa.ATOMS, 'D': {}, 'sigma2P': {}, 'C1': {}, 'C2': {}}
    for bp in oa.ATOMS:
        Dexpr = sp.factor(D[bp].as_expr())
        result['D'][bp] = str(Dexpr)
        print(f'row {bp:>4}: clearing factor {Dexpr}  (degree {sp.Poly(D[bp].as_expr(), u, v, w).total_degree()})')
    tot = 0
    for name, blk, mult in (('sigma2P', O, sigma**2), ('C1', C1, FLD(1)), ('C2', C2, FLD(1))):
        for (bp, b), op in sorted(blk.items()):
            Dr = FLD.from_expr(D[bp].as_expr())
            sym_op = {kk: sp.expand(to_sym(Dr * mult * cf)) for kk, cf in op.items()}
            rec = stencil_to_recurrence(operator_to_stencil(sym_op))
            dN = [sum(o) for o in rec]
            tot += len(rec)
            print(f'{name:7s} {bp:>4} <- {b:<4}: {len(rec):5d} offsets, dN in [{min(dN)}, {max(dN)}]', flush=True)
            result[name][str((bp, b))] = {str(o): str(cf) for o, cf in rec.items()}
    print('total offsets over all blocks:', tot)
    json.dump(result, open('order2_recurrence.json', 'w'))
