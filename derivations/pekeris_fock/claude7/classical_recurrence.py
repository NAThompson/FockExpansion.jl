"""Global Laguerre recurrence for the classical order-3 block (Codex's 18 atoms x {1, P}), through the phase-2 pipeline
(claude2/order2_recurrence.py: product rule, log(sigma) tower, row clearing with Laguerre finiteness). Stage 1 builds operators;
stage 2 (sharded) expands Laguerre stencils."""
import sys, pickle, json
import sympy as sp
sys.path.insert(0, '../claude2'); sys.path.insert(0, '../sym')
import order2_atoms as oa
import order2_recurrence as orr
from classical_gamma import r1, r2, c, P
from recurrence import operator_to_stencil, stencil_to_recurrence

G = pickle.load(open('classical_gamma.pkl', 'rb'))
P2 = 2 * (r1**2 + r2**2) - c**2
names = [f'g{i}' for i in range(18)]
# sector atoms: g_i and P*g_i
def sector_gamma():
    S = {}
    for i in range(18):
        for rho in (0, 1):
            nm = ('P' if rho else '') + names[i]
            S[nm] = [dict() for _ in range(3)]
            for k in range(3):
                d = S[nm][k]
                if rho:
                    dPP = sp.diff(P2, (r1, r2, c)[k]) / (2 * P2)       # dP/P
                    d[nm] = d.get(nm, 0) + dPP
                for (ii, j, rad), v in G[k].items():
                    if ii != i: continue
                    sig = rad[1]
                    if rho + sig == 2: tgt, cf = names[j], v * P2
                    else: tgt, cf = ('P' if rho + sig else '') + names[j], v
                    d[tgt] = d.get(tgt, 0) + cf
    # keep only sectors reachable from the plain atoms
    reach = set(names); frontier = list(names)
    while frontier:
        a = frontier.pop()
        for k in range(3):
            for b in S[a][k]:
                if b not in reach: reach.add(b); frontier.append(b)
    return {a: S[a] for a in S if a in reach}

if __name__ == '__main__' and sys.argv[1] == 'build':
    S = sector_gamma()
    ATOMS = sorted(S)
    print('sector atoms:', len(ATOMS), ATOMS)
    # map to the order-2 pipeline's expectations: oa.ATOMS and a Gamma table in r-coordinates with keys oa.r1...
    oa.ATOMS = ATOMS
    PI = sp.Symbol('PI')
    orr.FLD, orr.fu, orr.fv, orr.fw, orr.fe, _ = sp.field([orr.u, orr.v, orr.w, orr.eps, PI], sp.QQ)
    orr.FX = (orr.fu, orr.fv, orr.fw); orr.sigma = orr.fu + orr.fv + orr.fw
    def to_sym(fe_):
        q = sp.cancel(fe_.numer.as_expr() / fe_.denom.as_expr())
        assert sp.fraction(q)[1].free_symbols <= {orr.eps, PI}, sp.fraction(q)[1]
        return q
    orr.to_sym = to_sym
    sub = {r1: oa.r1, r2: oa.r2, c: oa.r12, sp.pi: PI}
    Gr = {a: [{t: sp.sympify(v).subs(sub, simultaneous=True) for t, v in S[a][k].items()} for k in range(3)] for a in ATOMS}
    oa.gamma_table = lambda: Gr
    Gm = orr.gamma_uvw()
    O, C1, C2 = orr.build(Gm)
    D = orr.row_denominators(O, C1)
    ops = {}
    for name, blk, mult in (('sigma2P', O, orr.sigma**2), ('C1', C1, orr.FLD(1)), ('C2', C2, orr.FLD(1))):
        for (bp, b), op in blk.items():
            Dr = orr.FLD.from_expr(D[bp].as_expr())
            ops[(name, bp, b)] = {str(kk): str(sp.expand(orr.to_sym(Dr * mult * cf))) for kk, cf in op.items()}
    json.dump(dict(atoms=ATOMS, D={a: str(sp.factor(D[a].as_expr())) for a in ATOMS},
                   ops={str(k): v for k, v in ops.items()}), open('classical_ops.json', 'w'))
    degs = {a: sp.Poly(D[a].as_expr(), orr.u, orr.v, orr.w).total_degree() for a in ATOMS}
    print('blocks:', {n: sum(1 for k in ops if k[0] == n) for n in ('sigma2P', 'C1', 'C2')}, 'row factor degrees:', sorted(set(degs.values())))

if __name__ == '__main__' and sys.argv[1] == 'stencil':
    shard, nshard = int(sys.argv[2]), int(sys.argv[3])
    d = json.load(open('classical_ops.json'))
    keys = sorted(d['ops'])[shard::nshard]
    out = {}
    for k in keys:
        op = {eval(kk): sp.sympify(v, locals={'epsilon': orr.eps, 'u': orr.u, 'v': orr.v, 'w': orr.w, 'PI': sp.Symbol('PI')}) for kk, v in d['ops'][k].items()}
        rec = stencil_to_recurrence(operator_to_stencil(op))
        dN = [sum(o) for o in rec]
        out[k] = dict(n=len(rec), dN=[min(dN), max(dN)], rec={str(o): str(cf) for o, cf in rec.items()})
        print(k, len(rec), [min(dN), max(dN)], flush=True)
    json.dump(out, open(f'classical_stencils_{shard}.json', 'w'))
