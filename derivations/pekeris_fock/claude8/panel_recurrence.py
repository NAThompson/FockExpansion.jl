"""Global Laguerre recurrence for one elliptic panel (Codex's reduced 23-state rational distance connection, codex/phase8/x16_connection.txt),
through the phase-2 pipeline (claude2/order2_recurrence.py), pi kept symbolic as PI. Writes panel{j}_ops.json and exact term lists."""
import sys, json
import sympy as sp
sys.path.insert(0, '../claude2'); sys.path.insert(0, '../sym')
import order2_atoms as oa
import order2_recurrence as orr
j = sys.argv[1]
x, y, c = sp.symbols('r1 r2 r12')
loc = {'r1': x, 'r2': y, 'r12': c, 'pi': sp.Symbol('PI')}
G = {}
for line in open('../codex/phase8/x16_connection.txt'):
    pnl, axis, row, col, num, den = line.rstrip('\n').split('|')
    if pnl != j: continue
    k = {'r1': 0, 'r2': 1, 'r12': 2}[axis]
    G.setdefault(row, [dict(), dict(), dict()])
    G.setdefault(col, [dict(), dict(), dict()])
    G[row][k][col] = sp.sympify(num, locals=loc) / sp.sympify(den, locals=loc)
ATOMS = sorted(G)
print('panel', j, 'atoms', len(ATOMS))
oa.ATOMS = ATOMS
PI = sp.Symbol('PI')
orr.FLD, orr.fu, orr.fv, orr.fw, orr.fe, _ = sp.field([orr.u, orr.v, orr.w, orr.eps, PI], sp.QQ)
orr.FX = (orr.fu, orr.fv, orr.fw); orr.sigma = orr.fu + orr.fv + orr.fw
def to_sym(fe_):
    q = sp.cancel(fe_.numer.as_expr() / fe_.denom.as_expr())
    assert sp.fraction(q)[1].free_symbols <= {orr.eps, PI}, sp.fraction(q)[1]
    return q
orr.to_sym = to_sym
sub = {x: oa.r1, y: oa.r2, c: oa.r12}
Gr = {a: [{t: sp.sympify(v).subs(sub, simultaneous=True) for t, v in G[a][k].items()} for k in range(3)] for a in ATOMS}
oa.gamma_table = lambda: Gr
Gm = orr.gamma_uvw()
O, C1, C2 = orr.build(Gm)
D = orr.row_denominators(O, C1)
lines = []
u, v, w = orr.u, orr.v, orr.w
nblk = {}
for name, blk, mult in (('sigma2P', O, orr.sigma**2), ('C1', C1, orr.FLD(1)), ('C2', C2, orr.FLD(1))):
    for (bp, b), op in blk.items():
        Dr = orr.FLD.from_expr(D[bp].as_expr())
        key = str((name, f'p{j}:{bp}', f'p{j}:{b}'))
        nblk[name] = nblk.get(name, 0) + 1
        for kk, cf in op.items():
            poly = sp.Poly(sp.expand(orr.to_sym(Dr * mult * cf)), u, v, w, orr.eps, PI)
            for mono, cc in poly.terms():
                cc = sp.Rational(cc)
                lines.append(f"{key}|{' '.join(map(str, kk))}|{' '.join(map(str, mono))}|{cc.p}|{cc.q}")
open(f'panel{j}_terms.txt', 'w').write('\n'.join(lines) + '\n')
json.dump({a: str(sp.factor(D[a].as_expr())) for a in ATOMS}, open(f'panel{j}_rowfactors.json', 'w'))
print('blocks', nblk, 'row factor degrees', sorted({sp.Poly(D[a].as_expr(), u, v, w).total_degree() for a in ATOMS}))
