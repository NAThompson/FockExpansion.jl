"""Phase 7: per-block (classical; panel j) cyclic bases and the line degrees / branch behaviour of their closure relations.
Reuses Codex's exact finite-field code read-only (codex/phase5/goal2_combined_exact.py) and claude6 helpers."""
import sys, json, time
sys.path.insert(0, '../codex/phase5'); sys.path.insert(0, '../claude6')
import goal2_combined_exact as G
from guess_probe import solve_mod, SIGNS, CALL   # also installs the signed sqrt hook
A = G.A; P = A.P
D = int(os.environ.get('D', '6')) if (os := __import__('os')) else 6
A.N = D + 1
A.IJ = [(i, k - i) for k in range(A.N + 1) for i in range(k + 1)]
A.IX = {v: k for k, v in enumerate(A.IJ)}
A.PAIRS = [[(A.IX[(a, b)], A.IX[(i - a, j - b)]) for a in range(i + 1) for b in range(j + 1)] for i, j in A.IJ]

def piece_rows(t, q, D):
    L, C, T, pieces = G.build(t, q, True)
    out = {}
    for piece in pieces:
        prev = {(0, 0): piece['rows']}; rows_by = {}
        for order in range(D + 1):
            G.WORK = A.N - order; cur = {}
            for a in range(order + 1):
                b = order - a
                if not order: rows = piece['rows']
                else:
                    axis = 0 if a else 1; old = prev[a - 1, b] if a else prev[a, b - 1]; rows = []
                    for row in old:
                        o = {i: v.diff(axis) for i, v in row.items()}
                        for i, v in row.items():
                            for j, w in C[axis].get(i, []): o[j] = o.get(j, G.J()) + v * w
                        rows.append(o)
                cur[a, b] = rows
                for src, row in enumerate(rows): rows_by[(src, a, b)] = G.project(row, T)
            prev = cur
        out[piece['name']] = rows_by
    return out

def select_basis(rows_by):
    basis = []; keys = []
    for k in sorted(rows_by, key=lambda k: (k[1] + k[2], k)):
        if A.rankadd(basis, list(rows_by[k])): keys.append(k)
    return keys

def closure_targets(keys, D):
    s = set(keys); tg = []
    for (src, a, b) in keys:
        for nk in ((src, a + 1, b), (src, a, b + 1)):
            if nk not in s and nk not in tg and nk[1] + nk[2] <= D: tg.append(nk)
    return tg

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'branch':
        t, q = 3, 4; res = {}
        for signs in [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)]:
            SIGNS[:] = signs; CALL[0] = 0
            pr = piece_rows(t, q, D)
            res[signs] = pr
        ref = res[(1, 1, 1)]
        for name in ref:
            keys = select_basis(ref[name]); tg = closure_targets(keys, D)
            Bm = [ref[name][k] for k in keys]
            c0 = {k: solve_mod(Bm, ref[name][k]) for k in tg}
            for signs, pr in res.items():
                try:
                    Bs = [pr[name][k] for k in keys]
                    cs = {k: solve_mod(Bs, pr[name][k]) for k in tg}
                    ch = sum(1 for k in tg for a, b in zip(c0[k], cs[k]) if a != b)
                except AssertionError:
                    ch = 'basis-dependent'
                print(name, 'rank', len(keys), 'targets', len(tg), 'branch', signs, 'changed', ch, flush=True)
    else:  # line data: python block_line.py line NPTS OUT start step
        npts, outf, start, step = int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
        ref = piece_rows(3, 4, D)
        meta = {name: (select_basis(ref[name]),) for name in ref}
        meta = {name: (m[0], closure_targets(m[0], D)) for name, m in meta.items()}
        xs = []; vals = {name: [] for name in ref}; lam = start; t0 = time.time()
        while len(xs) < npts:
            lam += step
            t, q = (3 + lam) % P, (4 + 2 * lam) % P
            try:
                pr = piece_rows(t, q, D)
                cv = {name: [solve_mod([pr[name][k] for k in meta[name][0]], pr[name][k]) for k in meta[name][1]] for name in pr}
            except (ValueError, AssertionError, ZeroDivisionError):
                continue
            xs.append(lam)
            for name in pr: vals[name].append(cv[name])
            print('point', len(xs), f'{time.time()-t0:.0f}s', flush=True)
            json.dump(dict(xs=xs, vals=vals, meta={n: [[list(k) for k in m[0]], [list(k) for k in m[1]]] for n, m in meta.items()}), open(outf, 'w'))
