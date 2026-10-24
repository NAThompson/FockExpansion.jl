"""Goal 2 probe: are the order-8 relation coefficients (boundary derivative rows in terms of the 45 basis rows) RATIONAL in (t,q)?
Test: recompute them at the same point with flipped branches of the panel radicals (disc = sqrt D, z1, z3). Rational (Galois-invariant)
coefficients must not change. Uses Codex's exact finite-field Taylor code (codex/phase5/goal2_combined_exact.py), read-only import."""
import sys, json, itertools
sys.path.insert(0, '../codex/phase5')
import goal2_combined_exact as G
A = G.A; P = A.P
_orig_sqrtj = G.sqrtj
SIGNS = [1, 1, 1]; CALL = [0]
def sqrtj_signed(x):
    r = _orig_sqrtj(x)
    s = SIGNS[CALL[0] % 3]; CALL[0] += 1
    return r if s == 1 else -r
G.sqrtj = sqrtj_signed

def rows_upto(t, q, D):
    L, C, T = G.build(t, q)
    previous = {(0, 0): L}; out = {}
    for order in range(D + 1):
        G.WORK = A.N - order; current = {}
        for a in range(order + 1):
            b = order - a
            if not order: rows = L
            else:
                axis = 0 if a else 1; old = previous[a - 1, b] if a else previous[a, b - 1]; rows = []
                for row in old:
                    o = {i: v.diff(axis) for i, v in row.items()}
                    for i, v in row.items():
                        for j, w in C[axis].get(i, []): o[j] = o.get(j, G.J()) + v * w
                    rows.append(o)
            current[a, b] = rows
            for src, row in enumerate(rows): out[(src + 1, a, b)] = G.project(row, T)
        previous = current
    return out

def solve_mod(Bmat, v):
    """coefficients c with sum_j c_j Bmat[j] = v (rows of length 240), mod P; Bmat rows independent."""
    m = len(Bmat); n = len(v)
    M = [list(Bmat[j]) for j in range(m)]
    # solve M^T c = v by elimination on augmented [M^T | v]
    Aug = [[M[j][i] for j in range(m)] + [v[i]] for i in range(n)]
    piv = []; r = 0
    for col in range(m):
        pr = next((i for i in range(r, n) if Aug[i][col] % P), None)
        assert pr is not None, 'basis rows dependent'
        Aug[r], Aug[pr] = Aug[pr], Aug[r]
        inv = pow(Aug[r][col], -1, P); Aug[r] = [x * inv % P for x in Aug[r]]
        for i in range(n):
            if i != r and Aug[i][col]:
                f = Aug[i][col]; Aug[i] = [(x - f * y) % P for x, y in zip(Aug[i], Aug[r])]
        r += 1
    assert all(Aug[i][m] == 0 for i in range(r, n)), 'v not in span'
    return [Aug[i][m] for i in range(m)]

if __name__ == '__main__':
    sel = json.load(open('../codex/phase5/goal2_combined_exact.json'))['results'][0]['selected']
    basis_keys = [tuple(x) for x in sel]
    t, q = 3, 4
    res = {}
    for signs in [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1), (-1, -1, -1)]:
        SIGNS[:] = signs; CALL[0] = 0
        try:
            rows = rows_upto(t, q, 9)
        except ValueError as e:
            print('skip', e); continue
        B = [rows[k] for k in basis_keys]
        targets = [k for k in rows if k not in basis_keys and k[1] + k[2] in (8, 9)][:6]
        res[signs] = {k: solve_mod(B, rows[k]) for k in targets}
        print('signs', signs, 'targets', targets[:3], flush=True)
    ref = res[(1, 1, 1)]
    for signs, r in res.items():
        diff = sum(1 for k in ref for a, b in zip(ref[k], r[k]) if a != b)
        print('branch', signs, ': number of changed coefficients vs reference =', diff, 'of', sum(len(v) for v in ref.values()))
