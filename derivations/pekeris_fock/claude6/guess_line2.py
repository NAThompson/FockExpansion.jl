"""Goal 2: degree probe. Relation coefficients c(t,q) (rational, branch-invariant) along the line t = t0 + lam, q = q0 + 2 lam over F_P.
For each coefficient find the minimal (deg num, deg den) rational interpolant; report. Jets truncated to order D+1 for speed."""
import sys, json, time
sys.path.insert(0, '../codex/phase5')
import goal2_combined_exact as G
from guess_probe import rows_upto, solve_mod
A = G.A; P = A.P
D = 9
A.N = D + 1
A.IJ = [(i, k - i) for k in range(A.N + 1) for i in range(k + 1)]
A.IX = {v: k for k, v in enumerate(A.IJ)}
A.PAIRS = [[(A.IX[(a, b)], A.IX[(i - a, j - b)]) for a in range(i + 1) for b in range(j + 1)] for i, j in A.IJ]

def ratinterp(xs, ys, dn, dd):
    """find N (deg<=dn), Dn (deg<=dd, monic-ish) with N(x) = y Dn(x) at all points; return True if consistent with len(xs) > dn+dd+1."""
    n = dn + dd + 2
    rows = []
    for x, y in zip(xs, ys):
        rows.append([pow(x, k, P) for k in range(dn + 1)] + [(-y * pow(x, k, P)) % P for k in range(dd + 1)])
    # nullspace dimension of rows (len(xs) x n): rank < n means a solution exists
    M = [r[:] for r in rows]; rank = 0; cols = n
    for c in range(cols):
        pr = next((i for i in range(rank, len(M)) if M[i][c] % P), None)
        if pr is None: continue
        M[rank], M[pr] = M[pr], M[rank]; inv = pow(M[rank][c], -1, P); M[rank] = [v * inv % P for v in M[rank]]
        for i in range(len(M)):
            if i != rank and M[i][c]:
                f = M[i][c]; M[i] = [(a - f * b) % P for a, b in zip(M[i], M[rank])]
        rank += 1
    return rank < n

TARGETS = [(2, 2, 0), (2, 1, 1), (2, 1, 2), (2, 0, 3), (3, 4, 4), (3, 5, 3), (3, 6, 2), (3, 7, 1), (3, 8, 0),
           (3, 4, 5), (3, 3, 6), (3, 2, 7), (3, 1, 8), (3, 0, 9)]
if __name__ == '__main__':
    sel = json.load(open('../codex/phase5/goal2_combined_exact.json'))['results'][0]['selected']
    basis_keys = [tuple(x) for x in sel]
    t0, q0 = 3, 4
    npts = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    OUTF = sys.argv[2] if len(sys.argv) > 2 else 'guess_line.json'
    slope = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    xs, vals = [], []
    lam = 0; t_start = time.time()
    while len(xs) < npts:
        lam += 1
        t, q = (t0 + lam) % P, (q0 + slope * lam) % P
        try:
            rows = rows_upto(t, q, D)
        except (ValueError, AssertionError, ZeroDivisionError) as e:
            continue
        B = [rows[k] for k in basis_keys]
        targets = TARGETS
        try:
            c = [solve_mod(B, rows[k]) for k in targets]
        except AssertionError:
            continue
        xs.append(lam); vals.append(c)
        print(f'point {len(xs)} lam={lam} ({time.time()-t_start:.0f}s)', flush=True)
        json.dump(dict(xs=xs, vals=vals, targets=[list(k) for k in targets]), open(OUTF, 'w'))
