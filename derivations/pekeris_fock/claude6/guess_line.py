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

if __name__ == '__main__':
    sel = json.load(open('../codex/phase5/goal2_combined_exact.json'))['results'][0]['selected']
    basis_keys = [tuple(x) for x in sel]
    t0, q0 = 3, 4
    npts = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    xs, vals = [], []
    lam = 0; t_start = time.time()
    while len(xs) < npts:
        lam += 1
        t, q = (t0 + lam) % P, (q0 + 2 * lam) % P
        try:
            rows = rows_upto(t, q, D)
        except (ValueError, AssertionError, ZeroDivisionError) as e:
            continue
        B = [rows[k] for k in basis_keys]
        targets = [k for k in rows if k not in basis_keys and k[1] + k[2] == 8][:3]
        try:
            c = [solve_mod(B, rows[k]) for k in targets]
        except AssertionError:
            continue
        xs.append(lam); vals.append(c)
        print(f'point {len(xs)} lam={lam} ({time.time()-t_start:.0f}s)', flush=True)
        json.dump(dict(xs=xs, vals=vals, targets=[list(k) for k in targets]), open('guess_line.json', 'w'))
    # degree search for a few coefficients
    for ti in range(len(vals[0])):
        for j in (0, 10, 20, 44):
            ys = [v[ti][j] for v in vals]
            found = None
            for tot in range(0, len(xs) - 2):
                for dd in range(0, tot + 1):
                    dn = tot - dd
                    if dn + dd + 2 < len(xs) and ratinterp(xs, ys, dn, dd):
                        found = (dn, dd); break
                if found: break
            print(f'target {ti} coeff {j}: minimal (deg num, deg den) along line = {found}', flush=True)
