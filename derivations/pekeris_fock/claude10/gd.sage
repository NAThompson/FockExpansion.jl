# Griffiths-Dwork / creative-telescoping ansatz over GF(p) at a fixed target-line parameter lam0:
#   sum_{m<=r} c_m d_lam^m F^{-1/2} = d_X(U/F^{r-1/2}) + d_q(V/F^{r-1/2}),  F = X^3 + 2B X^2 + A C X.
# Finds the minimal r with c != 0 and returns c (normalised c_r = 1).
import sys, time, os
which = sys.argv[1]; r_max = int(sys.argv[2]); p = Integer(os.environ.get("PRIME", 2^31 - 1))
def setup(which):
    R.<X,q,l> = GF(p)[]
    D = 1 + q^4; P1, P2 = q^4 - 1, 2*q^2
    a, line = which.split(':')        # e.g. 'e1:x1=l,x2=1/3' ; lines are affine in l
    env = {'l': l, 'R': R}
    x1 = eval(line.split(',')[0].split('=')[1].replace('/', '*R(1)/'), {}, env) if True else None
    x2 = eval(line.split(',')[1].split('=')[1].replace('/', '*R(1)/'), {}, env)
    ax, aP = {'e2': (x2, P2), 'e1': (x1, P1), 'me2': (-x2, -P2), 'me1': (-x1, -P1)}[a]
    Px = P1*x1 + P2*x2
    A = 1 - ax^2; B = aP - ax*Px; C = D^2 - Px^2
    return R, X^3 + 2*B*X^2 + A*C*X
R, F = setup(which)
X, q, l = R.gens()
N = [R(1)]
for m in range(r_max):
    N.append(N[m].derivative(l)*F - (m + R(1)/2)*N[m]*F.derivative(l))
def solve_at(lam0, r, dUX, dUq):
    k = GF(p); S.<X,q> = k[]
    phi = R.hom([X, q, k(lam0)], S)
    sub = lambda f: phi(f)
    Fs = sub(F); FX = Fs.derivative(X); Fq = Fs.derivative(q)
    cols = []
    for m in range(r + 1):
        cols.append(sub(N[m]) * Fs^(r - m))
    s = k(r) - k(1)/2
    mons = [(i, j) for i in range(dUX + 1) for j in range(dUq + 1)]
    for (i, j) in mons:
        u = X^i*q^j
        cols.append(-(u.derivative(X)*Fs - s*u*FX))
    for (i, j) in mons:
        v = X^i*q^j
        cols.append(-(v.derivative(q)*Fs - s*v*Fq))
    allm = sorted(set(sum([c.monomials() for c in cols], [])))
    idx = {mm: t for t, mm in enumerate(allm)}
    M = matrix(k, len(allm), len(cols), sparse=True)
    for jc, c in enumerate(cols):
        for coef, mm in zip(c.coefficients(), c.monomials()):
            M[idx[mm], jc] = coef
    K = M.right_kernel_matrix()
    # project kernel to the c-part
    Cpart = K[:, :r + 1]
    if Cpart.rank() == 0:
        return None
    E = Cpart.echelon_form()
    return [row for row in E.rows() if row != 0], M.dimensions()
EXTRA = int(sys.argv[4]) if len(sys.argv) > 4 else 0
if __name__ == '__main__' and not globals().get('NO_MAIN'):
    lam0 = QQ(sys.argv[3]) if len(sys.argv) > 3 else QQ(2)/7
    for r in range(1, r_max + 1):
        t = time.time()
        res = solve_at(lam0, r, 3*r + EXTRA, 8*r + 2*EXTRA)
        print('r =', r, 'result:', None if res is None else (res[1], len(res[0]), res[0][-1]), ' %.1fs' % (time.time() - t)); sys.stdout.flush()
        if res is not None: break
