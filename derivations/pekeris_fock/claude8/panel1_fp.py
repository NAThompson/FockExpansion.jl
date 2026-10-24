"""Panel 1 (and panel-type 1/2) rational presentation, pointwise over F_p.

Integrand family: f(q; t, s) = exp(t A(q) + s D(q)) / W(q), A = atan q, D = atan(q/h), W = sqrt(T),
T(q) = (1-b) q^4 + (4 Y2 - 2) q^2 + 1 + b,  b = (h^2-1)/(h^2+1),  Y2 = y^2.   (panel 1: h = r1/r2, Y2 = r12^2/(2R^2)).
Fixed upper endpoint 1 (Codex, phase 7: physical N0, N2 equal their fixed-endpoint-1 integrals).
States (per family F in {A: t-jets, D: s-jets}, jets n = 0,1,2 as Taylor coefficients):
  J = int f, K = int q f, J_i = int f/(q - z_i) (i = 1..8), Js = FP int f/q, e = f(1), f0 = f(0).
Poles z_1..4 = (z, -z, h/z, -h/z) roots of T, z_5..8 = (i, -i, i h, -i h); exponents a = (-1/2 x4, t/2i, -t/2i, s/2i, -s/2i).
Connection (Codex connection.md with u = 1 fixed, u_x = 0), for x in {h, Y2}:
  dJ = k J - sum a_i v_i J_i,  dK = k K - sum a_i v_i (J + z_i J_i),
  dJ_i = k J_i + v_i H_ii - sum_j a_j v_j H_ij,  dJs = k Js + sum (a_i v_i / z_i)(Js - J_i),
  de = [k - sum a_i v_i/(1 - z_i)] e,  df0 = l0 f0,
  k = l0 - sum a_i v_i / z_i,  H_ij = (J_i - J_j)/(z_i - z_j),  H_ii = [e/(1-z_i) + f0/z_i - sum_{j!=i} a_j H_ij]/(a_i - 1).
Galois-symmetric rebasing: S_k = sum_{i=1..4} z_i^k J_i (k = 0..3). Everything else unchanged. New connection C' = (dV) V^-1 + V C V^-1.
"""
import sys
P = 998244353

def inv(a): return pow(a % P, P - 2, P)
def sqrt_mod(a):
    import sympy
    r = sympy.sqrt_mod(a % P, P, all_roots=True)
    if not r: raise ValueError('nonsplit')
    return min(r)

I = pow(3, (P - 1) // 4, P); I0 = I
assert I * I % P == P - 1
# --- jets: Taylor polynomials of length 3 in the family parameter (t or s) ---
def jmul(x, y): return [(x[0]*y[0]) % P, (x[0]*y[1] + x[1]*y[0]) % P, (x[0]*y[2] + x[1]*y[1] + x[2]*y[0]) % P]
def jadd(x, y): return [(a + b) % P for a, b in zip(x, y)]
def jsc(x, c): return [(a * c) % P for a in x]
def jinv(x):
    v0 = inv(x[0]); v1 = (-x[1] * v0 * v0) % P; v2 = (-(x[1] * v1 + x[2] * v0) * v0) % P
    return [v0, v1, v2]
def jc(c): return [c % P, 0, 0]

NS = 13  # J, K, J1..J8, Js, e, f0

def geometry(h, Y2, root_choice=0):
    b = (h*h - 1) * inv(h*h + 1) % P
    c4, c2, c0 = (1 - b) % P, (4*Y2 - 2) % P, (1 + b) % P
    disc = sqrt_mod(c2*c2 - 4*c4*c0)
    Z2 = [(-c2 + disc) * inv(2*c4) % P, (-c2 - disc) * inv(2*c4) % P]
    z = sqrt_mod(Z2[root_choice % 2])
    if root_choice >= 2: z = (-z) % P
    zs = [z, (-z) % P, h * inv(z) % P, (-h * inv(z)) % P]
    assert all((c4*pow(q, 4, P) + c2*q*q + c0) % P == 0 for q in zs)
    return b, (c4, c2, c0), zs + [I, (-I) % P, I*h % P, (-I*h) % P]

def connection(h, Y2, axis, family, root_choice=0):
    """Matrix C (39x39, states x jets; index = state*3 + jet) with dX/dx = C X for one family; axis 0 = h, 1 = Y2."""
    b, (c4, c2, c0), z = geometry(h, Y2, root_choice)
    # parameter derivatives of the quartic coefficients
    if axis == 0:   # d/dh: b_h = 4h/(h^2+1)^2
        bx = 4*h * inv(pow(h*h + 1, 2, P)) % P
        T_x = lambda q: (bx * (1 - pow(q, 4, P))) % P
        hx = 1
        l0 = (-inv(2*(1 + b)) * bx) % P          # f0 = (1+b)^-1/2
    else:
        bx = 0; T_x = lambda q: (4*q*q) % P; hx = 0; l0 = 0
    Tq = lambda q: (4*c4*pow(q, 3, P) + 2*c2*q) % P
    v = [(-T_x(q) * inv(Tq(q))) % P for q in z[:4]] + [0, 0, I*hx % P, (-I*hx) % P]
    tj = [0, 1, 0] if family == 'A' else [0, 0, 0]
    sj = [0, 1, 0] if family == 'D' else [0, 0, 0]
    i2 = inv(2*I)
    a = [jc(-inv(2)) for _ in range(4)] + [jsc(tj, i2), jsc(tj, -i2), jsc(sj, i2), jsc(sj, -i2)]
    kap = jc(l0)
    for ai, zi, vi in zip(a, z, v): kap = jadd(kap, jsc(ai, -vi * inv(zi)))
    n = NS * 3
    C = [[0]*n for _ in range(n)]
    def addjet(row_state, col_state, coef):   # contribution coef (a jet) * X[col_state] to d X[row_state], jet-convolved
        for jr in range(3):
            for jc_ in range(jr + 1):
                C[row_state*3 + jr][col_state*3 + jc_] = (C[row_state*3 + jr][col_state*3 + jc_] + coef[jr - jc_]) % P
    J, K, Js, E, F0 = 0, 1, 10, 11, 12
    Ji = lambda i: 2 + i
    # dJ
    addjet(J, J, kap)
    for i in range(8): addjet(J, Ji(i), jsc(a[i], -v[i]))
    # dK
    addjet(K, K, kap)
    for i in range(8):
        addjet(K, J, jsc(a[i], -v[i])); addjet(K, Ji(i), jsc(a[i], -v[i] * z[i]))
    # dJ_i: H_ij = (J_i - J_j)/(z_i - z_j); H_ii = [E/(1-z_i) + F0/z_i - sum_j a_j H_ij] / (a_i - 1)
    for i in range(8):
        addjet(Ji(i), Ji(i), kap)
        w = jc(-1)   # (1 - a_i) * 1/(a_i - 1): the -a_i v_i H_ii term is included (bug fix)
        # v_i H_ii
        addjet(Ji(i), E, jsc(w, v[i] * inv(1 - z[i])))
        addjet(Ji(i), F0, jsc(w, v[i] * inv(z[i])))
        for j in range(8):
            if j == i: continue
            d = inv(z[i] - z[j])
            # - v_i/(a_i-1) * a_j H_ij   and   - a_j v_j H_ij
            coef = jadd(jsc(jmul(w, a[j]), -v[i]), jsc(a[j], -v[j]))
            addjet(Ji(i), Ji(i), jsc(coef, d)); addjet(Ji(i), Ji(j), jsc(coef, -d))
    # dJs
    addjet(Js, Js, kap)
    for i in range(8):
        f = jsc(a[i], v[i] * inv(z[i]))
        addjet(Js, Js, f); addjet(Js, Ji(i), jsc(f, -1))
    # de, df0
    ep = kap
    for i in range(8): ep = jadd(ep, jsc(a[i], -v[i] * inv(1 - z[i])))
    addjet(E, E, ep)
    addjet(F0, F0, jc(l0))
    return C, z, v

def sym_transform(h, Y2, axis, family, root_choice=0):
    """C' = dV V^-1 + V C V^-1 in the basis (J, K, S0..S3, J5..J8, Js, e, f0) x jets."""
    C, z, v = connection(h, Y2, axis, family, root_choice)
    n = NS * 3
    V = [[0]*n for _ in range(n)]; dV = [[0]*n for _ in range(n)]
    for st in range(NS):
        for jt in range(3):
            r = st*3 + jt
            if 2 <= st <= 5:   # S_k, k = st-2, combines J_1..J_4 at same jet
                k = st - 2
                for i in range(4):
                    V[r][(2 + i)*3 + jt] = pow(z[i], k, P)
                    dV[r][(2 + i)*3 + jt] = (k * pow(z[i], k - 1, P) * v[i]) % P if k else 0
            else:
                V[r][r] = 1
    Vi = matinv(V)
    A1 = matmul(dV, Vi); A2 = matmul(matmul(V, C), Vi)
    return [[(x + y) % P for x, y in zip(r1, r2)] for r1, r2 in zip(A1, A2)]

def matmul(A, B):
    n, m, k = len(A), len(B[0]), len(B)
    Bt = list(zip(*B))
    return [[sum(A[i][t] * Bt[j][t] for t in range(k)) % P for j in range(m)] for i in range(n)]
def matinv(A):
    n = len(A); M = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        pr = next(i for i in range(c, n) if M[i][c] % P)
        M[c], M[pr] = M[pr], M[c]; iv = inv(M[c][c]); M[c] = [x * iv % P for x in M[c]]
        for i in range(n):
            if i != c and M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % P for x, y in zip(M[i], M[c])]
    return [row[n:] for row in M]

if __name__ == '__main__':
    # Galois test: the symmetric connection must not depend on which quartic root is called z
    h, Y2 = 5, 7
    for axis in (0, 1):
        for fam in ('A', 'D'):
            ref = None
            for rc in range(4):
                try:
                    Cs = sym_transform(h, Y2, axis, fam, rc)
                except ValueError:
                    print('nonsplit point'); break
                if ref is None: ref = Cs
                else: print('axis', axis, 'family', fam, 'root choice', rc, 'entries changed:', sum(1 for r1, r2 in zip(ref, Cs) for x, y in zip(r1, r2) if x != y))
