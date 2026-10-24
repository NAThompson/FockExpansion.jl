"""Exact reconstruction of the Galois-symmetric panel-1 connection over Q(i)(h, Y2):
real/imaginary parts via the two embeddings of i, bivariate dense rational interpolation per entry at each prime,
CRT over several primes, rational number reconstruction, validation on a fresh prime. Output: panel1_connection.json
{axis: {family: {"r,c": [re_num, re_den, im_num, im_den]}}} with polynomials as {"i j": "p/q"} (h^i Y2^j)."""
import json, random, sys
import sympy as sp
import panel1_fp as F

PRIMES = [998244353, 1004535809, 469762049, 167772161, 754974721]   # all = 1 mod 4

def set_prime(p):
    F.P = p
    I = pow(3, (p - 1) // 4, p)
    if I * I % p != p - 1:
        g = next(g for g in range(2, 50) if pow(g, (p - 1) // 2, p) == p - 1)
        I = pow(g, (p - 1) // 4, p)
    assert I * I % p == p - 1
    F.I0 = I

def eval_parts(h, Y2, p):
    """dict (axis, family) -> (re, im) matrices over F_p at (h, Y2)."""
    set_prime(p)
    out = {}
    for ax in (0, 1):
        for fam in ('A', 'D'):
            F.I = F.I0; M1 = F.sym_transform(h, Y2, ax, fam)
            F.I = (p - F.I0) % p; M2 = F.sym_transform(h, Y2, ax, fam)
            inv2 = pow(2, p - 2, p); inv2i = pow(2 * F.I0 % p, p - 2, p)
            re = [[(a + b) * inv2 % p for a, b in zip(r1, r2)] for r1, r2 in zip(M1, M2)]
            im = [[(a - b) * inv2i % p for a, b in zip(r1, r2)] for r1, r2 in zip(M1, M2)]
            out[(ax, fam)] = (re, im)
    return out

def monos(d): return [(i, t - i) for t in range(d + 1) for i in range(t + 1)]

def nullvec(rows, ncols, p):
    M = [r[:] for r in rows]; piv = []; r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(M)) if M[i][c] % p), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]; iv = pow(M[r][c], p - 2, p); M[r] = [x * iv % p for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
        piv.append(c); r += 1
    free = [c for c in range(ncols) if c not in piv]
    if len(free) != 1: return None, len(free)
    f = free[0]; vec = [0] * ncols; vec[f] = 1
    for i, c in enumerate(piv): vec[c] = (-M[i][f]) % p
    return vec, 1

def fit(points, values, dn, dd, p):
    mn, md = monos(dn), monos(dd)
    rows = []
    for (h, y), val in zip(points, values):
        rows.append([pow(h, a, p) * pow(y, b, p) % p for a, b in mn] + [(-val * pow(h, a, p) * pow(y, b, p)) % p for a, b in md])
    vec, k = nullvec(rows, len(mn) + len(md), p)
    if vec is None: return None
    num, den = vec[:len(mn)], vec[len(mn):]
    # normalize: first nonzero denominator coefficient (in monomial order) = 1
    j = next(i for i, x in enumerate(den) if x)
    iv = pow(den[j], p - 2, p)
    return [x * iv % p for x in num], [x * iv % p for x in den], j

def ratrec_int(a, m):
    """rational number reconstruction of a mod m (Wang)."""
    import math
    r0, r1, s0, s1 = m, a % m, 0, 1
    bound = math.isqrt(m // 2)
    while r1 > bound:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound: return None
    from fractions import Fraction
    return Fraction(r1 * (1 if s1 > 0 else -1), abs(s1))

if __name__ == '__main__':
    random.seed(1)
    npts = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    data = {}
    for p in PRIMES:
        pts, vals = [], []
        while len(pts) < npts:
            h, y = random.randrange(2, p - 1), random.randrange(2, p - 1)
            try:
                vals.append(eval_parts(h, y, p)); pts.append((h, y))
            except (ValueError, StopIteration, ZeroDivisionError):
                continue
        data[p] = (pts, vals)
        print('prime', p, 'points', len(pts), flush=True)
    import pickle; pickle.dump(data, open('panel1_points.pkl', 'wb'))
