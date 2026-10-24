"""Univariate rational reconstruction mod P via the extended Euclidean algorithm (Pade / Wang), with held-out validation."""
P = 998244353

def trim(a):
    while a and a[-1] == 0: a.pop()
    return a
def padd(a, b):
    n = max(len(a), len(b)); return trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % P for i in range(n)])
def psub(a, b): return padd(a, [(-x) % P for x in b])
def pmul(a, b):
    if not a or not b: return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): out[i + j] = (out[i + j] + x * y) % P
    return trim(out)
def pdivmod(a, b):
    a = a[:]; q = [0] * max(0, len(a) - len(b) + 1); inv = pow(b[-1], -1, P)
    while len(a) >= len(b) and a:
        c = a[-1] * inv % P; d = len(a) - len(b); q[d] = c
        for i, y in enumerate(b): a[i + d] = (a[i + d] - c * y) % P
        trim(a)
    return trim(q), a
def peval(a, x):
    r = 0
    for c in reversed(a): r = (r * x + c) % P
    return r
def interp(xs, ys):
    """Newton interpolation -> coefficient list."""
    n = len(xs); c = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) * pow((xs[i] - xs[i - j]) % P, -1, P) % P
    poly = [c[-1]]
    for i in range(n - 2, -1, -1):
        poly = padd(pmul(poly, [(-xs[i]) % P, 1]), [c[i]])
    return poly
def ratrec(xs, ys, holdout=8):
    """Return (deg N, deg D, N, D) of the minimal rational function through the first n-holdout points that also fits the holdout, or None."""
    xa, ya = xs[:-holdout], ys[:-holdout]
    M = [1]
    for x in xa: M = pmul(M, [(-x) % P, 1])
    r0, r1 = M, interp(xa, ya); t0, t1 = [], [1]
    best = None
    while True:
        # candidate N = r1, D = t1
        if t1 and all(peval(t1, x) for x in xs):
            ok = all(peval(r1, x) * pow(peval(t1, x), -1, P) % P == y % P for x, y in zip(xs, ys))
            if ok:
                cand = (len(r1) - 1 if r1 else -1, len(t1) - 1, r1, t1)
                if best is None or cand[0] + cand[1] < best[0] + best[1]: best = cand
        if not r1: break
        q, r = pdivmod(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, psub(t0, pmul(q, t1))
    return best
