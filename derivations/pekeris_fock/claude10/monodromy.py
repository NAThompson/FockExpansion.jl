# Numerical analytic continuation for a linear ODE sum_m c_m(l) D^m y = 0 with polynomial coefficients (Taylor stepping, mpmath).
from mpmath import mp, mpf, mpc, matrix, binomial, factorial, exp, pi, sqrt
mp.dps = 50
def polyshift(c, z0):
    # coefficients of c(z0 + t) in powers of t; c = list low->high
    n = len(c); out = [mpc(0)] * n
    for i, a in enumerate(c):
        for k in range(i + 1):
            out[k] += a * binomial(i, k) * z0**(i - k)
    return out
def series_solution(ops, z0, jet, N):
    # ops[m] = coefficient list of c_m; jet = [y, y', ..., y^(r-1)] at z0; returns Taylor coeffs a_n (n < N)
    r = len(ops) - 1
    sh = [polyshift(c, z0) for c in ops]
    a = [jet[k] / factorial(k) for k in range(r)] + [mpc(0)] * (N - r)
    # sum_m sum_i sh[m][i] t^i * sum_n a_n n!/(n-m)! t^(n-m) = 0 ; coefficient of t^K: solve for a_{K+r}
    lead = sh[r][0]
    for K in range(0, N - r):
        s = mpc(0)
        for m in range(r + 1):
            for i, ci in enumerate(sh[m]):
                n = K - i + m
                if n < 0 or n >= N or (m == r and i == 0): continue
                if n < m: continue
                s += ci * a[n] * factorial(n) / factorial(n - m)
        a[K + r] = -s / (lead * factorial(K + r) / factorial(K))
    return a
def eval_jet(a, t, r):
    out = []
    for m in range(r):
        s = mpc(0)
        for n in range(m, len(a)):
            s += a[n] * factorial(n) / factorial(n - m) * t**(n - m)
        out.append(s)
    return out
def transport(ops, path, jet, N=80, step_frac=0.4, sing=()):
    # move jet along polyline 'path' (list of complex points), steps limited by distance to singularities
    r = len(ops) - 1
    z = mpc(path[0]); y = [mpc(v) for v in jet]
    for target in path[1:]:
        target = mpc(target)
        while abs(target - z) > mpf(10)**(-40):
            dist = min(abs(z - s) for s in sing) if sing else mpf(1)
            h = target - z
            if abs(h) > step_frac * dist: h = h / abs(h) * step_frac * dist
            a = series_solution(ops, z, y, N)
            y = eval_jet(a, h, r); z = z + h
    return y
def monodromy_matrix(ops, base, around, radius, sing, N=80):
    # loop: base -> around - radius (along straight line) -> circle ccw -> back
    r = len(ops) - 1
    M = matrix(r, r)
    start = mpc(around) + radius * (mpc(base) - mpc(around)) / abs(mpc(base) - mpc(around))
    circle = [mpc(around) + (start - mpc(around)) * exp(2j * pi * k / 16) for k in range(17)]
    for col in range(r):
        jet = [mpc(1) if k == col else mpc(0) for k in range(r)]
        y = transport(ops, [base, start] + circle[1:] + [base], jet, N, sing=sing)
        for row in range(r): M[row, col] = y[row]
    return M
