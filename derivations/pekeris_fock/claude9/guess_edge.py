# Guess a linear ODE sum_{i<=r, j<=d} c_ij tau^j f^(i)(tau) = 0 for an edge restriction, from series data (edge_series.txt).
import sys
from mpmath import mp, mpf, matrix, svd, factorial, nstr
mp.dps = 64
name = sys.argv[1]; maxunk = int(sys.argv[2]) if len(sys.argv) > 2 else 120
rows = []
for l in open('edge_series.txt'):
    p = l.split()
    if p[0] == name:
        rows.append((mpf(p[1]), [mpf(x) * factorial(i) for i, x in enumerate(p[2:])]))
best = None
for tot in range(2, maxunk + 1):
    for r in range(1, 13):
        if tot % (r + 1): continue
        d = tot // (r + 1) - 1
        if d < 0: continue
        A = matrix([[t**j * fd[i] for i in range(r + 1) for j in range(d + 1)] for t, fd in rows])
        S = svd(A, compute_uv=False)
        s = sorted([S[k] for k in range(len(S))])
        ratio = s[0] / s[1]
        if ratio < mpf(10)**-20:
            print(f'{name}: ODE found r={r} d={d} unknowns={tot} smin={nstr(s[0],3)} next={nstr(s[1],3)}'); sys.exit()
print(name, 'no ODE with <=', maxunk, 'unknowns')
