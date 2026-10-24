exec(open('verify41.py').read().split("c = {l.split")[0])
from mpmath import mpf, sqrt, sin, cos, matrix, lu_solve, mp
import mpmath
D = [list(map(mpf, l.split())) for l in open('data2.txt')]
def basis(a, t, kind):
    x1, x2 = cos(a / 2), sin(a / 2); x12 = sqrt(2 * sin(mpi / 4 - a / 2)**2 + 2 * sin(a) * sin(t / 2)**2)
    Y40 = 4 * (x1**2 - x2**2)**2 - 1; Y42 = sin(a)**2 * (3 * cos(t)**2 - 1) / 2
    if kind == 'harm': return [Y40, Y42]
    if kind == 'harm+1': return [Y40, Y42, mpf(1)]
    return [x1**i * x2**j * x12**k for j in (0, 1) for i in range(5) for k in range(5) if i + j + k <= 4]
y = [v - F(a, t, 0, 0) for a, t, v in D]
for kind in ('harm', 'harm+1', 'poly4'):
    rows = [basis(a, t, kind) for a, t, v in D]
    M = matrix(rows); Y = matrix(y)
    s = lu_solve(M.T * M, M.T * Y)
    r = max(abs(y[i] - sum(s[j] * rows[i][j] for j in range(len(rows[0])))) for i in range(len(y)))
    print(kind, mp.nstr(r, 4), [mp.nstr(s[j], 12) for j in range(min(3, len(rows[0])))])
