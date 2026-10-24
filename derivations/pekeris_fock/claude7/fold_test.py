"""Is a single elliptic panel moment N0(b,y), N2(b,y) analytic across b = 0 (where the endpoint u = min(1,h) switches)?
Fit a polynomial in b on b>0 and extrapolate to b<0 (and vice versa); compare with direct values."""
import sys
sys.path.insert(0, '../codex/phase3')
import mpmath as m
from panel_connection import direct_moments
m.mp.dps = 40
y = m.mpf('-0.6')
bs_pos = [m.mpf(k) / 400 for k in range(1, 13)]
bs_neg = [-b for b in bs_pos]
for idx, name in ((0, 'N0'), (1, 'N2')):
    fp = [direct_moments(b, y)[idx] for b in bs_pos]
    fn = [direct_moments(b, y)[idx] for b in bs_neg]
    # polynomial fit (degree 9) on positive side, evaluate at negative points

    def fit_eval(xs, ys, xe):
        V = m.matrix([[x**k for k in range(10)] for x in xs]); c = m.lu_solve(V.T * V, V.T * m.matrix(ys))
        return [sum(c[k] * x**k for k in range(10)) for x in xe]
    ext = fit_eval(bs_pos, fp, bs_neg[:4])
    print(name, 'extrapolation error across b=0:', [m.nstr(e - f, 3) for e, f in zip(ext, fn[:4])])
    fd = [(direct_moments(m.mpf('1e-8'), y)[idx] - direct_moments(-m.mpf('1e-8'), y)[idx]) / m.mpf('2e-8')]
    l1 = (direct_moments(m.mpf('2e-8'), y)[idx] - direct_moments(m.mpf('1e-8'), y)[idx]) / m.mpf('1e-8')
    r1 = (direct_moments(-m.mpf('1e-8'), y)[idx] - direct_moments(-m.mpf('2e-8'), y)[idx]) / m.mpf('1e-8')
    print(name, 'one-sided slopes at b=0:', m.nstr(l1, 12), m.nstr(r1, 12))
