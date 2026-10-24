"""Fit every nonzero real/imag entry of the symmetric panel-1 connection as a bivariate rational function of (h, Y2):
degrees from a univariate line probe (first prime), dense fits at 4 primes, CRT + rational reconstruction, validation at the 5th prime."""
import pickle, json, sys
from fractions import Fraction
sys.path.insert(0, '../claude6')
import ratrec as RR
from panel1_reconstruct import eval_parts, fit, monos, ratrec_int, PRIMES, set_prime

data = pickle.load(open('panel1_points.pkl', 'rb'))
keys = [(ax, fam) for ax in (0, 1) for fam in ('A', 'D')]
p0 = PRIMES[0]
# --- line probe for degrees
RR.P = p0
xs, line = [], []
lam = 0
while len(xs) < 40:
    lam += 1
    try:
        line.append(eval_parts((5 + lam) % p0, (7 + 2 * lam) % p0, p0)); xs.append(lam)
    except (ValueError, StopIteration, ZeroDivisionError):
        pass
entries = []
for k in keys:
    for part in (0, 1):
        for i in range(39):
            for j in range(39):
                ys = [v[k][part][i][j] for v in line]
                if not any(ys): continue
                r = RR.ratrec(xs, ys, holdout=6)
                assert r is not None, (k, part, i, j)
                entries.append((k, part, i, j, r[0], r[1]))
print('nonzero entry parts:', len(entries), 'max degrees', max(e[4] for e in entries), max(e[5] for e in entries), flush=True)
result = {}
for (k, part, i, j, dn, dd) in entries:
    fits = []
    for p in PRIMES[:4]:
        pts, vals = data[p]
        f = fit(pts, [v[k][part][i][j] for v in vals], dn, dd, p)
        assert f is not None, ('fit failed', k, part, i, j, dn, dd, p)
        fits.append(f)
    assert len({f[2] for f in fits}) == 1
    # CRT each coefficient
    mod = 1
    for p in PRIMES[:4]: mod *= p
    def crt(vals):
        x = 0
        for v, p in zip(vals, PRIMES[:4]):
            m = mod // p; x += v * m * pow(m, p - 2, p)
        return x % mod
    num = [ratrec_int(crt([f[0][t] for f in fits]), mod) for t in range(len(fits[0][0]))]
    den = [ratrec_int(crt([f[1][t] for f in fits]), mod) for t in range(len(fits[0][1]))]
    assert all(x is not None for x in num + den), ('reconstruction failed', k, part, i, j)
    result[f'{k[0]}{k[1]}|{part}|{i}|{j}'] = dict(dn=dn, dd=dd, num=[str(x) for x in num], den=[str(x) for x in den])
# validation at the 5th prime on its sample points
p = PRIMES[4]; pts, vals = data[p]; bad = 0; checked = 0
for key, r in result.items():
    kk, part, i, j = key.split('|'); k = (int(kk[0]), kk[1]); part, i, j = int(part), int(i), int(j)
    mn, md = monos(r['dn']), monos(r['dd'])
    for (h, y), v in list(zip(pts, vals))[:20]:
        N = sum(Fraction(c).numerator * pow(Fraction(c).denominator, p - 2, p) * pow(h, a, p) * pow(y, b, p) for c, (a, b) in zip(r['num'], mn)) % p
        D = sum(Fraction(c).numerator * pow(Fraction(c).denominator, p - 2, p) * pow(h, a, p) * pow(y, b, p) for c, (a, b) in zip(r['den'], md)) % p
        checked += 1
        if (N - v[k][part][i][j] * D) % p: bad += 1
print('validation at fresh prime', p, ':', checked, 'checks,', bad, 'failures')
json.dump(result, open('panel1_connection.json', 'w'))
