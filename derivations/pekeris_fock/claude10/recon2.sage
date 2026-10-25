# Multi-prime reconstruction: per prime, canonical c_m/c_r = num/den (den monic, reduced); CRT + rational reconstruction over Q.
import sys, glob
tag, r = sys.argv[1], int(sys.argv[2])
primes = {2^31 - 1: 's_%s_*.txt' % tag, 2147483629: 't_%s_2147483629_*.txt' % tag, 2147483587: 't_%s_2147483587_*.txt' % tag}
def canon(p, pattern):
    k = GF(p); Rl.<l> = k[]
    pts = []
    for fn in sorted(glob.glob(pattern)):
        for line in open(fn):
            a = line.split()
            if len(a) > 1 and a[1] != 'FAIL': pts.append((k(a[0]), [k(x) for x in a[1:]]))
    hold = 3; fit, chk = pts[:-hold], pts[-hold:]
    M = prod(l - x for x, _ in fit)
    res = []
    for m in range(r + 1):
        f = Rl.lagrange_polynomial([(x, v[m]) for x, v in fit])
        found = None
        for tot in range(0, len(fit) - 1):
            for dd in range(0, tot + 1):
                try: num, den = f.rational_reconstruction(M, tot - dd, dd)
                except Exception: continue
                if den != 0 and all(num(x) == v[m] * den(x) for x, v in chk):
                    g = gcd(num, den); num, den = num // g, den // g
                    lc = den.leading_coefficient(); found = (num / lc, den / lc); break
            if found: break
        res.append(found)
    return len(pts), res
data = {}
for p, pat in primes.items():
    n, res = canon(p, pat)
    print('prime', p, 'samples', n, 'degrees', [(a.degree(), b.degree()) if a is not None and (a, b) != (None, None) else None for a, b in [x if x else (None, None) for x in res]])
    data[p] = res
Rq.<l> = QQ[]
ps = list(data)
def lift(coeff_lists):
    out = []
    for cs in zip(*coeff_lists):
        a = crt([ZZ(c) for c in cs], ps)
        out.append(rational_reconstruction(a, prod(ps)))
    return Rq(out)
ops = []
for m in range(r + 1):
    nums = [data[p][m][0].list() for p in ps]; dens = [data[p][m][1].list() for p in ps]
    L = max(len(x) for x in nums); nums = [x + [0] * (L - len(x)) for x in nums]
    L = max(len(x) for x in dens); dens = [x + [0] * (L - len(x)) for x in dens]
    ops.append(lift(nums) / lift(dens))
den = lcm([c.denominator() for c in ops])
cs = [Rq(c * den) for c in ops]
g = gcd([x for c in cs for x in c.coefficients()]); cs = [c / g for c in cs]
cz = lcm([x.denominator() for c in cs for x in c.coefficients()]); cs = [c * cz for c in cs]
print('operator over Q (coefficient of D^m):')
for m, c in enumerate(cs): print('  D^%d:' % m, c.factor())
save(cs, 'op_%s_Q.sobj' % tag)
