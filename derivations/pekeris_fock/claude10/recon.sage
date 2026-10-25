# Rational reconstruction in lam (mod p) of the normalised telescoper coefficients from sample files.
import sys, glob
p = 2^31 - 1; k = GF(p); Rl.<l> = k[]
tag, r = sys.argv[1], int(sys.argv[2])
pts = []
for fn in sorted(glob.glob('s_%s_*.txt' % tag)):
    for line in open(fn):
        a = line.split()
        if a[1] != 'FAIL': pts.append((k(a[0]), [k(x) for x in a[1:]]))
n = len(pts); hold = 6
print(tag, 'samples', n)
fit, chk = pts[:n-hold], pts[n-hold:]
M = prod(l - x for x, _ in fit)
ops = []
for m in range(r + 1):
    f = Rl.lagrange_polynomial([(x, v[m]) for x, v in fit])
    found = None
    for tot in range(0, len(fit) - 1):
        for dd in range(0, tot + 1):
            try:
                num, den = f.rational_reconstruction(M, tot - dd, dd)
            except Exception:
                continue
            if den != 0 and all(num(x) == v[m] * den(x) for x, v in chk):
                found = (num, den); break
        if found: break
    print('c_%d:' % m, 'FAILED' if found is None else 'deg num %d, den %d' % (found[0].degree(), found[1].degree()))
    ops.append(found)
if all(o is not None for o in ops):
    L = lcm([o[1] for o in ops])
    cs = [o[0] * (L // o[1]) for o in ops]
    g = gcd(cs)
    cs = [c // g for c in cs]
    print('common-denominator operator: coefficient degrees', [c.degree() for c in cs])
    print('leading coefficient factors:', cs[r].factor())
    save(cs, 'op_%s_modp.sobj' % tag)
