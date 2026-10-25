import sys
p = 2^31 - 1
for tag in sys.argv[1:]:
    cs = load('op_%s_modp.sobj' % tag)
    lc = cs[-1].leading_coefficient()
    Rq.<l> = QQ[]
    out = []
    for c in cs:
        out.append(Rq([rational_reconstruction(ZZ(a / lc), p) for a in c.list()]))
    # clear denominators
    den = lcm([x.denominator() for c in out for x in c.coefficients()])
    out = [den * c for c in out]
    g = gcd([ZZ(x) for c in out for x in c.coefficients()])
    out = [c / g for c in out]
    print('==', tag, 'max coefficient height', max(abs(ZZ(x)) for c in out for x in c.coefficients()))
    for m, c in enumerate(out):
        print('  D^%d :' % m, c.factor())
    save(out, 'op_%s_Q.sobj' % tag)
