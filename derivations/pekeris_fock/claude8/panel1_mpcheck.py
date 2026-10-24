"""Debug: the raw (unsymmetrized) connection formula of panel1_fp in mpmath complex arithmetic vs finite differences of actual integrals."""
import mpmath as m
m.mp.dps = 40
def jmul(x, y): return [x[0]*y[0], x[0]*y[1] + x[1]*y[0], x[0]*y[2] + x[1]*y[1] + x[2]*y[0]]
def jadd(x, y): return [a + b for a, b in zip(x, y)]
def jsc(x, c): return [a * c for a in x]
def jinv(x):
    v0 = 1/x[0]; v1 = -x[1]*v0*v0; v2 = -(x[1]*v1 + x[2]*v0)*v0; return [v0, v1, v2]
def jc(c): return [m.mpc(c), 0, 0]
I = 1j
def geom(h, Y2):
    b = (h*h - 1)/(h*h + 1); c4, c2, c0 = 1 - b, 4*Y2 - 2, 1 + b
    disc = m.sqrt(c2*c2 - 4*c4*c0 + 0j); z = m.sqrt((-c2 + disc)/(2*c4))
    return b, (c4, c2, c0), [z, -z, h/z, -h/z, I, -I, I*h, -I*h]
def connection(h, Y2, axis, family):
    b, (c4, c2, c0), z = geom(h, Y2)
    if axis == 0:
        bx = 4*h/(h*h + 1)**2; T_x = lambda q: bx*(1 - q**4); hx = 1; l0 = -bx/(2*(1 + b))
    else:
        bx = 0; T_x = lambda q: 4*q*q; hx = 0; l0 = 0
    Tq = lambda q: 4*c4*q**3 + 2*c2*q
    v = [-T_x(q)/Tq(q) for q in z[:4]] + [0, 0, I*hx, -I*hx]
    tj = [0, 1, 0] if family == 'A' else [0, 0, 0]; sj = [0, 1, 0] if family == 'D' else [0, 0, 0]
    i2 = 1/(2*I)
    a = [jc(-0.5) for _ in range(4)] + [jsc(tj, i2), jsc(tj, -i2), jsc(sj, i2), jsc(sj, -i2)]
    kap = jc(l0)
    for ai, zi, vi in zip(a, z, v): kap = jadd(kap, jsc(ai, -vi/zi))
    n = 39; C = [[0]*n for _ in range(n)]
    def addjet(rs, cs, coef):
        for jr in range(3):
            for jc_ in range(jr + 1): C[rs*3 + jr][cs*3 + jc_] += coef[jr - jc_]
    J, K, Js, E, F0 = 0, 1, 10, 11, 12; Ji = lambda i: 2 + i
    addjet(J, J, kap)
    for i in range(8): addjet(J, Ji(i), jsc(a[i], -v[i]))
    addjet(K, K, kap)
    for i in range(8): addjet(K, J, jsc(a[i], -v[i])); addjet(K, Ji(i), jsc(a[i], -v[i]*z[i]))
    for i in range(8):
        addjet(Ji(i), Ji(i), kap); w = jc(-1)
        addjet(Ji(i), E, jsc(w, v[i]/(1 - z[i]))); addjet(Ji(i), F0, jsc(w, v[i]/z[i]))
        for j in range(8):
            if j == i: continue
            d = 1/(z[i] - z[j]); coef = jadd(jsc(jmul(w, a[j]), -v[i]), jsc(a[j], -v[j]))
            addjet(Ji(i), Ji(i), jsc(coef, d)); addjet(Ji(i), Ji(j), jsc(coef, -d))
    addjet(Js, Js, kap)
    for i in range(8):
        f = jsc(a[i], v[i]/z[i]); addjet(Js, Js, f); addjet(Js, Ji(i), jsc(f, -1))
    ep = kap
    for i in range(8): ep = jadd(ep, jsc(a[i], -v[i]/(1 - z[i])))
    addjet(E, E, ep); addjet(F0, F0, jc(l0))
    return C
def states(h, Y2, fam):
    b, (c4, c2, c0), zs = geom(h, Y2)
    W = lambda q: m.sqrt(c4*q**4 + c2*q*q + c0)
    ang = (lambda q: m.atan(q)) if fam == 'A' else (lambda q: m.atan(q/h))
    vec = [0]*39
    for nn in range(3):
        f = lambda q: ang(q)**nn/(m.factorial(nn)*W(q))
        base = 1/m.sqrt(c0) if nn == 0 else 0
        vals = [m.quad(f, [0, 1]), m.quad(lambda q: q*f(q), [0, 1])] + [m.quad(lambda q, zi=zi: f(q)/(q - zi), [0, 1]) for zi in zs] + \
               [m.quad(lambda q: (f(q) - base)/q, [0, 1]), f(1), base]
        for st, val in enumerate(vals): vec[st*3 + nn] = val
    return vec
h0, Y0 = m.mpf('1.7'), m.mpf('0.21')
for ax in (0, 1):
    for fam in ('A', 'D'):
        d = m.mpf('1e-12')
        Xp = states(h0 + (d if ax == 0 else 0), Y0 + (d if ax == 1 else 0), fam)
        Xm = states(h0 - (d if ax == 0 else 0), Y0 - (d if ax == 1 else 0), fam)
        X0 = states(h0, Y0, fam); C = connection(h0, Y0, ax, fam)
        errs = []
        for r in range(39):
            fd = (Xp[r] - Xm[r])/(2*d); an = sum(C[r][c]*X0[c] for c in range(39))
            errs.append((abs(fd - an)/(1 + abs(fd)), r))
        print(ax, fam, 'worst', m.nstr(max(errs)[0], 3), 'at state', max(errs)[1]//3, 'jet', max(errs)[1] % 3)
print('--- per-state errors, axis 1, family A')
d = m.mpf('1e-12'); fam = 'A'
Xp = states(h0, Y0 + d, fam); Xm = states(h0, Y0 - d, fam); X0 = states(h0, Y0, fam); C = connection(h0, Y0, 1, fam)
b, cs, zs = geom(h0, Y0); print('roots', [m.nstr(z, 5) for z in zs[:4]])
for r in range(0, 39, 3):
    fd = (Xp[r] - Xm[r])/(2*d); an = sum(C[r][c]*X0[c] for c in range(39))
    print('state', r//3, 'fd', m.nstr(fd, 8), 'an', m.nstr(an, 8))
