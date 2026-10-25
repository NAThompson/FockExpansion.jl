import numpy as np, itertools
def polys(a, x1, x2):
    q = np.polynomial.Polynomial([0, 1])
    D = 1 + q**4; P1, P2 = q**4 - 1, 2 * q**2
    Px = P1 * x1 + P2 * x2
    if a == 'e2': ax, aP = x2, P2
    elif a == 'e1': ax, aP = x1, P1
    else: ax, aP = -x1, -P1
    A = 1 - ax**2; B = aP - ax * Px; C = D**2 - Px**2
    return A, B, C
for a in ('e2', 'e1'):
    for x1, x2 in ((0.3, 0.2), (-0.4, 0.5), (0.1, -0.6), (0.6, 0.6), (-0.2, -0.3)):
        A, B, C = polys(a, x1, x2)
        disc = B**2 - A * C
        r = disc.roots(); rr = sorted(r[abs(r.imag) < 1e-9].real)
        print(a, (x1, x2), 'real roots of B^2-AC:', np.round(rr, 4), ' C real roots:', np.round(sorted(C.roots()[abs(C.roots().imag)<1e-9].real),4))
