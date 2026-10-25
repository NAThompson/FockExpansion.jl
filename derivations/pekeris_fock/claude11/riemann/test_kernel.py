# Test candidate Riemann/flux kernels K(Y;X,eta) against the exact-mode identity
#   int_I K(Y;X,eta) h(eta) d eta = h(Y) h(X) / h(b),  b = 3pi/2, I = (pi/2, 3pi/2),
# for the exact c=32 modes h0 = sin, h+ = 1+sqrt2 cos, h- = 1-sqrt2 cos.
import mpmath as mp
mp.mp.dps = 25
pi = mp.pi; b = 3*pi/2; s2 = mp.sqrt(2)
modes = {'h0': lambda z: mp.sin(z), 'h+': lambda z: 1 + s2*mp.cos(z), 'h-': lambda z: 1 - s2*mp.cos(z)}
def z0(Y, X, eta):
    ab, bb = (b + Y)/2, (Y - b)/2
    a, be = (X + eta)/2, (eta - X)/2
    return mp.cos(ab)*mp.cos(a) + mp.cos(bb)*mp.cos(be)
def lhs(Y, X, h, N):
    f = lambda eta: (lambda z: N(X, eta) * (4*z**3 - 3*z) / mp.sqrt(z*z - 1) if z > 1 else mp.mpf(0))(z0(Y, X, eta)) * h(eta)
    # find support endpoints in eta where z0 = 1
    grid = [pi/2 + pi*k/400 for k in range(401)]
    vals = [z0(Y, X, e) - 1 for e in grid]
    pts = [pi/2]
    for i in range(400):
        if vals[i]*vals[i+1] < 0:
            pts.append(mp.findroot(lambda e: z0(Y, X, e) - 1, (grid[i], grid[i+1]), solver='bisect'))
    pts.append(3*pi/2)
    return mp.quad(f, pts)
Ns = {'1': lambda X, e: 1, 'sin a sin b': lambda X, e: mp.sin((X+e)/2)*mp.sin((e-X)/2), '-cosX cos eta': lambda X, e: -mp.cos(X)*mp.cos(e)}
for (X, Y) in [(mp.mpf('2.2'), mp.mpf('2.9')), (mp.mpf('3.5'), mp.mpf('2.0')), (mp.mpf('4.3'), mp.mpf('4.0'))]:
    for nn, N in Ns.items():
        r = [lhs(Y, X, h, N) / (h(Y)*h(X)/h(b)) for h in modes.values()]
        print(float(X), float(Y), nn, [mp.nstr(x, 10) for x in r])
