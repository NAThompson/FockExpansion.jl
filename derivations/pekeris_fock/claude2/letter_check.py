import mpmath as mp
mp.mp.dps = 30
for (a, b, c) in [(0.7, 0.3, 0.5), (0.4, 0.6, 0.3), (0.9, 0.2, 0.95), (0.5, 0.5, 0.1)]:
    a, b, c = map(mp.mpf, (a, b, c))
    R = mp.sqrt(a*a + b*b); P = mp.sqrt(2*R*R - c*c); sh, d = a + b, a - b
    L1, L2, L3, L4, L5 = sh*P + d*c, sh*c - d*P, sh*P - d*c, sh*c + d*P, c*P + sh*d
    A = (sh*sh - c*c) * L4 / (L5 * L3)
    B = (c*c - d*d) * L1 / (L5 * L2)
    print(mp.nstr(A, 12), mp.nstr(B, 12), ' L1L3/(..)', mp.nstr(L1*L3, 8), mp.nstr(L2*L4, 8), mp.nstr(L5*(c*P - sh*d), 8))
