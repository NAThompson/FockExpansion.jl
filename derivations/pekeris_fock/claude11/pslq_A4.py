import mpmath as mp
mp.mp.dps = 55
x = mp.mpf('-0.01028083791780141522795259479153765743261843693754253802')
s2 = mp.sqrt(2); pi = mp.pi; L = mp.log(1 + s2); l2 = mp.log(2); G = mp.catalan
w = mp.gamma(0.25)**2 / (2 * mp.sqrt(2 * pi))
li = mp.polylog(2, s2 - 1)
wt = [('1',1),('pi',pi),('L',L),('l2',l2),('pi^2',pi**2),('piL',pi*L),('pil2',pi*l2),('L^2',L**2),('Ll2',L*l2),('l2^2',l2**2),('G',G),('Li2',li)]
B = [(n, v) for n, v in wt] + [('s2 ' + n, s2 * v) for n, v in wt]
for nm, BB in (('w<=2 (x and x*pi)', B), ('plus lemniscate', B + [('w', w), ('s2w', s2 * w), ('w^2', w**2), ('pi^2/w^2', pi**2 / w**2)])):
    for scale, sn in ((1, 'x'), (pi, 'pi*x')):
        r = mp.pslq([x * scale] + [b for _, b in BB], maxcoeff=2000, maxsteps=10**6)
        print(nm, sn, None if r is None else (r[0], [(c, n) for c, (n, _) in zip(r[1:], BB) if c]))
