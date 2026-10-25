import mpmath as mp
mp.mp.dps = 74
x = mp.mpf('-0.010280837917801415227952594791537657432618436882542490235847238933562546690')
s2 = mp.sqrt(2); pi = mp.pi; L = mp.log(1 + s2); l2 = mp.log(2); G = mp.catalan
w = mp.gamma(0.25)**2 / (2 * mp.sqrt(2 * pi))
li = mp.polylog(2, s2 - 1); z3 = mp.zeta(3)
w1 = [('1',1),('pi',pi),('L',L),('l2',l2)]
w2 = [('pi^2',pi**2),('piL',pi*L),('pil2',pi*l2),('L^2',L**2),('Ll2',L*l2),('l2^2',l2**2),('G',G),('Li2',li)]
def run(name, B, scales=((1,'x'),)):
    for sc, sn in scales:
        r = mp.pslq([x*sc] + [b for _, b in B], maxcoeff=10**5, maxsteps=10**7)
        if r is None: print(name, sn, 'none'); continue
        terms = [(c, n) for c, (n, _) in zip(r[1:], B) if c]
        print(name, sn, 'x-coef', r[0], 'max|c|', max(abs(c) for c in r), terms)
B1 = w1 + [('s2 '+n, s2*v) for n, v in w1]
B2 = B1 + w2 + [('s2 '+n, s2*v) for n, v in w2]
run('weight<=1', B1, ((1,'x'),(pi,'pi x'),(pi**2,'pi^2 x')))
run('weight<=2', B2, ((1,'x'),(pi,'pi x')))
lem = [('w', w), ('s2 w', s2*w), ('pi/w', pi/w), ('s2 pi/w', s2*pi/w), ('w^2', w**2), ('pi^2/w^2', pi**2/w**2), ('s2 w^2', s2*w**2), ('s2 pi^2/w^2', s2*pi**2/w**2)]
run('weight<=1 + lemniscate', B1 + lem, ((1,'x'),(pi,'pi x')))
