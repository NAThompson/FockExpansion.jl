from monodromy import *
from mpmath import mpf, sqrt, eig, nstr, matrix, eye
ops = [[mpf(eval(x.replace('/', '*1.0/'))) if False else mp.mpf(x) if '/' not in x else mp.mpf(x.split('/')[0]) / mp.mpf(x.split('/')[1]) for x in l.split()] for l in open('op_e2_A.txt')]
l0 = sqrt(mpf(8) / 9)
sing = [l0, -l0, mpf(1), mpf(-1)]
base = mpf(1) / 5
for pt in (l0, -l0, mpf(1)):
    T = monodromy_matrix(ops, base, pt, mpf(1) / 20, sing)
    E, _ = eig(T)
    from mpmath import svd
    s = svd(T - eye(5), compute_uv=False)
    print('around', nstr(pt, 6), 'eigs', [nstr(e, 8) for e in E], 'sv(T-1)', [nstr(x, 3) for x in s])
