"""Admissibility probe for single panel moments as functions of the angular coordinates (x,y): fit tau-expansions
(tau^k, tau^k log tau, tau^(k+1/2)) along normal paths into faces/corners; real non-analytic terms persist when K grows."""
import sys
sys.path.insert(0, '../codex/phase3')
import mpmath as m
from panel_connection import direct_moments
m.mp.dps = 90
def maps(x, y):
    r1, r2, c = (1 + x) / 2, (1 - x) / 2, y
    R = m.sqrt(r1**2 + r2**2); v = (R**2 - c**2) / R**2
    return [((r1**2 - r2**2) / R**2, -c / (m.sqrt(2) * R)), (-(r1**2 - r2**2) / R**2, -c / (m.sqrt(2) * R)), (-v, -r1 / R), (-v, -r2 / R)]
paths = [('w-face y=1, x=0.3', lambda t: (m.mpf('0.3'), 1 - t)), ('u-face y=x, x=0.4', lambda t: (m.mpf('0.4'), m.mpf('0.4') + t)),
         ('e-e corner', lambda t: (t / 2, t)), ('nucleus corner', lambda t: (1 - 2 * t, 1 - t))]
K = int(sys.argv[1]) if len(sys.argv) > 1 else 5
for name, path in paths:
    ncol = 3 * K + 2
    taus = [m.mpf(2)**(-j) for j in range(14, 14 + ncol)]
    for pj in range(4):
        for idx in (0, 1):
            ys = [direct_moments(*maps(*path(t))[pj])[idx] for t in taus]
            cols = []
            labs = []
            for k in range(K + 1):
                cols.append(lambda t, k=k: t**k); labs.append(f't^{k}')
                if k >= 1: cols.append(lambda t, k=k: t**k * m.log(t)); labs.append(f't^{k}log')
                cols.append(lambda t, k=k: t**k * m.sqrt(t)); labs.append(f't^{k}.5')
            V = m.matrix([[c(t) for c in cols] for t in taus]); cf = m.lu_solve(V, m.matrix(ys))
            nonan = [(l, m.nstr(cf[i], 3)) for i, l in enumerate(labs) if ('log' in l or '.5' in l) and abs(cf[i]) > 1e-8]
            print(name, f'panel{pj+1}', 'N0' if idx == 0 else 'N2', 'non-analytic terms > 1e-8:', nonan[:3], flush=True)
